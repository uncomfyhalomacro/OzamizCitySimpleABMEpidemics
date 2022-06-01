module OzamizCitySimpleABMEpidemics

using Agents
using CSV
using DataFrames
using Distributions: Poisson, DiscreteNonParametric
using DrWatson: @dict
using GLMakie
using Graphs
using LinearAlgebra: diagind
using Random
using Statistics
using Dates
using InteractiveDynamics

# Initialization
export Person
export COVID19DATA
export initialize_parameters
export initialize_model

# Agent step
export agent_step!
export migrate!
export transmit!
export recover_or_die!
export update!

# Agent utils
export generate_agents!
export generate_agents_from_csv!
export generate_random_infected!
export generate_agents_from_csv!

# Model step
export model_step!
export update_numbers!
export vaccinate!

# Model utils
export vaccine_information_bias!
export vaccinated_bias

# General utils
export total_susceptible
export total_infected
export total_recovered
export total_partially_vaccinated
export total_dead

include("./docs.jl")
include("./data.jl")
include("./utils.jl")
include("./model_steps.jl")
include("./agent_steps.jl")

"""
    mutable struct Person <: AbstractAgent

An agent representing a person within a `GraphSpace`.
It has the following parameters:
  `id` → self-explanatory. each agent has a unique id
  `pos` → position in the `GraphSpace`
  `days_infected`
  `days_vaccinated`
  `days_recovered`
  `num_vac_shots`
  `age`
  `status`
  `vaccination_status`
  `address`
  `date_infected`
  `from_csv_data`

The GraphSpace represents the all the barangays in the city

"""
mutable struct Person <: AbstractAgent
    id::Int
    pos::Int
    days_infected::Int64
    days_vaccinated::Int64
    days_recovered::Int64
    num_vac_shots::Int64
    age::Float64
    status::Symbol
    vaccination_status::Symbol
    address::Symbol
    date_infected::Dates.Date
    from_csv_data::Bool
    individual_reproduction_number::Int64
end

@inline function initialize_parameters(;
                                       number_of_barangays,
                                       max_travel_rate = 0.02,
                                       infection_period = 30,
                                       infection_probability = 0.08,
                                       reinfection_probability = 0.05,
                                       vaccination_probability = 0.12,
                                       death_rate = 0.012,
                                       detection_time = 14,
                                       Is = [zeros(Int, number_of_barangays - 1)..., 1],
                                       undetected_probability_range,
                                       vaccination_hesitancy_probability_range)
    Random.seed!(RandomDevice())
    Ns = [OZAMIZ_BARANGAY_POPULATION...]
    number_of_barangays = Int64(length(Ns))
    β_und = rand(undetected_probability_range, number_of_barangays)
    β_det = β_und ./ 10.0
    vaccination_hesitancy_lowered = rand(vaccination_hesitancy_probability_range,
                                         number_of_barangays)

    vaccination_hesitancy = vaccination_hesitancy_lowered ./ 2.0
    Random.seed!(RandomDevice())
    migration_rates = zeros(number_of_barangays, number_of_barangays)
    for b in 1:number_of_barangays
        for b2 in 1:number_of_barangays
            migration_rates[b, b2] = (Ns[b] + Ns[b2]) / Ns[b]
        end
    end
    max_migration_rate = maximum(migration_rates)
    migration_rates = (migration_rates .* max_travel_rate) ./ max_migration_rate
    migration_rates[diagind(migration_rates)] .= 1.0
    days_passed = Dates.Date("2020-04-22", "yyyy-mm-dd")
    params = @dict(Ns,
                   migration_rates,
                   β_und,
                   β_det,
                   vaccination_hesitancy,
                   vaccination_hesitancy_lowered,
                   vaccination_probability,
                   infection_period,
                   infection_probability,
                   reinfection_probability,
                   detection_time,
                   death_rate,
                   Is,
                   days_passed)

    return params
end

@inline function initialize_model(;
                                  Ns,
                                  migration_rates,
                                  β_und,
                                  β_det,
                                  vaccination_probability,
                                  vaccination_hesitancy,
                                  vaccination_hesitancy_lowered,
                                  infection_period,
                                  infection_probability,
                                  reinfection_probability,
                                  death_rate,
                                  detection_time,
                                  Is = [zeros(Int, Ns - 1)..., 1],
                                  days_passed)
    rng = RandomDevice()
    @assert length(Ns) ==
            length(Is) ==
            size(migration_rates, 1) ==
            length(β_und) ==
            length(β_det) ==
            length(vaccination_hesitancy) ==
            length(vaccination_hesitancy_lowered)

    @assert size(migration_rates, 1) == size(migration_rates, 2)

    starting_date = days_passed
    vaccine_informed_bias = 0
    number_of_barangays = Int64(length(Ns))
    migration_rates_sum = sum(migration_rates; dims = 2)
    for c in 1:number_of_barangays
        migration_rates[c, :] ./= migration_rates_sum[c]
    end

    mean_R₀ = 0.0
    median_R₀ = 0.0

    susceptible, infected, recovered, partially_vaccinated, fully_vaccinated_without_booster, fully_vaccinated_with_booster, dead = fill(Int64(0),
                                                                                                                                         7)
    properties = @dict(Ns,
                       Is,
                       β_und,
                       β_det,
                       vaccination_hesitancy,
                       vaccination_hesitancy_lowered,
                       vaccine_informed_bias,
                       vaccination_probability,
                       migration_rates,
                       infection_period,
                       infection_probability,
                       reinfection_probability,
                       detection_time,
                       number_of_barangays,
                       death_rate,
                       days_passed,
                       starting_date,
                       susceptible,
                       infected,
                       recovered,
                       partially_vaccinated,
                       fully_vaccinated_without_booster,
                       fully_vaccinated_with_booster,
                       dead,
                       mean_R₀,
                       median_R₀,)

    space = GraphSpace(complete_digraph(number_of_barangays))
    @info space
    model = ABM(Person, space; properties, rng)
    model = generate_agents_from_csv!(model)

    model = generate_agents!(model, number_of_barangays)
    ## Comment this out if you want
    # model = generate_random_infected!(model, number_of_barangays)
    ##
    model.susceptible = total_susceptible(model)
    model.infected = total_infected(model)
    model.partially_vaccinated = total_partially_vaccinated(model)
    model.dead = total_dead(model)
    model.vaccine_informed_bias = vaccine_information_bias!(model)
    @info model
    return model
end

@inline function generate_agents_from_csv!(model::AgentBasedModel)
    for row in eachrow(COVID19DATA)
        address = Symbol(row.ADDRESS)
        if address ∈ keys(OZAMIZ_BARANGAY_POPULATION)
            age = row.AGE
            id = findfirst(i -> ==(i, address), keys(OZAMIZ_BARANGAY_POPULATION))
            date_infected = row.DATES
            status = date_infected == model.starting_date ? (:I) : (:S)
            add_agent!(id,
                       model,
                       fill(Int64(0), 4)...,
                       age,
                       status,
                       :Z,
                       address,
                       date_infected,
                       true,
                       0)
        end
    end
    return model
end

@inline function generate_agents!(model::AgentBasedModel, total_barangays::Int64)
    OZAMIZ_BARANGAY_KEYNAMES = keys(OZAMIZ_BARANGAY_POPULATION)
    dummy_date = model.starting_date - Day(1)
    for barangay_id in 1:total_barangays
        negate_counter = count(c -> OZAMIZ_BARANGAY_KEYNAMES[barangay_id] == Symbol(c),
                               COVID19DATA.ADDRESS)

        for _ in 1:(OZAMIZ_BARANGAY_POPULATION[barangay_id] - negate_counter)
            add_agent!(barangay_id,
                       model,
                       fill(Int64(0), 4)...,
                       rand(1.0:100.0),
                       :S,
                       :Z,
                       OZAMIZ_BARANGAY_KEYNAMES[barangay_id],
                       dummy_date,
                       false,
                       0)
        end
    end
    return model
end

@inline function generate_random_infected!(model::ABM, total_barangays::Int64)
    dummy_date = model.starting_date - Day(1)
    for barangay in 1:(total_barangays - 34)
        inds = ids_in_position(barangay, model)
        for n in 1:OZAMIZ_BARANGAY_POPULATION[barangay]
            agent = model[inds[n]]
            if model.starting_date > dummy_date
                if rand(model.rng) ≤ model.infection_probability
                    agent.status = (:I)
                    agent.days_infected = 1
                end
            else
                continue
            end
        end
    end
    return model
end

end
