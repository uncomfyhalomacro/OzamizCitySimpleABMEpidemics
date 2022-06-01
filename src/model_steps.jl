@inline function model_step!(model::AgentBasedModel)
    recover_or_die!(model)
    if model.days_passed ≥ model.starting_date + Dates.Day(313)
        vaccinate!(model)
        update_infection_probability!(model)
    end
    vaccine_information_bias!(model)
    update_vaccination_hesitancy!(model)
    update_migration_rate!(model)
    update_β_infected!(model)
    update_numbers!(model)
    if model.days_passed ≥ model.starting_date + Dates.Day(23)
        update_infection_probability!(model)
    end
    model.days_passed += Day(1)
    return model
end

function update_numbers!(model::ABM)
    model.susceptible = total_susceptible(model)
    model.infected = total_infected(model)
    model.recovered = total_recovered(model)
    model.partially_vaccinated = total_partially_vaccinated(model)
    model.fully_vaccinated_without_booster = total_fully_vaccinated_without_booster(model)
    model.fully_vaccinated_with_booster = total_fully_vaccinated_with_booster(model)
    model.dead = total_dead(model)
    return model
end

@inline function vaccinate!(model::ABM)
    # NOTE: Vaccination happens in a crowded place, thus, there is a high possibility of transmission
    vaccinate_the_unvaccinated!(model)
    vaccinate_the_vaccinated!(model)
    return model
end

function vaccinate_the_unvaccinated!(model::ABM)
    persons = filter(a -> a.vaccination_status == :Z, collect(allagents(model)))
    length(persons) == 0 && return model
    for to_vaccinate in 1:min(1500, length(persons))
        person = persons[to_vaccinate]
        agent = model[person.id]
        rate = model.vaccination_hesitancy[agent.pos]
        d = Poisson(rate)
        n = rand(model.rng, d)

        if n ≠ 0
            if rand(model.rng) ≤ model.vaccination_probability + model.vaccine_informed_bias
                agent.vaccination_status = :V
                agent.num_vac_shots += 1
            end
        end
        transmit!(agent, model; from_model_step = true)
        model[person.id] = agent
    end
    return model
end

function vaccinate_the_vaccinated!(model::ABM)
    vaccinated_persons = filter(a -> a.vaccination_status == :V && a.num_vac_shots < 3,
                                collect(allagents(model)))

    # NOTE: this is when all people are fully vaccinated if num_vac_shots == 3, hence the filter function having `a.num_vac_shots < 3`
    length(vaccinated_persons) == 0 && return model

    for to_vaccinate in 1:min(900, length(vaccinated_persons))
        person = vaccinated_persons[to_vaccinate]
        agent = model[person.id]
        rate = model.vaccination_hesitancy_lowered[agent.pos]
        d = Poisson(rate)
        n = rand(model.rng, d)
        if n ≠ 0
            # for second vaccine shot
            if (agent.num_vac_shots == 1 &&
                agent.days_vaccinated ≥ 90) &&
               rand(model.rng) ≤ model.vaccination_probability + model.vaccine_informed_bias
                agent.num_vac_shots += 1
                agent.days_vaccinated = 0

                # for third vaccine shot (booster)
            elseif (agent.num_vac_shots == 2 &&
                    agent.days_vaccinated ≥ 90) &&
                   rand(model.rng) ≤
                   model.vaccination_probability + model.vaccine_informed_bias
                agent.num_vac_shots += 1
                agent.days_vaccinated = 0
            else
                continue
            end
        end
        transmit!(agent, model; from_model_step = true)
        model[person.id] = agent
    end
    return model
end

function vaccinated_bias()
    Random.seed!(RandomDevice())
    return rand(0.0:0.0002:0.01)
end

function vaccine_information_bias!(model::ABM)
    Random.seed!(RandomDevice())
    model.vaccine_informed_bias = rand(0.02:0.002:0.6)
    return model
end

function update_vaccination_hesitancy!(model::ABM)
    Random.seed!(RandomDevice())
    vaccination_hesitancy_lowered = rand(0.001:0.02:0.6, model.number_of_barangays)
    vaccination_hesitancy = vaccination_hesitancy_lowered ./ 2.0

    model.vaccination_hesitancy = vaccination_hesitancy
    model.vaccination_hesitancy_lowered = vaccination_hesitancy_lowered
    return model
end

function update_migration_rate!(model::ABM)
    Random.seed!(RandomDevice())
    max_travel_rate = 0.02
    migration_rates = zeros(model.number_of_barangays, model.number_of_barangays)
    for b in 1:(model.number_of_barangays)
        for b2 in 1:(model.number_of_barangays)
            migration_rates[b, b2] = (model.Ns[b] + model.Ns[b2]) / model.Ns[b]
        end
    end

    max_migration_rate = maximum(migration_rates)
    migration_rates = (migration_rates .* max_travel_rate) ./ max_migration_rate
    migration_rates[diagind(migration_rates)] .= 1.0
    migration_rates_sum = sum(migration_rates; dims = 2)
    for c in 1:(model.number_of_barangays)
        migration_rates[c, :] ./= migration_rates_sum[c]
    end
    model.migration_rates = migration_rates
    return model
end

function update_β_infected!(model::ABM)
    Random.seed!(RandomDevice())
    model.β_und = rand(0.3:0.02:0.6, model.number_of_barangays)
    model.β_det = model.β_und ./ 10.0
    return model
end

function update_R₀!(model::ABM)
    #= NOTE: naive calculation of reproduction number.
    		The basic reproduction number is the number of secondary cases or infections **one** individual can infect
        within their infection period from a given point in time and space during an epidemic.
    		Here we find the mean and median of all individual reproduction number at each step
    =#
    all_infected = filter(a -> a.status == (:I) && a.days_infected ≥ model.infection_period,
                          collect(allagents(model)))
    model.mean_R₀ = !isempty(all_infected) ?
                    mean([agent.individual_reproduction_number for agent in all_infected]) :
                    model.mean_R₀

    model.median_R₀ = !isempty(all_infected) ?
                      median([agent.individual_reproduction_number
                              for agent in all_infected]) : model.median_R₀
    return model
end

@inline function recover_or_die!(model::AgentBasedModel)
    all_infected = filter(a -> a.status == (:I) && a.days_infected ≥ model.infection_period,
                          collect(allagents(model)))
    update_R₀!(model)
    for agent in all_infected
        if rand(model.rng) ≤
           model.death_rate + ifelse(agent.vaccination_status == :V,
                  ifelse(agent.num_vac_shots == 3, 0.025, 0.015),
                  0.0)
            Agents.kill_agent!(agent, model)
        else
            agent.status = :R
            agent.days_infected = 0
        end
    end
    return model
end

function update_infection_probability!(model::ABM)
    Random.seed!(RandomDevice())
    model.infection_probability = rand(0.08:0.002:0.24)
    return model
end

function update_vaccination_probability!(model::ABM)
    Random.seed!(RandomDevice())
    model.vaccination_probability = rand(0.02:0.002:0.09788)
end
