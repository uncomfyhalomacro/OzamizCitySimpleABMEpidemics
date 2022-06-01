@inline function migrate!(agent::AbstractAgent, model::AgentBasedModel)
    pid = agent.pos
    d = DiscreteNonParametric(1:(model.number_of_barangays), model.migration_rates[pid, :])
    pos_rng = rand(model.rng, d)
    if model ≠ pid
        move_agent!(agent, pos_rng, model)
    end
    return model
end

# NOTE: reproduction
@inline function transmit!(agent::AbstractAgent, model::AgentBasedModel;
                           from_model_step = false)
    agent.status == :S && return nothing
    rate = if agent.days_infected < model.detection_time
        model.β_und[agent.pos]
    else
        model.β_det[agent.pos]
    end

    d = Poisson(rate)
    n = rand(model.rng, d)
    n == 0 && return nothing
    individual_reproduction_number = 0
    for contactID in ids_in_position(agent, model)
        contact = model[contactID]
        waning_protection_per_day = if contact.num_vac_shots == 3
            contact.days_vaccinated * -0.000095
        elseif contact.days_vaccinated == 2
            contact.days_vaccinated * -0.0008
        elseif contact.days_vaccinated == 1
            contact.days_vaccinated * -0.002
        else
            0.0
        end

        protection_from_vaccine = contact.vaccination_status == (:V) ?
                                  ifelse(contact.num_vac_shots == 1,
                                         rand(0.002:0.0006:0.08),
                                         ifelse(contact.num_vac_shots == 2,
                                                rand(0.08:0.002:0.12),
                                                # NOTE: number of vaccine shots == 3
                                                rand(0.2:0.04:0.65))) +
                                  waning_protection_per_day :
                                  0.0

        if ((contact.status == :S &&
             (rand(model.rng) + protection_from_vaccine) ≤ model.infection_probability) ||
            (contact.status == :R &&
             (rand(model.rng) + protection_from_vaccine) ≤ model.reinfection_probability))
            contact.status = :I
            individual_reproduction_number += 1
        elseif contact.date_infected == model.days_passed
            contact.status = :I
            individual_reproduction_number += 1
        end

        # NOTE: Record the agent's reproduction number at this point
        agent.individual_reproduction_number += individual_reproduction_number
        n -= 1
        n == 0 && return nothing
    end
end

# @inline function recover_or_die!(agent::AbstractAgent, model::AgentBasedModel)
#     if agent.days_infected ≥ model.infection_period
#         if rand(model.rng) ≤
#            model.death_rate + ifelse(agent.vaccination_status == :V,
#                   ifelse(agent.num_vac_shots == 3, -0.0025, -0.0015),
#                   0.0)
#             Agents.kill_agent!(agent, model)
#         else
#             agent.status = :R
#             agent.days_infected = 0
#         end
#     end
# end

@inline function update!(agent::AbstractAgent, model::AgentBasedModel)
    if agent.status === :I
        agent.days_infected += 1
        agent.days_recovered = 0
    end
    if agent.status === :R
        agent.days_recovered += 1
        agent.days_infected = 0
    end
    if agent.vaccination_status == :V
        agent.days_vaccinated += 1
    end
    return model
end

@inline function agent_step!(agent::AbstractAgent, model::AgentBasedModel)
    migrate!(agent, model)
    transmit!(agent, model)
    update!(agent, model)
end
