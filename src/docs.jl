
#=
    Docstrings for the step functions:
    - `migrate!`
    - `transmit!`
    - `recover_or_die!`
    - `update!`
=#

"""
    migrate!
    migrate!(agent::AbstractAgent, model::AgentBasedModel)

A step function that simulates migration of individuals from place to place.
The real implementation takes two arguments:
    - agent → `AbstractAgent`
    - model → `AgentBasedModel` or `ABM`
"""
function migrate! end

"""
    transmit!
    transmit!(agent::AbstractAgent, model::AgentBasedModel)

A step function that simulates infection/reinfections.
Close individuals get infected through `ids_in_position`

The real implementation takes two arguments:
    - agent → `AbstractAgent`
    - model → `AgentBasedModel`
"""
function transmit! end

"""
    recover_or_die!
    recover_or_die!(agent::AbstractAgent, model::AgentBasedModel)

A step function that simulates recoveries and deaths

The real implementation takes two arguments:
    - agent → `AbstractAgent`
    - model → `AgentBasedModel`
"""
function recover_or_die! end

"""
    update!
    update!(agent::AbstractAgent, model::AgentBasedModel)

Updates the agents statuses e.g. `days_infected`, `days_recovered`, `days_vaccinated`

The real implementation takes two arguments:
    - agent → `AbstractAgent`
    - model → `AgentBasedModel`
"""
function update! end

# Docstrings for the model step functions:

"""
	model_step!(model::ABM)

The model step function. This updates the whole model with each step.
Each step corresponds as **one day** in the model.
The `model_step!` function takes other model step functions and
other utility functions that are relevant when updating the whole
model e.g. number of susceptible, infected, recovered, dead, and vaccinated

"""
function model_step! end

"""
	update_numbers!(model::ABM)
Updates the total number of susceptible, infected, recovered, vaccinated, and dead
in the model.
"""
function update_numbers! end

"""
	vaccinate!(model::ABM)
Starts the vaccination for both the *not fully* vaccinated
and the unvaccinated. This function does this by holding two
other functions `vaccinate_the_vaccinated!` and `vaccinate_the_unvaccinated!`
Both functions are self-explanatory.
NOTE: there is an increased probability of transmission during vaccinations
as the happen in crowded places
"""
function vaccinate! end

"""
	update_vaccination_hesitancy!(model::ABM)

Update the vaccination hesitancy probability in each barangay.
Both "normal" and "lowered" are updated.
"""
function update_vaccination_hesitancy! end

"""
	update_R₀!(model::ABM)

Update the mean and median R-nought of the current model step.
"""
function update_R₀! end
