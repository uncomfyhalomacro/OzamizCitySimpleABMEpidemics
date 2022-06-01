# An Agent-based model: Simulating available infection data from Ozamiz City

## Synopsis

This model uses ABM to simulate disease dynamics in Ozamiz City using available data provided by the City Health Office of Ozamiz City. The following goals of this 
repository is to create an environment that uses infection data to simulate infections,
recoveries, deaths, and vaccinations.

## Goals

- Observe relationship between infection rates and vaccination rates. 
  - Effectiveness of mass vaccinations
  - Vaccination hesitancy on vaccinations and infections
- Determine the speed of infection rates and vaccination rates from previous objective.
- Determine if there is observed herd immunity.

## Functions and Design

Because this project uses Agents.jl, and uses the concept of agent-based modeling,
we have to create model step functions that update the environment's parameters.
The model step functions along with the agent step functions.

* Agent step functions update an agent's status:
  * susceptibility
  * infection status
  * recovery
* Model step functions update both model and agent parameters/statuses:
  * vaccinations
  * update probabilities of infection, vaccinations, deaths, recoveries
  * update vaccine hesitancy probabilities
  * update migration rates 

The model uses a graph space to having 51 positions or "nodes" that
represent barangays of Ozamiz City. Caculating the mean R-nought or *R₀* uses
a simplistic naive approach - just count the number of people (agents) from one
infected individual within his infection period.

NOTE: Some of the code is copied from another [repository](https://github.com/JohannesNakayama/EpidemicModel.jl) (but just the CSV writing codes) 
and some ideas were taken 
from examples from the Agents.jl repository

## Limitations

The model lacks some reliable data because of the lack of contact tracing personnel
and also insufficient data collection from the City Health Office of Ozamiz City
because of the priorities on vaccinations and border controls to limit the 
spread of COVID-19. Also, the model needs some improvements on the 
analytical side since the owner lacks the experience on understanding disease dynamics.

## AUTHOR COMMENT

For me, I really find this simulation uncanny because of the lack of information to use.
Not to mention that the assumptions are more on guess estimation and not actual calculation.
Quite a limitation I guess.
