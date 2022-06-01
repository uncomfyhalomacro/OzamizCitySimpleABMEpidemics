total_susceptible(m::AgentBasedModel) = count(a.status == (:S) for a in allagents(m))
total_infected(m::AgentBasedModel) = count(a.status == (:I) for a in allagents(m))
total_recovered(m::AgentBasedModel) = count(a.status == (:R) for a in allagents(m))
function total_partially_vaccinated(m::AgentBasedModel)
    count(a.vaccination_status == (:V) && a.num_vac_shots == 1 for a in allagents(m))
end
function total_fully_vaccinated_without_booster(model::ABM)
    count(a.vaccination_status == (:V) && a.num_vac_shots == 2 for a in allagents(model))
end
function total_fully_vaccinated_with_booster(model::ABM)
    count(a.vaccination_status == (:V) && a.num_vac_shots == 3 for a in allagents(model))
end
total_dead(m::AgentBasedModel) = OZAMIZ_CITY_POPULATION - length(m.agents)

function run_data_collect!(model::ABM, steps::Int64)
    @info "Starting data collection!"
    agent_data, model_data = Agents.run!(model,
                                         agent_step!,
                                         model_step!,
                                         steps;
                                         adata = [
                                         # :status,
                                         # :vaccination_status,
                                         # :days_vaccinated,
                                         # :num_vac_shots,
                                         # :date_infected,
                                         # :days_infected,
                                         # :days_recovered,
                                         # :from_csv_data,
                                         ],
                                         mdata = [
                                             :days_passed,
                                             :susceptible,
                                             :infected,
                                             :recovered,
                                             :partially_vaccinated,
                                             :fully_vaccinated_without_booster,
                                             :fully_vaccinated_with_booster,
                                             :dead,
                                             :mean_R₀,
                                             :median_R₀,
                                             :infection_probability,
                                             :vaccination_probability,
                                         ],
                                         agents_first = true)
    return agent_data, model_data
end

function run_model()
    for n_experiments in 1:3
        @info "Running the model! - $(n_experiments)"
        params = initialize_parameters(;
                                       number_of_barangays = 51,
                                       max_travel_rate = 0.06,
                                       vaccination_hesitancy_probability_range = 0.001:0.02:0.6,
                                       undetected_probability_range = 0.3:0.02:0.6,
                                       infection_probability = 0.46,
                                       vaccination_probability = 0.09788)
        model = initialize_model(; params...)
        @info "Running $(1000 * n_experiments) steps/days of simulation"
        _, model_dataframe = run_data_collect!(model, 1000 * n_experiments)
        #########
        @info "Generating log plot!"
        @show model_dataframe[1:10, :]
        x = model_dataframe.step
        fig = Figure(; resolution = (1080, 720))
        ax = fig[1, 1] = Axis(fig; xlabel = "steps or days", ylabel = "log10(count)")
        ls = lines!(ax, x, log10.(model_dataframe[:, :susceptible]); color = :orange)
        li = lines!(ax, x, log10.(model_dataframe[:, :infected]); color = :yellow)
        lr = lines!(ax, x, log10.(model_dataframe[:, :recovered]); color = :red)
        lv = lines!(ax, x, log10.(model_dataframe[:, :fully_vaccinated_with_booster]);
                    color = :blue)
        ld = lines!(ax, x, log10.(model_dataframe[:, :dead]); color = :green)
        lmean_R₀ = lines!(ax, x, log10.(model_dataframe[:, :mean_R₀]); color = :black)
        lmedian_R₀ = lines!(ax, x, log10.(model_dataframe[:, :median_R₀]); color = :grey)
        fig[1, 2] = Legend(fig,
                           [ls, li, lr, lv, ld, lmean_R₀, lmedian_R₀],
                           [
                               "susceptible",
                               "infected",
                               "recovered",
                               "vaccinated",
                               "dead",
                               "mean R₀",
                               "median R₀",
                           ];
                           textsize = 12)
        save_plot_path = joinpath(dirname(Base.current_project(@__DIR__)),
                                  "plots_and_animations")
        mkpath(save_plot_path)
        save(joinpath(save_plot_path, "generated_log_plot_$(n_experiments).png"), fig)
        @info "Done generating log plot!"
        ############
        @info "Generating normal plot!"
        fig2 = Figure(; resolution = (1080, 720))
        ax = fig2[1, 1] = Axis(fig2; xlabel = "steps or days", ylabel = "count")
        ls = lines!(ax, x, model_dataframe[:, :susceptible]; color = :orange)
        li = lines!(ax, x, model_dataframe[:, :infected]; color = :yellow)
        lr = lines!(ax, x, model_dataframe[:, :recovered]; color = :red)
        lv = lines!(ax, x, model_dataframe[:, :fully_vaccinated_with_booster];
                    color = :blue)
        # dead = log10.(N .- model_dataframe[:, aggname(:status, length)])
        ld = lines!(ax, x, model_dataframe[:, :dead]; color = :green)
        # lmean_R₀ = lines!(ax, x, model_dataframe[:, :mean_R₀]; color = :black)
        # lmedian_R₀ = lines!(ax, x, model_dataframe[:, :median_R₀]; color = :grey)
        fig2[1, 2] = Legend(fig2,
                            [ls, li, lr, lv, ld],
                            [
                                "susceptible",
                                "infected",
                                "recovered",
                                "vaccinated",
                                "dead",
                            ];
                            textsize = 13)
        save(joinpath(dirname(Base.current_project()),
                      "plots_and_animations",
                      "generated_plot_$(n_experiments).png"),
             fig2)
        @info "Done generating normal plot!"
        # NOTE: Commenting this out because of the large memory it eats during the runs
        # @info "Writing agent data!"
        # # write_dataframe(;
        #     df=agent_dataframe,
        #     filename="agent_data_$(n_experiments).csv",
        #     path=joinpath(dirname(Base.current_project()), "data"),
        # )
        @info "Writing model data!"
        write_dataframe(;
                        df = model_dataframe,
                        filename = "model_data_$(n_experiments).csv",
                        path = joinpath(dirname(Base.current_project()), "data"))
        @info "Done running the simulation! - $(n_experiments)"
    end
    @info "Done running all simulations!"
end

function write_dataframe(; df::DataFrames.DataFrame, filename::String, path::String)
    symbol_to_string!(df)
    if endswith(filename, ".csv")
        CSV.write(joinpath(path, filename), df)
    else
        CSV.write(joinpath(path, filename * ".csv"), df)
    end
end

function symbol_to_string!(df::DataFrames.DataFrame)
    for colname in DataFrames.names(df)
        if typeof(df[!, colname]) == Vector{Symbol}
            df[!, colname] = string.(df[!, colname])
        end
    end
    return df
end
