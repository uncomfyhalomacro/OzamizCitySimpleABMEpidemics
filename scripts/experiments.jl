import Pkg
Pkg.add([
            "CSV",
            "DataFrames",
            "Plots",
            "StatsPlots",
            "Statistics",
            "Measures",
            "Distributions",
        ])
using CSV
using Plots
using DataFrames
using StatsPlots
using Statistics
using Distributions
using Measures

gr()

df_file_path = "/home/uncomfy/Projects/OzamizCitySimpleABMEpidemics/data/model_data_3.csv"

df = CSV.read(df_file_path, DataFrame)

function plot_infected(d::DataFrame)
    params = (;
              title = "Simulated infected cases in Ozamiz City",
              xlabel = "Dates",
              ylabel = "Number of infected cases",
              xticks = df.days_passed[begin:100:end],
              label = "infected",
              minorgrid = true,
              minorgridstyle = :dot,
              minorgridwidth = 3,
              size = (900, 900),
              xrot = 45,
              lw = 3)

    infected_plot = plot(d.days_passed, d.infected; params...)
    savefig(infected_plot, "infected_plot_1.png")
    log10_data = log10.(d.infected)
    log_infected_plot = plot(d.days_passed, log10_data; params...)
    savefig(log_infected_plot, "infected_plot_log_1.png")
    combine_plot_infected = plot(infected_plot, log_infected_plot, layout = (2, 1))
    savefig(combine_plot_infected, "combine_plot_infected.png")
    return d
end

function plot_recovered(d::DataFrame)
    params = (;
              title = "Simulated recoveries in Ozamiz City",
              xlabel = "Dates",
              ylabel = "Number of recoveries",
              xticks = df.days_passed[begin:100:end],
              label = "recovered",
              minorgrid = true,
              minorgridstyle = :dot,
              minorgridwidth = 3,
              size = (900, 900),
              xrot = 45,
              lw = 3)

    recoveries_plot = plot(d.days_passed, d.recovered; params...)
    savefig(recoveries_plot, "recoveries_plot_1.png")
    log10_data = log10.(d.recovered)
    log_recoveries_plot = plot(d.days_passed, log10_data; params...)
    savefig(log_recoveries_plot, "recoveries_plot_log_1.png")
    combine_plot_recovered = plot(recoveries_plot, log_recoveries_plot,
                                  layout = (2, 1))
    savefig(combine_plot_recovered,
            "combine_plot_recovered.png")
    return d
end

function plot_dead(d::DataFrame)
    params = (;
              title = "Simulated dead in Ozamiz City",
              xlabel = "Dates",
              ylabel = "Number of dead",
              xticks = df.days_passed[begin:100:end],
              label = "dead",
              minorgrid = true,
              minorgridstyle = :dot,
              minorgridwidth = 3,
              size = (900, 900),
              xrot = 45,
              lw = 3)

    dead_plot = plot(d.days_passed, d.dead; params...)
    savefig(dead_plot, "dead_plot_1.png")
    log10_data = log10.(d.dead)
    dead_plot = plot(d.days_passed, log10_data; params...)
    savefig(dead_plot, "dead_plot_log_1.png")
    return d
end

function plot_vaccinated(d::DataFrame)
    vaccination_list = [d.partially_vaccinated d.fully_vaccinated_without_booster d.fully_vaccinated_with_booster]
    vaccination_labels = ["partial" "full without booster" "full with booster"]
    params = (;
              title = "Simulated vaccinations in Ozamiz City",
              xlabel = "Dates",
              ylabel = "Number of vaccinations",
              xticks = df.days_passed[begin:100:end],
              label = vaccination_labels,
              minorgrid = true,
              minorgridstyle = :dot,
              minorgridwidth = 3,
              size = (900, 900),
              xrot = 45,
              lw = 3)
    vaccinated_plot = plot(d.days_passed, vaccination_list; params...)
    savefig(vaccinated_plot, "vaccinated_plot_1.png")
    log10_data = log10.(vaccination_list)
    log_vaccinated_plot = plot(d.days_passed, log10_data; params...)
    savefig(log_vaccinated_plot, "vaccinated_plot_log_1.png")
    combine_plot_vaccinated = plot(vaccinated_plot, log_vaccinated_plot, layout = (2, 1))
    savefig(combine_plot_vaccinated, "combine_plot_vaccinated.png")
    return d
end

function plot_sourced_data(d::DataFrame)
    another_df_path = "/home/uncomfy/Projects/OzamizCitySimpleABMEpidemics/data/cumulative_sum_you_should_groupby_date_later.csv"
    another_df = CSV.read(another_df_path, DataFrame)
    another_df_groupby_date = groupby(another_df, [:DATES]; sort = true)
    length(another_df_groupby_date)
    collect_infected_cases = [maximum(i.CUM_CASES) for i in another_df_groupby_date]
    length(collect_infected_cases)
    set_dates = Set(another_df.DATES[begin:end]) |> collect |> sort

    maximum(collect_infected_cases)
    log_collect_infection_cases = log10.(collect_infected_cases)

    params = (;
              title = "Plot of sourced infected data in Ozamiz City",
              xlabel = "Dates",
              ylabel = "Number of infected",
              label = "infected",
              minorgrid = true,
              minorgridstyle = :dot,
              minorgridwidth = 3,
              size = (900, 900),
              xrot = 45,
              lw = 3)
    infected_plot = plot(set_dates, collect_infected_cases; params...)
    savefig(infected_plot, "infected_plot_from_source_data_1.png")
    infected_plot = plot(set_dates, log_collect_infection_cases; params...)
    savefig(infected_plot, "infected_plot_log_from_source_data_1.png")
    params = (;
              title = "Overlay plot of sourced and simulated infected data in Ozamiz City",
              xlabel = "Dates",
              ylabel = "Number of infected",
              label = ["simulated" "sourced"],
              minorgrid = true,
              minorgridstyle = :dot,
              minorgridwidth = 3,
              size = (900, 900),
              xrot = 45,
              lw = 3)
    infected_plot = plot(set_dates,
                         [d.infected[begin:length(collect_infected_cases)] collect_infected_cases];
                         params...)

    savefig(infected_plot, "overlay_infected_plot_from_source_data_and_simulated_1.png")
    infected_plot = plot(set_dates,
                         log10.([d.infected[begin:length(collect_infected_cases)] collect_infected_cases]);
                         params...)
    savefig(infected_plot, "overlay_infected_log_plot_from_source_data_and_simulated_1.png")
end

function plot_mean_R₀(df::DataFrame)
    params = (;
              title = "Plotted R₀ in Ozamiz City",
              xlabel = "Dates",
              ylabel = "Mean R₀",
              label = "R₀",
              minorgrid = true,
              minorgridstyle = :dot,
              minorgridwidth = 3,
              size = (900, 900),
              xrot = 45,
              lw = 3)
    brn_plot = scatter(df.days_passed, df.mean_R₀,
                       smooth = :true, linecolor = :red, linewidth = 3; params...)
    # brn_plot = Plots.abline!(bhat_brn..., label = "trendline")
    savefig(brn_plot, "brn_plot_1.png")
    log_brn_plot = scatter(df.days_passed, log10.(df.mean_R₀);
                           params...)
    savefig(log_brn_plot, "log_brn_plot_1.png")
    combine_plot_mean_R₀ = plot(brn_plot, log_brn_plot, layout = (2, 1))
    savefig(combine_plot_mean_R₀, "combine_plot_mean_R.png")
end

function plot_mean_number_of_infections_per_day(df::DataFrame)
    another_df_path = "/home/uncomfy/Projects/OzamizCitySimpleABMEpidemics/data/cumulative_sum_you_should_groupby_date_later.csv"
    another_df = CSV.read(another_df_path, DataFrame)
    another_df_groupby_date = groupby(another_df, [:DATES]; sort = true)
    length(another_df_groupby_date)
    collect_infected_cases = [maximum(i.CUM_CASES) for i in another_df_groupby_date]
    length(collect_infected_cases)
    set_dates = Set(another_df.DATES[begin:end]) |> collect |> sort
    mean_new_infections_per_day = mean(diff(vec(collect_infected_cases)))
    maximum(collect_infected_cases)

    params = (;
              title = "Plotted new number of infections per day from source data in Ozamiz City",
              xlabel = "Days",
              ylabel = "Mean new infections",
              label = "new infections",
              minorgrid = false,
              size = (900, 900),
              xrot = 45)

    uwu = bar(1:361,
              collect(a ≤ 0.0 ? 0.0 : a for a in [4, diff(collect_infected_cases)...]);
              params...)
    savefig(uwu, "number_of_new_cases_per_day.png")
    return uwu
end

function plot_everything(df::DataFrame)
    column_list = [df.infected df.recovered df.partially_vaccinated df.fully_vaccinated_without_booster df.fully_vaccinated_with_booster]
    labels = ["infected" "recoveries" "partially vaccinated" "fully vaccinated without booster" "fully vaccinated with booster"]
    params = (;
              title = "Simulated infected, recoveries, and  vaccinations in Ozamiz City",
              xlabel = "Dates",
              ylabel = "Count",
              xticks = df.days_passed[begin:100:end],
              legend_position = :topleft,
              label = labels,
              minorgrid = true,
              minorgridstyle = :dot,
              minorgridwidth = 3,
              size = (900, 900),
              xrot = 45,
              lw = 3)
    plot_val = plot(df.days_passed, column_list; params...)
    log_plot_val = plot(df.days_passed,
                        log10.(column_list);
                        params...)

    savefig(plot_val, "all_plot.png")
    savefig(log_plot_val, "log_plot_val.png")
end

plot_infected(df), plot_recovered(df), plot_dead(df), plot_vaccinated(df),
plot_sourced_data(df), plot_mean_R₀(df), plot_mean_number_of_infections_per_day(df),
plot_everything(df)
