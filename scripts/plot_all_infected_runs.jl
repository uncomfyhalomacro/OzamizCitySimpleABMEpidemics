using Distributions: StatsBase
import Pkg
Pkg.add([
            "CSV",
            "DataFrames",
            "Plots",
            "StatsPlots",
            "Statistics",
            "Measures",
            "Distributions",
            "StatsBase",
        ])
using CSV
using Plots
using DataFrames
using StatsPlots
using Statistics
using Distributions
using Measures

gr()

datadir = joinpath(@__DIR__, "../data")
another_df_path = joinpath(datadir, "cumulative_sum_you_should_groupby_date_later.csv")
another_df = CSV.read(another_df_path, DataFrame)
another_df_groupby_date = groupby(another_df, [:DATES]; sort = true)
length(another_df_groupby_date)
collect_infected_cases = [maximum(i.CUM_CASES) for i in another_df_groupby_date]
length(collect_infected_cases)
set_dates = Set(another_df.DATES[begin:end]) |> collect |> sort
df_file_path_1 = joinpath(datadir, "model_data_1.csv")
df_file_path_2 = joinpath(datadir, "model_data_2.csv")
df_file_path_3 = joinpath(datadir, "model_data_3.csv")
df1 = CSV.read(df_file_path_1, DataFrame)
df2 = CSV.read(df_file_path_2, DataFrame)
df3 = CSV.read(df_file_path_3, DataFrame)

vaccination_labels = ["infected" "partial vaccination" "full without booster" "full with booster"]
params = (;
          title = ["Simulated infected cases in Ozamiz City (at 1000 steps)" "at 2000 steps" "at 3000 steps"],
          xlabel = ["" "" "Dates"],
          ylabel = ["" "Count" ""],
          labels = vaccination_labels,
          legend_position = :bottomright,
          minorgrid = true,
          minorgridstyle = :dot,
          minorgridwidth = 3,
          size = (900, 900),
          xrot = 45,
          lw = 3)
plt1 = plot(df1.days_passed,
            log10.([df1.infected df1.partially_vaccinated df1.fully_vaccinated_without_booster df1.fully_vaccinated_with_booster]))
plt2 = plot(df2.days_passed,
            log10.([df2.infected df2.partially_vaccinated df2.fully_vaccinated_without_booster df2.fully_vaccinated_with_booster]))
plt3 = plot(df3.days_passed,
            log10.([df3.infected df3.partially_vaccinated df3.fully_vaccinated_without_booster df3.fully_vaccinated_with_booster]))
combine_all = plot(plt1, plt2, plt3, layout = (3, 1); params...)
savefig(combine_all, joinpath(@__DIR__, "combine_all_infected_vaccinatio_simulations.png"))
@info length(set_dates)
plt1 = plot(set_dates, [collect_infected_cases, df1.infected[begin:length(set_dates)]])
plt2 = plot(set_dates, [collect_infected_cases, df2.infected[begin:length(set_dates)]])
plt3 = plot(set_dates, [collect_infected_cases, df3.infected[begin:length(set_dates)]])
stats = StatsBase.mean_and_var.([
                                    filter(>(0),
                                           diff(df1.infected[begin:length(set_dates)])),
                                    filter(>(0),
                                           diff(df2.infected[begin:length(set_dates)])),
                                    filter(>(0),
                                           diff(df3.infected[begin:length(set_dates)])),
                                    filter(>(0), diff(collect_infected_cases)),
                                ])
stats2 = StatsBase.mean_and_std.([
                                     filter(>(0),
                                            diff(df1.infected[begin:length(set_dates)])),
                                     filter(>(0),
                                            diff(df2.infected[begin:length(set_dates)])),
                                     filter(>(0),
                                            diff(df3.infected[begin:length(set_dates)])),
                                     filter(>(0), diff(collect_infected_cases)),
                                 ])
stats3 = StatsBase.mean_and_std.([
                                     filter(>(0),
                                            diff(df1.partially_vaccinated)),
                                     filter(>(0),
                                            diff(df2.partially_vaccinated)),
                                     filter(>(0),
                                            diff(df3.partially_vaccinated)),
                                 ])
stats4 = StatsBase.mean_and_std.([
                                     filter(>(0),
                                            diff(df1.fully_vaccinated_without_booster)),
                                     filter(>(0),
                                            diff(df2.fully_vaccinated_without_booster)),
                                     filter(>(0),
                                            diff(df3.fully_vaccinated_without_booster)),
                                 ])
stats5 = StatsBase.mean_and_std.([
                                     filter(>(0),
                                            diff(df1.fully_vaccinated_with_booster)),
                                     filter(>(0),
                                            diff(df2.fully_vaccinated_with_booster)),
                                     filter(>(0),
                                            diff(df3.fully_vaccinated_with_booster)),
                                 ])
@info stats
@info stats2
@info "partial" stats3
@info "fw/ob" stats4
@info "fullwith_booster" stats5
params = (;
          title = ["Plot of sourced and simulated infected cases in Ozamiz City (vs 1st run)" "source data vs 2nd run" "source data vs 3rd run"],
          xlabel = ["" "" "Dates (361 days)"],
          ylabel = ["" "Number of infected cases" ""],
          labels = ["sourced" "simulated"],
          legend_position = :bottomright,
          minorgrid = true,
          minorgridstyle = :dot,
          minorgridwidth = 3,
          size = (900, 900),
          xrot = 45,
          lw = 3)
combine_all = plot(plt1, plt2, plt3, layout = (3, 1); params...)
savefig(combine_all, joinpath(@__DIR__, "overlay_combine_all_infected_simulations.png"))
params = (;
          title = ["Plotted R₀ in Ozamiz City (1st run)" "2nd run" "3rd run"],
          xlabel = ["" "" "Dates"],
          ylabel = ["" "Mean R₀" ""],
          label = "R₀",
          minorgrid = true,
          minorgridstyle = :dot,
          minorgridwidth = 3,
          size = (900, 900),
          xrot = 45,
          lw = 3)

plt1 = scatter(df1.days_passed, df1.mean_R₀, smooth = :true, linecolor = :red, linewidth = 3)
plt2 = scatter(df2.days_passed, df2.mean_R₀,
            smooth = :true, linecolor = :red, linewidth = 3)
plt3 = scatter(df3.days_passed, df3.mean_R₀,
            smooth = :true, linecolor = :red, linewidth = 3)

combine_all = plot(plt1, plt2, plt3, layout = (3, 1); params...)

savefig(combine_all, joinpath(@__DIR__,"overlay_combine_all_brn_simulations.png"))
