include("compobs.jl")
include("plotobs.jl")
# using TTauUtils
# using Dierckx
# using Interpolations
# using Statistics
# using LsqFit
# using Plots
# using LaTeXStrings
# using SpecialFunctions
# using Printf


star = Star("RZPsc")

v_obs, r_obs = readobservation("spec/RZPsc_16-11-2013_proc.dat")
# stat_pars, stat_names = readmodels(star, "spec/RZPsc_16-11-2013_proc.dat", "stat_nonlocal", prof_suffix = "phot3crude")
# nonstat_pars, nonstat_names = readmodels(star, "spec/RZPsc_16-11-2013_proc.dat", "nonstat_nonlocal", prof_suffix = "phot3crude")
# δs = pars[:,8]

@time stat_pars, stat_names = loadparameters("paper-grid_RZPsc_stat.dat", 4, 4)
@time nonstat_pars, nonstat_names = loadparameters("paper-grid_RZPsc_nonstat.dat", 4, 4)

r_mis = [2.0:1:10.0;]
Ws = [1:0.2:4;]
T_maxs = [10000:1000:15000;]
lgṀs = [-11:0.1:-9.5;]
angs = [30:5:60;]

bound_stat_pars, bound_stat_names = boundpars(stat_pars, stat_names, (1, 10.0^(-11), 10.0^(-9.5)), (2, 1e4, 15e3), (3, 2, 11), (4, 1, 5), (5, 30, 70))
bound_nonstat_pars, bound_nonstat_names = boundpars(nonstat_pars, nonstat_names, (1, 10.0^(-11), 10.0^(-9.5)), (2, 1e4, 15e3), (3, 2, 11), (4, 1, 5), (5, 35, 60))

gridded_stat_pars, gridded_stat_names = putongrid(stat_pars, stat_names, (lgṀs, (sigfig = 3, axis = :log)), T_maxs, r_mis, Ws, angs); ""
gridded_nonstat_pars, gridded_nonstat_names = putongrid(nonstat_pars, nonstat_names, (lgṀs, (sigfig = 3, axis = :log)), T_maxs, r_mis, Ws, angs); ""
x = correctnonstatgrid!(gridded_nonstat_pars, gridded_nonstat_names, gridded_stat_pars, gridded_stat_names)
x = correctgridforcorotation!(gridded_stat_pars, corotationradius(star))
x = correctgridforcorotation!(gridded_nonstat_pars, corotationradius(star))

grid_stat_pars, grid_stat_names = flattengrid(gridded_stat_pars, gridded_stat_names)
grid_nonstat_pars, grid_nonstat_names = flattengrid(gridded_nonstat_pars, gridded_nonstat_names)
savepars("paper-grid_RZPsc_stat", grid_stat_pars, grid_stat_names)
savepars("paper-grid_RZPsc_nonstat", grid_nonstat_pars, grid_nonstat_names)