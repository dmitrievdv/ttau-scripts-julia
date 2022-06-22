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


star = Star("V2129Oph")


obs_file = "spec/V2129Oph/Ha58312.dat"
v_obs, r_obs = readobservation(obs_file)

stat_pars, stat_names = readmodels(star, obs_file, "stat_nonlocal", prof_suffix = "")
nonstat_pars, nonstat_names = readmodels(star, obs_file, "nonstat_nonlocal", prof_suffix = "")


r_mis = [4.0:1:8.0;]
Ws = [0.5:0.5:4;]
T_maxs = [8000:1000:12000;]#:1000:15000;]
lgṀs = [-10:0.2:-8;]
angs = [40:5:80;]

gridded_stat_pars, gridded_stat_names = putongrid(stat_pars, stat_names, (lgṀs, (sigfig = 3, axis = :log)), T_maxs, r_mis, Ws, angs); ""
gridded_nonstat_pars, gridded_nonstat_names = putongrid(nonstat_pars, nonstat_names, (lgṀs, (sigfig = 3, axis = :log)), T_maxs, r_mis, Ws, angs); ""
x = correctnonstatgrid!(gridded_nonstat_pars, gridded_nonstat_names, gridded_stat_pars, gridded_stat_names)
# x = correctgridforcorotation!(gridded_stat_pars, corotationradius(star))
# x = correctgridforcorotation!(gridded_nonstat_pars, corotationradius(star))

grid_stat_pars, grid_stat_names = flattengrid(gridded_stat_pars, gridded_stat_names)
grid_nonstat_pars, grid_nonstat_names = flattengrid(gridded_nonstat_pars, gridded_nonstat_names)
# savepars("V2129Oph_stat", grid_stat_pars, grid_stat_names)
# savepars("V2129Oph_nonstat", grid_nonstat_pars, grid_nonstat_names)

i = findmin(grid_nonstat_pars[:,9])[2]
println(grid_nonstat_pars[i,:])
plt = plotmodel(grid_nonstat_pars, grid_nonstat_names, i)