include("RZPsc.jl")

# include("compobs.jl")
# include("plotobs.jl")
# using StaticArrays
# using LinearAlgebra
# # using TTauUtils
# # using Dierckx
# # using Interpolations
# # using Statistics
# # using LsqFit
# # using Plots
# # using LaTeXStrings
# # using SpecialFunctions
# # using Printf

conf_path = "/home/coloboquito/work/confs/conf-crao-2022/"

# v_obs, r_obs = readobservation("spec/RZPsc_16-11-2013_proc.dat")
# stat_pars, stat_names = readmodels(star, "spec/RZPsc_16-11-2013_proc.dat", "stat_nonlocal", prof_suffix = "phot3crude")
# # stat_pars = addfluxconstant(0.01, stat_pars, stat_names)
# vsini_stat_pars, vsini_stat_names = readmodels(star, "spec/RZPsc_16-11-2013_proc.dat", "stat_nonlocal", prof_suffix = "phot3crude-sini")

# x = 1
# gridded_stat_pars, gridded_stat_names = putongrid(stat_pars, stat_names, (lgṀs, (sigfig = 3, axis = :log)), T_maxs, r_mis, Ws, angs); ""
# gridded_vsini_stat_pars, gridded_vsini_stat_names = putongrid(vsini_stat_pars, vsini_stat_names, (lgṀs, (sigfig = 3, axis = :log)), T_maxs, r_mis, Ws, angs); ""
# x = correctgridforcorotation!(gridded_stat_pars, corotationradius(star))
# x = correctgridforcorotation!(gridded_vsini_stat_pars, corotationradius(star))
plotheat(Matrix(grid_stat_pars), (T_maxs, (index = 2,)), (lgṀs, (index = 1, axis = :log, sigfig = 3)), 
         clims = (0, 0.1), xlabel = L"T_\mathrm{max},\ \mathrm{K}", ylabel = L"\lg\dot{M},\ \mathrm{M_\odot/yr}")
savefig(conf_path*"stat_all_MT.pdf")

plotheat(Matrix(grid_nonstat_pars), (T_maxs, (index = 2,)), (lgṀs, (index = 1, axis = :log, sigfig = 3)), 
         clims = (0, 0.1), xlabel = L"T_\mathrm{max},\ \mathrm{K}", ylabel = L"\lg\dot{M},\ \mathrm{M_\odot/yr}")
savefig(conf_path*"nonstat_all_MT.pdf")

plotheat(Matrix(best_nonstat_pars), (T_maxs, (index = 2,)), (lgṀs, (index = 1, axis = :log, sigfig = 3)), 
clims = (0, 0.1), xlabel = L"T_\mathrm{max},\ \mathrm{K}", ylabel = L"\lg\dot{M},\ \mathrm{M_\odot/yr}")
savefig(conf_path*"nonstat_best_MT.pdf")

plotheat(Matrix(best_stat_pars), (T_maxs, (index = 2,)), (lgṀs, (index = 1, axis = :log, sigfig = 3)), 
clims = (0, 0.1), xlabel = L"T_\mathrm{max},\ \mathrm{K}", ylabel = L"\lg\dot{M},\ \mathrm{M_\odot/yr}")
savefig(conf_path*"stat_best_MT.pdf")

plotheat(Matrix(CaNa_best_nonstat_pars), (T_maxs, (index = 2,)), (lgṀs, (index = 1, axis = :log, sigfig = 3)), 
clims = (0, 0.1), xlabel = L"T_\mathrm{max},\ \mathrm{K}", ylabel = L"\lg\dot{M},\ \mathrm{M_\odot/yr}")
savefig(conf_path*"nonstat_cana_MT.pdf")

plotheat(Matrix(CaNa_best_stat_pars), (T_maxs, (index = 2,)), (lgṀs, (index = 1, axis = :log, sigfig = 3)), 
clims = (0, 0.1), xlabel = L"T_\mathrm{max},\ \mathrm{K}", ylabel = L"\lg\dot{M},\ \mathrm{M_\odot/yr}")
savefig(conf_path*"stat_cana_MT.pdf")

plotheat(Matrix(grid_nonstat_pars), (r_mis, (index = 3,)), (angs, (index = 5,)), 
         clims = (0, 0.1), xlabel = L"r_\mathrm{mi},\ \mathrm{R_{star}}", ylabel = L"i,\ \degree", xlims = (1.5,8.5))
savefig(conf_path*"nonstat_all_RI.pdf")

plotheat(Matrix(best_nonstat_pars), (r_mis, (index = 3,)), (angs, (index = 5,)), 
         clims = (0, 0.1), xlabel = L"r_\mathrm{mi},\ \mathrm{R_{star}}", ylabel = L"i,\ \degree", xlims = (1.5,8.5))
savefig(conf_path*"nonstat_best_RI.pdf")

plotheat(Matrix(CaNa_best_nonstat_pars), (r_mis, (index = 3,)), (angs, (index = 5,)), 
         clims = (0, 0.1), xlabel = L"r_\mathrm{mi},\ \mathrm{R_{star}}", ylabel = L"i,\ \degree", xlims = (1.5,8.5))
savefig(conf_path*"nonstat_cana_RI.pdf")

plotheat(Matrix(grid_nonstat_pars), (r_mis, (index = 3,)), (angs, (index = 5,)), 
         clims = (0, 0.1), xlabel = L"r_\mathrm{mi},\ \mathrm{R_{star}}", ylabel = L"i,\ \degree")
savefig(conf_path*"stat_all_RI.pdf")

plotheat(Matrix(best_nonstat_pars), (r_mis, (index = 3,)), (angs, (index = 5,)), 
         clims = (0, 0.1), xlabel = L"r_\mathrm{mi},\ \mathrm{R_{star}}", ylabel = L"i,\ \degree")
savefig(conf_path*"stat_best_RI.pdf")

plotheat(Matrix(CaNa_best_nonstat_pars), (r_mis, (index = 3,)), (angs, (index = 5,)), 
         clims = (0, 0.1), xlabel = L"r_\mathrm{mi},\ \mathrm{R_{star}}", ylabel = L"i,\ \degree")
savefig(conf_path*"stat_cana_RI.pdf")

δ, id = findmin(stat_pars[:,9])
plotmodel(stat_pars, stat_names, id, xlabel = L"v,\ \mathrm{km/s}", ylabel = L"I/I_c")
savefig(conf_path*"best_stat_prof.pdf")

δ, id = findmin(CaNa_best_nonstat_pars[:,9])
plotmodel(CaNa_best_nonstat_pars, CaNa_best_nonstat_names, id, xlabel = L"v,\ \mathrm{km/s}", ylabel = L"I/I_c")
savefig(conf_path*"best_prof.pdf")