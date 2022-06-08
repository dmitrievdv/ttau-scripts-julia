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
stat_pars, stat_names = readmodels(star, "spec/RZPsc_16-11-2013_proc.dat", "stat_nonlocal", prof_suffix = "phot3crude")
nonstat_pars, nonstat_names = readmodels(star, "spec/RZPsc_16-11-2013_proc.dat", "nonstat_nonlocal", prof_suffix = "phot3crude")
# δs = pars[:,8]

function calcprofile(star, model_name, profile_name, angle, line; args...)
    does_model_exist = checkmodel(model_name, star)
    model = if does_model_exist
        loadmodel(star, model_name)
    else
        throw(ErrorException("model does not exist!"))
    end
    profile = HydrogenProfileDoppler(model, line[1], line[2], angle, 0.03, 0.03, 0.01, 50, blue_v_max = 300, red_v_max = 600)
    saveprofile(profile, linename(line)*'_'*string(floor(Int, angle)))
    profile = TTauUtils.addphotosphespecdoppler(profile, 0.03, "spec/RZ_Psc_Ha_syn_unwid_corr.dat")
    saveprofile(profile, profile_name)
end

@time stat_pars, stat_names = loadparameters("paper-grid_RZPsc_stat.dat", 4, 4)
@time nonstat_pars, nonstat_names = loadparameters("paper-grid_RZPsc_nonstat.dat", 4, 4)

r_mis = [2.0:1:10.0;]
Ws = [1:0.2:4;]
T_maxs = [10000:1000:15000;]
lgṀs = [-11:0.1:-9.5;]
angs = [35:5:60;]

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

best_stat_pars, best_stat_names = bestmodels(stat_pars, stat_names)
best_nonstat_pars, best_nonstat_names = bestmodels(nonstat_pars, nonstat_names)

computeltemodels(star, best_stat_pars, best_stat_names)
min_nh, max_nh = findminmaxnh(best_stat_pars, best_stat_names)
mean_nh = findmeannh(best_stat_pars, best_stat_names); ""

plt = scatter(best_stat_pars[:,2], min_nh, yaxis = :log, label = "min NH")#, xaxis = :log, xlims = (2e-11, 1.2e-10))
scatter!(plt, best_stat_pars[:,2], max_nh, label = "max NH")#, xticks = ([1e-11:1e-11:1e-10;], ["1⋅10^{-11}", "", "3⋅10^{-11}", 
                                                                                            #"", "5⋅10^{-11}", "", "7⋅10^{-11}",
                                                                                            #"", "", "1⋅10^{-10}"]))
plt = scatter!(plt, best_stat_pars[:,2], mean_nh, yaxis = :log, label = "mean NH")                                                                                      
plot!(plt, yticks = ([1e8,3e8,5e8,7e8,1e9,3e9,5e9,7e9,1e10,3e10,5e10,7e10], ["10^{8}", "3⋅10^{8}", "5⋅10^{8}", "7⋅10^{8}",
                                                                            "10^{9}", "3⋅10^{9}", "5⋅10^{9}", "7⋅10^{9}",
                                                                            "10^{10}","3⋅10^{10}","5⋅10^{10}","7⋅10^{10}"]))