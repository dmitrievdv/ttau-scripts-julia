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

@time stat_pars, stat_names = loadparameters("paper-grid_lowT_RZPsc_stat.dat", 4, 4)
@time nonstat_pars, nonstat_names = loadparameters("paper-grid_lowT_RZPsc_nonstat.dat", 4, 4)

r_mis = [2.0:1:10.0;]
Ws = [1:0.2:4;]
T_maxs = [7000:1000:15000;]
lgṀs = [-11:0.2:-8.4;]
angs = [35:5:60;]

# bound_stat_pars, bound_stat_names = boundpars(stat_pars, stat_names, (1, 10.0^(-11), 10.0^(-9.5)), (2, 1e4, 15e3), (3, 2, 11), (4, 1, 5), (5, 30, 70))
# bound_nonstat_pars, bound_nonstat_names = boundpars(nonstat_pars, nonstat_names, (1, 10.0^(-11), 10.0^(-9.5)), (2, 1e4, 15e3), (3, 2, 11), (4, 1, 5), (5, 35, 60))

gridded_stat_pars, gridded_stat_names = putongrid(stat_pars, stat_names, (lgṀs, (sigfig = 3, axis = :log)), T_maxs, r_mis, Ws, angs); ""
gridded_nonstat_pars, gridded_nonstat_names = putongrid(nonstat_pars, nonstat_names, (lgṀs, (sigfig = 3, axis = :log)), T_maxs, r_mis, Ws, angs); ""
x = correctnonstatgrid!(gridded_nonstat_pars, gridded_nonstat_names, gridded_stat_pars, gridded_stat_names)
x = correctgridforcorotation!(gridded_stat_pars, corotationradius(star))
x = correctgridforcorotation!(gridded_nonstat_pars, corotationradius(star))

grid_stat_pars, grid_stat_names = flattengrid(gridded_stat_pars, gridded_stat_names)
grid_nonstat_pars, grid_nonstat_names = flattengrid(gridded_nonstat_pars, gridded_nonstat_names)
# savepars("paper-grid_lowT_RZPsc_stat", grid_stat_pars, grid_stat_names)
# savepars("paper-grid_lowT_RZPsc_nonstat", grid_nonstat_pars, grid_nonstat_names)

best_stat_pars, best_stat_names = bestmodels(grid_stat_pars, grid_stat_names)
best_nonstat_pars, best_nonstat_names = bestmodels(grid_nonstat_pars, grid_nonstat_names)

function plotnh(pars, names)
    computeltemodels(star, pars, names)
    min_nh, max_nh = findminmaxnh(pars, names)
    mean_nh = findmeannh(pars, names); ""
    plt = scatter(pars[:,2], min_nh, yaxis = :log, label = "min NH")#, xaxis = :log, xlims = (2e-11, 1.2e-10))
    scatter!(plt, pars[:,2], max_nh, label = "max NH")#, xticks = ([1e-11:1e-11:1e-10;], ["1⋅10^{-11}", "", "3⋅10^{-11}", 
                                                                                            #"", "5⋅10^{-11}", "", "7⋅10^{-11}",
                                                                                            #"", "", "1⋅10^{-10}"]))
    plt = scatter!(plt, pars[:,2], mean_nh, yaxis = :log, label = "mean NH", xlims = (minimum(T_maxs)-500, maximum(T_maxs)+500), ylims = (7e7, 1e13))                                                                                      
    plot!(plt, yticks = ([1e8, 3e8, 5e8, 7e8,
                          1e9, 3e9, 5e9, 7e9,
                          1e10,3e10,5e10,7e10,
                          1e11,3e11,5e11,7e11,
                          1e12,3e12,5e12,7e12], 
                         ["10^{8}", "3⋅10^{8}", "5⋅10^{8}", "7⋅10^{8}",
                          "10^{9}", "3⋅10^{9}", "5⋅10^{9}", "7⋅10^{9}",
                          "10^{10}","3⋅10^{10}","5⋅10^{10}","7⋅10^{10}",
                          "10^{11}","3⋅10^{11}","5⋅10^{11}","7⋅10^{11}",
                          "10^{12}","3⋅10^{12}","5⋅10^{12}","7⋅10^{12}"]))
end

plotnh(best_nonstat_pars, best_nonstat_names)

Ca_τ = [2.624 2.624 2.646 2.572 2.333 1.937 1.472 1.041;
        83.41 91.32 86.50 69.17 48.88 32.40 20.39 12.48;
        900.6 929.8 807.7 526.0 284.9 146.6 76.86 41.25]

Na_τ = [3.657e-2 7.795e-2 9.010e-2 6.814e-2 4.492e-2 2.906e-02 1.918e-02 1.304e-02;  
        7.504    5.012    2.523    1.292    0.708    0.415     0.258     0.169;
        180.7    80.91    36.22    17.80    9.519    5.474     3.351     2.162]

τ_T = [8000:1e3:15000;]; τ_NH = [9:1.0:11;]

lgCa_τ_spl2d = Spline2D(τ_NH, τ_T, log10.(Ca_τ), kx = 2)
lgNa_τ_spl2d = Spline2D(τ_NH, τ_T, log10.(Na_τ), kx = 2)

τNa(T, nh) = 10^lgNa_τ_spl2d(log10(nh), T)/1e11*TTauUtils.Models.hydrogenthermalvelocity(T)/1e-4sqrt(23)
τCa(T, nh) = 10^lgCa_τ_spl2d(log10(nh), T)/1e11*TTauUtils.Models.hydrogenthermalvelocity(T)/1e-4/sqrt(40)

function plotCa(pars, names)
    computeltemodels(star, pars, names)
    min_nh, max_nh = findminmaxnh(pars, names)
    mean_nh = findmeannh(pars, names)
    mean_Ca_τ = τCa.(pars[:,2], mean_nh)
    min_Ca_τ = τCa.(pars[:,2], min_nh)
    max_Ca_τ = τCa.(pars[:,2], max_nh)
    plt = scatter(pars[:,2], min_Ca_τ, yaxis = :log, label = "min Ca τ")#, xaxis = :log, xlims = (2e-11, 1.2e-10))
    scatter!(plt, pars[:,2], max_Ca_τ, label = "min Ca τ")#, xticks = ([1e-11:1e-11:1e-10;], ["1⋅10^{-11}", "", "3⋅10^{-11}", 
                                                                                            #"", "5⋅10^{-11}", "", "7⋅10^{-11}",
                                                                                            #"", "", "1⋅10^{-10}"]))
    plt = scatter!(plt, pars[:,2], mean_Ca_τ, yaxis = :log, label = "mean Ca τ", xlims = (minimum(T_maxs)-500, maximum(T_maxs)+500))                                                                                      
    # plot!(plt, yticks = ([1e8, 3e8, 5e8, 7e8,
    #                       1e9, 3e9, 5e9, 7e9,
    #                       1e10,3e10,5e10,7e10,
    #                       1e11,3e11,5e11,7e11,
    #                       1e12,3e12,5e12,7e12], 
    #                      ["10^{8}", "3⋅10^{8}", "5⋅10^{8}", "7⋅10^{8}",
    #                       "10^{9}", "3⋅10^{9}", "5⋅10^{9}", "7⋅10^{9}",
    #                       "10^{10}","3⋅10^{10}","5⋅10^{10}","7⋅10^{10}",
    #                       "10^{11}","3⋅10^{11}","5⋅10^{11}","7⋅10^{11}",
    #                       "10^{12}","3⋅10^{12}","5⋅10^{12}","7⋅10^{12}"]))
end

function plotNa(pars, names)
    computeltemodels(star, pars, names)
    min_nh, max_nh = findminmaxnh(pars, names)
    mean_nh = findmeannh(pars, names)
    mean_Na_τ = τNa.(pars[:,2], mean_nh)
    min_Na_τ = τNa.(pars[:,2], min_nh)
    max_Na_τ = τNa.(pars[:,2], max_nh)
    # plt = scatter(pars[:,2], min_Na_τ, yaxis = :log, label = "min Na τ")#, xaxis = :log, xlims = (2e-11, 1.2e-10))
    # scatter!(plt, pars[:,2], max_Na_τ, label = "min Na τ")#, xticks = ([1e-11:1e-11:1e-10;], ["1⋅10^{-11}", "", "3⋅10^{-11}", 
                                                                                            #"", "5⋅10^{-11}", "", "7⋅10^{-11}",
                                                                                            #"", "", "1⋅10^{-10}"]))
    plt = scatter(pars[:,2], mean_Na_τ, yaxis = :log, label = L"τ_\mathrm{Na}", xlims = (minimum(T_maxs)-500, maximum(T_maxs)+500))
    yticks = [0.001, 0.005, 0.01, 0.05, 0.1, 5, 1, 5, 10]                                                                                      
    plot!(plt, yticks = (yticks, 
                         yticks))
    plot!(plt, xticks = (T_maxs, T_maxs))
end

function plotCaNa(pars, names)
    computeltemodels(star, pars, names)
    mean_nh = findmeannh(pars, names)
    mean_Na_τ = τNa.(pars[:,2], mean_nh)
    computeltemodels(star, pars, names)
    mean_nh = findmeannh(pars, names)
    mean_Ca_τ = τCa.(pars[:,2], mean_nh)
    # plt = scatter(pars[:,2], min_Na_τ, yaxis = :log, label = "min Na τ")#, xaxis = :log, xlims = (2e-11, 1.2e-10))
    # scatter!(plt, pars[:,2], max_Na_τ, label = "min Na τ")#, xticks = ([1e-11:1e-11:1e-10;], ["1⋅10^{-11}", "", "3⋅10^{-11}", 
                                                                                            #"", "5⋅10^{-11}", "", "7⋅10^{-11}",
                                                                                            #"", "", "1⋅10^{-10}"]))
    plt = scatter(pars[:,2], mean_Na_τ, yaxis = :log, label = L"\mathrm{Na}", xlims = (minimum(T_maxs)-500, maximum(T_maxs)+500))
    scatter!(plt, pars[:,2], mean_Ca_τ, yaxis = :log, label = L"\mathrm{Ca}", xlims = (minimum(T_maxs)-500, maximum(T_maxs)+500))
    yticks = [0.001, 0.005, 0.01, 0.05, 0.1, 5, 1, 5, 10]                                                                                      
    plot!(plt, yticks = (yticks, 
                         yticks))
    plot!(plt, xticks = (T_maxs, T_maxs))
    plot!(plt, xlabel = L"T,\ [K]", ylabel = L"\tau")
end

Ca_plt = plotCa(best_nonstat_pars, best_nonstat_names)
Na_plt = plotNa(best_nonstat_pars, best_nonstat_names)
CaNa_plt = plotCaNa(best_stat_pars, best_stat_names)

ind_T9000 = findmodels(grid_nonstat_pars, [9000], [2], [100])
ind_T8000 = findmodels(grid_nonstat_pars, [8000], [2], [100])
ind_T10000 = findmodels(grid_nonstat_pars, [10000], [2], [100])

best_ind_T9000 = ind_T9000[findmin(grid_nonstat_pars[ind_T9000,9])[2]]
best_ind_T8000 = ind_T8000[findmin(grid_nonstat_pars[ind_T8000,9])[2]]
best_ind_T10000 = ind_T10000[findmin(grid_nonstat_pars[ind_T10000,9])[2]]

plt_nonstat_T9000 = plotmodel(grid_nonstat_pars, grid_nonstat_names, best_ind_T9000)
plt_stat_T9000 = plotmodel(grid_stat_pars, grid_stat_names, best_ind_T9000)

plt_nonstat_T8000 = plotmodel(grid_nonstat_pars, grid_nonstat_names, best_ind_T8000)
plt_stat_T8000 = plotmodel(grid_stat_pars, grid_stat_names, best_ind_T8000)

plt_nonstat_T10000 = plotmodel(grid_nonstat_pars, grid_nonstat_names, best_ind_T10000)
plt_stat_T10000 = plotmodel(grid_stat_pars, grid_stat_names, best_ind_T10000)