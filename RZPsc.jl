include("compobs.jl")
include("plotobs.jl")
using StaticArrays
using LinearAlgebra
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
# stat_pars = addfluxconstant(0.01, stat_pars, stat_names)
nonstat_pars, nonstat_names = readmodels(star, "spec/RZPsc_16-11-2013_proc.dat", "nonstat_nonlocal", prof_suffix = "phot3crude")
# nonstat_pars = addfluxconstant(0.01, nonstat_pars, nonstat_names)
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

# @time stat_pars, stat_names = loadparameters("paper-grid_absonly_lowT_RZPsc_stat.dat", 4, 4)
# @time nonstat_pars, nonstat_names = loadparameters("paper-grid_absonly_lowT_RZPsc_nonstat.dat", 4, 4)

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
savepars("paper-grid_lowT_RZPsc_stat", grid_stat_pars, grid_stat_names)
savepars("paper-grid_lowT_RZPsc_nonstat", grid_nonstat_pars, grid_nonstat_names)

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

Ca_τ = [6.84   6.14   6.15   6.73   7.73   7.97   7.01   5.42   3.72   2.41;
        1.67e2 1.55e2 1.85e2 1.96e2 1.69e2 1.18e2 7.31e1 4.33e1 2.58e1 1.56e1;
        1.32e3 1.60e3 1.81e3 1.85e3 1.54e3 9.14e2 4.40e2 2.11e2 1.06e2 5.47e1;]

Na_τ = [5.48e-2 5.21e-2 1.08e-1 2.54e-1 2.57e-1 1.72e-1 1.05e-1 5.49e-2 4.17e-2 2.79e-2;
        1.86    1.08e1  1.91e1  1.09e1  5.21    2.64    1.44    8.46e-1 5.26e-1 3.44e-1;
        2.04e2  5.79e2  3.47e2  1.60e2  7.21e1  3.55e1  1.91e1  1.10e1  6.74    4.34;]

τ_T = [6000:1e3:15000;]; τ_NH = [9:1.0:11;]
T_min_τ = minimum(τ_T); T_max_τ = maximum(τ_T)
lgNH_min_τ = minimum(τ_NH); lgNH_max_τ = maximum(τ_NH)

lgCa_τ_spl2d = Spline2D(τ_NH, τ_T, log10.(Ca_τ), kx = 2)
lgNa_τ_spl2d = Spline2D(τ_NH, τ_T, log10.(Na_τ), kx = 2)

function τNa(T, nh) 
    lgNH = log10(nh)
    T_τ = 0.0
    lgNH_τ = 0.0
    if (T > T_min_τ) & (T < T_max_τ) & (lgNH > lgNH_min_τ) & (lgNH < lgNH_max_τ)
        T_τ = T
        lgNH_τ = lgNH
    else
        lgNH_τ = τ_NH[findmin(abs.(τ_NH .- lgNH))[2]]
        T_τ = τ_T[findmin(abs.(τ_T .- T))[2]]
    end
    10^lgNa_τ_spl2d(lgNH_τ, T_τ)/TTauUtils.starradiiincm(star)*TTauUtils.Models.hydrogenthermalvelocity(T)/1e-4/sqrt(23)
end

function τCa(T, nh) 
    lgNH = log10(nh)
    T_τ = 0.0
    lgNH_τ = 0.0
    if (T > T_min_τ) & (T < T_max_τ) & (lgNH > lgNH_min_τ) & (lgNH < lgNH_max_τ)
        T_τ = T
        lgNH_τ = lgNH
    else
        lgNH_τ = τ_NH[findmin(abs.(τ_NH .- lgNH))[2]]
        T_τ = τ_T[findmin(abs.(τ_T .- T))[2]]
    end
    10^lgCa_τ_spl2d(log10(nh), T)/TTauUtils.starradiiincm(star)*TTauUtils.Models.hydrogenthermalvelocity(T)/1e-4/sqrt(40)
end

function calcabsprofiledip(pars :: Matrix{Float64}, names, τfunc, atomic_mass, v_z, dr, max_dz; star_dir = "stars")
    n_models = length(names)
    abs_prof_dip = zeros(n_models)
    for i = 1:n_models
        abs_prof_dip[i] = calcabsprofiledip(pars[i,:], names[i], τfunc, atomic_mass, v_z, dr, max_dz, star_dir = star_dir)
    end
    abs_prof_dip
end

function calcabsprofiledip(par :: Vector{Float64}, name :: Vector{S}, τfunc, atomic_mass, v_z, dr, max_dz; star_dir = "stars") where S <: AbstractString
    Ṁ, T_max, R_in, W = par[1:4]
    model_name, profile_name = name
    println("$model_name $profile_name")
    lte_name = join(split(model_name, '_')[1:3], '_')*"_lte"
    model = if !isfile("$star_dir/$(star.name)/$lte_name/$lte_name.dat")
        model = TTauUtils.Models.SolidMagnetosphereNHCoolLTE(model_name, star, R_in, R_in + W, Ṁ, T_max, 10)
        savemodel(model)
        model
    else
        loadmodel(star, lte_name)
    end
    angle = par[5]
    calcabsprofiledip(model, angle, τfunc, atomic_mass, v_z, dr, max_dz)
end

function calcabsprofiledip(model :: TTauUtils.HydrogenModel, angle :: Real, τfunc, atomic_mass, v, dr, max_dz)
    grid = TTauUtils.Grids.polargrid(1, dr)
    star = model.star
    R_s = TTauUtils.starradiiincm(star)
    kin = model.kinematics
    orientation = TTauUtils.GeometryAndOrientations.Orientation(angle)
    S_star = 0.0
    S_tau = 0.0
    λ_s = 0.1
    for (x, y, dS) in zip(grid...)
        borders = TTauUtils.GeometryAndOrientations.calcborders(x, y, model.geometry, orientation)
        borders_n = length(borders)
        depth = 0.0
        τ_ray = 0.0
        for i=1:borders_n
            depth += (2*((i+1)%2 - 1) + 1)*borders[i]
        end
        min_dz = max_dz/1e2
        if depth > 1e-5
            for i=borders_n:-2:2
                # dz = depth/m_z
                z_in = borders[i]
                z_out = borders[i-1]
                last = false
                z_end = z_in
                z = z_end
                r = √(x^2 + y^2 + z^2)
                r̂ = SA[x/r, y/r, z/r]
                θ = acos(r̂ ⋅ orientation.dipole_axis)
                x_1, x_2 = TTauUtils.spheretogrid(model.geometry, r, θ)
                T_e = TTauUtils.gridTe(model, x_1, x_2)
                v_t = TTauUtils.Models.hydrogenthermalvelocity(T_e)/√(atomic_mass)
                r̂, θ̂, ϕ̂ = TTauUtils.Profiles.sphericalbasis(x, y, z, r, θ, orientation)
                R_to_star_axis = r*sin(acos(r̂ ⋅ orientation.star_axis))
                sphere_v = TTauUtils.Models.velocityfieldspherical(kin, r, θ)
                sphere_v += SA[0.0, 0.0, R_to_star_axis*star.v_eq*1e5]
                v_z = -(r̂[3]*sphere_v[1] + θ̂[3]*sphere_v[2] + ϕ̂[3]*sphere_v[3])
                l_s = TTauUtils.Profiles.sobolevlength(kin, r, θ, r̂, θ̂, ϕ̂, v_t)
                dz = min(λ_s*l_s/R_s, max_dz)
                dz = max(dz, min_dz)
                while !last
                    z_end -= dz
                    if z_end < z_out
                        last = true
                        z_end += dz
                        dz = z_end - z_out
                        z_end -= dz
                        z = z_end + dz/2
                    else
                        z = z_end + dz/2 
                    end
                    # println(z)
                    r = √(x^2 + y^2 + z^2)
                    r̂ = SA[x/r, y/r, z/r]
                    θ = acos(r̂ ⋅ orientation.dipole_axis)
                    r̂, θ̂, ϕ̂ = TTauUtils.Profiles.sphericalbasis(x, y, z, r, θ, orientation)
                    R_to_star_axis = r*sin(acos(r̂ ⋅ orientation.star_axis))
                    sphere_v = TTauUtils.Models.velocityfieldspherical(kin, r, θ)
                    sphere_v += SA[0.0, 0.0, R_to_star_axis*model.star.v_eq*1e5]
                    v_z = -(r̂[3]*sphere_v[1] + θ̂[3]*sphere_v[2] + ϕ̂[3]*sphere_v[3])
                    # v_z = rayvelocity(x, y, z, r, θ, orientation, sphere_v)
                    x_1, x_2 = TTauUtils.spheretogrid(model.geometry, r, θ)
                    T_e = TTauUtils.gridTe(model, x_1, x_2)
                    n_h = TTauUtils.gridnh(model, x_1, x_2)
                    v_t = TTauUtils.Models.hydrogenthermalvelocity(T_e)/√(atomic_mass)
                    if abs(v - v_z*1e-5) < v_t*1e-5
                        # println("$v $(v_z*1e-5) $(v_t*1e-5)")
                        abs_coef = τfunc(T_e, n_h)*1e-4/v_t#*TTauUtils.Models.hydrogenthermalvelocity(T)/1e-4/sqrt(23)
                        τ_ray += dz*abs_coef*R_s
                    end
                    l_s = TTauUtils.Profiles.sobolevlength(model.kinematics, r, θ, r̂, θ̂, ϕ̂, v_t)
                    dz = min(λ_s*l_s/R_s, max_dz)
                    dz = max(dz, min_dz)
                end
                # println(τ_ray[1:20:end])
            end 
            S_tau += dS*(1 - exp(-τ_ray))            
        end
        S_star += dS
    end
    return S_tau/S_star
end

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

function chooseCaNainds(Ca_abs_dip, Na_abs_dip, Ca_thres, Na_thres)
    n_models = length(Ca_abs_dip)
    inds = Int[]
    for i=1:n_models
        if (Ca_abs_dip[i] > Ca_thres) & (Na_abs_dip[i] < Na_thres)
            push!(inds, i)
        end
    end
    return inds
end

Ca_abs_dip_best_stat = calcabsprofiledip(best_stat_pars, best_stat_names, τCa, 40, 200, 0.05, 0.1)
Ca_abs_dip_best_nonstat = calcabsprofiledip(best_nonstat_pars, best_nonstat_names, τCa, 40, 200, 0.05, 0.1)

Na_abs_dip_best_stat = calcabsprofiledip(best_stat_pars, best_stat_names, τNa, 23, 200, 0.05, 0.1)
Na_abs_dip_best_nonstat = calcabsprofiledip(best_nonstat_pars, best_nonstat_names, τNa, 23, 200, 0.05, 0.1)

Ca_nonstat_plt = scatter(best_nonstat_pars[:,2], Ca_abs_dip_best_nonstat)
Ca_stat_plt = scatter(best_stat_pars[:,2], Ca_abs_dip_best_stat)

Na_nonstat_plt = scatter(best_nonstat_pars[:,2], Na_abs_dip_best_nonstat)
Na_stat_plt = scatter(best_stat_pars[:,2], Na_abs_dip_best_stat)

CaNa_best_nonstat_inds = chooseCaNainds(Ca_abs_dip_best_nonstat, Na_abs_dip_best_nonstat, 0.1, 0.1)
CaNa_best_stat_inds = chooseCaNainds(Ca_abs_dip_best_stat, Na_abs_dip_best_stat, 0.1, 0.1)

CaNa_best_nonstat_pars = best_nonstat_pars[CaNa_best_nonstat_inds, :]
CaNa_best_nonstat_names = best_nonstat_names[CaNa_best_nonstat_inds]

CaNa_best_stat_pars = best_stat_pars[CaNa_best_stat_inds, :]
CaNa_best_stat_names = best_stat_names[CaNa_best_stat_inds]

CaNa_gridded_stat_pars, CaNa_gridded_stat_names = putongrid(CaNa_best_stat_pars, CaNa_best_stat_names, (lgṀs, (sigfig = 3, axis = :log)), T_maxs, r_mis, Ws, angs); ""
CaNa_gridded_nonstat_pars, CaNa_gridded_nonstat_names = putongrid(CaNa_best_nonstat_pars, CaNa_best_nonstat_names, (lgṀs, (sigfig = 3, axis = :log)), T_maxs, r_mis, Ws, angs); ""
x = correctnonstatgrid!(CaNa_gridded_nonstat_pars, CaNa_gridded_nonstat_names, CaNa_gridded_stat_pars, CaNa_gridded_stat_names)
x = correctgridforcorotation!(CaNa_gridded_stat_pars, corotationradius(star))
x = correctgridforcorotation!(CaNa_gridded_nonstat_pars, corotationradius(star))



# Ca_plt = plotCa(best_nonstat_pars, best_nonstat_names)
# Na_plt = plotNa(best_nonstat_pars, best_nonstat_names)
# CaNa_plt = plotCaNa(best_stat_pars, best_stat_names)

# ind_T9000 = findmodels(grid_nonstat_pars, [9000], [2], [100])
# ind_T8000 = findmodels(grid_nonstat_pars, [8000], [2], [100])
# ind_T10000 = findmodels(grid_nonstat_pars, [10000], [2], [100])

# best_ind_T9000 = ind_T9000[findmin(grid_nonstat_pars[ind_T9000,9])[2]]
# best_ind_T8000 = ind_T8000[findmin(grid_nonstat_pars[ind_T8000,9])[2]]
# best_ind_T10000 = ind_T10000[findmin(grid_nonstat_pars[ind_T10000,9])[2]]

# plt_nonstat_T9000 = plotmodel(grid_nonstat_pars, grid_nonstat_names, best_ind_T9000)
# plt_stat_T9000 = plotmodel(grid_stat_pars, grid_stat_names, best_ind_T9000)

# plt_nonstat_T8000 = plotmodel(grid_nonstat_pars, grid_nonstat_names, best_ind_T8000)
# plt_stat_T8000 = plotmodel(grid_stat_pars, grid_stat_names, best_ind_T8000)

# plt_nonstat_T10000 = plotmodel(grid_nonstat_pars, grid_nonstat_names, best_ind_T10000)
# plt_stat_T10000 = plotmodel(grid_stat_pars, grid_stat_names, best_ind_T10000)


# plotNa(best_nonstat_pars, best_nonstat_names)
