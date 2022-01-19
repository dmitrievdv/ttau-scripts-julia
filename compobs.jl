using TTauUtils
using Dierckx
using Interpolations
using Statistics
using LsqFit
using Plots
using LaTeXStrings
using SpecialFunctions

function readobservation(filename :: AbstractString)
    r_obs = Float64[]
    v_obs = Float64[]
    for line in readlines(filename)
        if !isempty(line)
            push!.([v_obs, r_obs], parse.(Float64, split(line)))
        end
    end
    v_obs, r_obs
end

function absorbtionmeansquares(v_obs, r_obs, v_mod, r_mod) # Assuming v_mod and v_obs are sorted
    n_points = 0
    squares = 0.0
    model = interpolate((v_mod,), r_mod, Gridded(Linear()))
    for (i,v) in enumerate(v_obs)
        if r_obs[i] < 1.0
            squares += (model(v) - r_obs[i])^2
            n_points += 1
        end
    end
    return squares/n_points
end

function absorbtionmeanabs(v_obs, r_obs, v_mod, r_mod) # Assuming v_mod and v_obs are sorted
    n_points = 0
    squares = 0.0
    model = interpolate((v_mod,), r_mod, Gridded(Linear()))
    for (i,v) in enumerate(v_obs)
        if r_obs[i] < 1.0
            squares += abs(model(v) - r_obs[i])
            n_points += 1
        end
    end
    return squares/n_points
end

function velocitiesmodtoobs(v_obs, r_obs, v_mod, r_mod) # Assuming v_mod and v_obs are sorted
    n_points = 0
    n_obs = length(v_obs)
    squares = 0.0
    r_mod_obs = zeros(n_obs)
    model = interpolate((v_mod,), r_mod, Gridded(Linear()))
    v_min = minimum(v_mod)
    v_max = maximum(v_mod)
    for (i,v) in enumerate(v_obs)
        r_mod_obs[i] = if (v < v_min) | (v > v_max)
            r_obs[i]
        else
            model(v)
        end
    end
    return r_mod_obs
end

function velocitiesobstomod(v_mod, r_mod, v_obs, r_obs) # Assuming v_mod and v_obs are sorted
    r_obs_mod = velocitiesmodtoobs(v_mod, r_mod, v_obs, r_obs)
    return r_obs_mod
end

function getvandr(profile :: TTauUtils.HydrogenProfile)
    u = profile.upper_level; l = profile.lower_level
    ν_0 = TTauUtils.HydrogenPopulations.linefrequency(u,l)
    v_prof = (ν_0 .- profile.frequencies) / ν_0*3e5 
    r_prof = profile.profile
    if v_prof[1] > 0.0
        reverse!(v_prof)
        reverse!(r_prof)
    end
    return v_prof, r_prof
end

function voigt(x, σ, γ)
    z = (γ - 1im*x)/(√2*σ)
    real(erfcx(z))/(σ*√(2π))
end

function hotspot_voigt_model(x, p)
    @. voigt(x, p[2], p[3])*p[1]/voigt(0, p[2], p[3]) .+ 0.01
end

hotspot_gauss_model(x, p) = @. p[1]*exp(-x^2/(2*p[2]^2)) .+ 0.01

function readmodels(star :: TTauUtils.AbstractStar, obs_file, suffix; prof_suffix = "")
    # Assuming there is a grid
    # getting all file names
    star_name = star.name
    model_names = readdir("stars/$star_name")
    
    # deleting star file
    deleteat!(model_names, findall(name -> name == "RZPsc.dat", model_names))

    # clearing nonstat files
    deleteat!(model_names, findall(name -> (join(split(name, '_')[4:end], '_') != suffix), model_names))
    n_models = length(model_names)
    
    # counting profiles
    n_profiles = 0
    for model_name in model_names
        for profile_file in readdir("stars/$star_name/$model_name")
            profile_name = profile_file[1:end-4]
            if profile_name != model_name
                if prof_suffix != "" 
                    if split(profile_name, "_")[end] == prof_suffix
                        n_profiles += 1 
                    end
                else
                    n_profiles += 1 
                end
            end
        end
        
    end

    #reading observation profile
    v_obs, r_obs = readobservation(obs_file)

    # reading models and parameters
    models = []
    n_models = length(model_names)
    parameters = zeros(n_profiles, 9) # five parameters and error: lg Ṁ, T_max, r_mi, W (r_mo = r_mi + W), i, gauss spot, δ
    n_model = 0
    # gaussian
    # hotspot_gauss_model(x, p) = @. p[1]*exp(-x^2/(2*p[2]^2))
    fit_par = [0.0,0.0,0.0]
    names = Vector{String}[]
    for i = 1:n_models
        model_name = model_names[i]
        model = TTauUtils.Models.loadmodel(star, model_name)
        lg10Ṁ = model.Mdot; T_max = model.T_max; r_mi = model.geometry.r_mi; W = model.geometry.r_mo - r_mi
        profile_files = readdir("stars/$star_name/$model_name")
        deleteat!(profile_files, findall(name -> name == "$model_name.dat", profile_files))
        for k = 1:length(profile_files)
            println(model_name)
            profile_name = profile_files[k][1:end-4]
            if prof_suffix != "" 
                if split(profile_name, "_")[end] != prof_suffix
                    continue
                end
            end
            profile = HydrogenProfile(star, model, profile_name)
            v_mod, r_mod = getvandr(profile)
            r_mod_obs = velocitiesobstomod(v_obs, r_obs, v_mod, r_mod)
            r_to_fit = r_obs .- r_mod_obs 
            fit_par = [0.6,50.0,10.0]
            fit = curve_fit(hotspot_gauss_model, v_obs, r_to_fit, fit_par)
            fit_par = coef(fit) # [0.0,0.0]
            println(fit_par)
            res = abs.(r_mod_obs .- r_obs .+ hotspot_gauss_model(v_obs, fit_par))
            δ = sum(res[abs.(v_obs) .> 0])/length(abs.(v_obs) .> 0)
            i_ang = profile.orientation.i
            n_model += 1
            if n_model > n_profiles
                break
            end
            parameters[n_model, :] .= [lg10Ṁ, T_max, r_mi, W, i_ang/π*180, fit_par[1], fit_par[2], fit_par[3], δ]
            push!(names, [model_name, profile_name])
        end
    end

    # deleting not-on-grid points
    # n_par = length(parameters[1,1:end-1])
    # similar_points = zeros(Int, n_profiles)
    # similar_points_occurences = Int[]
    # similar_points_values = Int[]
    # for i=1:n_profiles
    #     similar_pars = zeros(n_par)
    #     for j=1:n_profiles
    #         if j == i
    #             similar_points[i] += 1
    #         else
    #             for k=1:n_par
    #                 if abs(parameters[i, k] - parameters[j, k])/parameters[i,k] < 1e-5
    #                     similar_points[i] += 1
    #                     similar_pars[k] += 1
    #                     break
    #                 end
    #             end
    #         end
    #     end
    #     loc = findfirst(n -> n == similar_points[i], similar_points_values)
    #     if isnothing(loc)
    #         push!(similar_points_occurences, 1)
    #         push!(similar_points_values, similar_points[i])
    #     else
    #         similar_points_occurences[loc] += 1
    #     end
    # end

    # for (val, occ) in zip(similar_points_values, similar_points_occurences)
    #     println(val, " ", occ)
    # end

    # # similar_points_summed = dropdims(sum(similar_points, dims = 2), dims = 2)
    # # n_grid = median(sum(similar_points_summed))

    # deleteat!(model_names, findall(n -> n < n_grid, similar_points_summed))
    # println(n_grid, " ", length(model_names))

    return parameters, names
end

v_obs, r_obs = readobservation("spec/RZPsc_16-11-2013_proc.dat")
star = Star("RZPsc")
stat_pars, stat_names = readmodels(star, "spec/RZPsc_16-11-2013_proc.dat", "stat_nonlocal", prof_suffix = "phot3crude")
nonstat_pars, nonstat_names = readmodels(star, "spec/RZPsc_16-11-2013_proc.dat", "nonstat_nonlocal", prof_suffix = "phot3crude")
# δs = pars[:,8]
min_stat_δ = minimum(stat_pars[:,end])
# min_nonstat_δ = minimum(nonstat_pars[:,end])

# pars = vcat(nonstat_pars, stat_pars)
# names = vcat(nonstat_names, stat_names)

function plotmodel(pars, names, id)
    plt = plot()
    model_name, profile_name = names[id]
    println(names[id], " ", pars[id,5:end])
    model = TTauUtils.Models.loadmodel(star, model_name)
    profile_nophot = HydrogenProfile(star, model, profile_name[1:5])
    profile = HydrogenProfile(star, model, profile_name)
    v_mod, r_mod = getvandr(profile)
    v_mag, r_mag = getvandr(profile_nophot)
    r_gauss_mod = r_mod .+ hotspot_gauss_model(v_mod, pars[id,6:end-1])
    plot!(plt, v_mod, r_gauss_mod, lc = :red, label = "all")
    plot!(plt, v_mod, hotspot_gauss_model(v_mod, pars[id,6:end-1]) .+ 1, lc = :red, la = 0.5, ls =:dash, label = "hotspot")
    plot!(plt, v_mod, r_mod, lc = :orange, ls = :dash, la = 0.5, label = "mag + phot")
    plot!(plt, v_mag, r_mag, lc = :blue, ls = :dash, la = 0.5, label = "mag")
    plot!(plt, v_obs, r_obs, xlims = (-500, 600), ylims = (0.4, 1.4), lc = :black, label = "obs")
    r_obs_mod = velocitiesobstomod(v_mod, r_gauss_mod, v_obs, r_obs)
    plot!(plt, v_mod, r_gauss_mod .- r_obs_mod .+ 1, lc = :black, la = 0.3, label = "mod - obs")
end

function plotmodelerr(pars, names, id)
    plt = plot()
    model_name, profile_name = names[id]
    println(names[id], " ", pars[id,5:end])
    model = TTauUtils.Models.loadmodel(star, model_name)
    profile_nophot = HydrogenProfile(star, model, profile_name[1:5])
    profile = HydrogenProfile(star, model, profile_name)
    v_mod, r_mod = getvandr(profile)
    v_mag, r_mag = getvandr(profile_nophot)
    r_gauss_mod = r_mod .+ hotspot_gauss_model(v_mod, pars[id,6:end-1])
    r_obs_mod = velocitiesobstomod(v_mod, r_gauss_mod, v_obs, r_obs)
    plot!(plt, v_mod, r_gauss_mod .- r_obs_mod, lc = :red, label = "mod - obs")
end

function plotδ(pars, names, δ_cut; la = 0.1)
    
    plt = plot(legend = false)
    # min_δ = min(min_stat_δ, min_nonstat_δ)
    best_names = []
    best_pars = []
    for i=1:length(names)
        model_name, profile_name = names[i]
        type = split(model_name, '_')[4]
        # if type == "stat"
        #     min_δ = min_stat_δ
        # elseif type == "nonstat"
        #     min_δ = min_nonstat_δ
        # end
        if pars[i,end] < δ_cut
            println(names[i], " ", pars[i,5:end])
            model = TTauUtils.Models.loadmodel(star, model_name)
            profile_nophot = HydrogenProfile(star, model, profile_name[1:5])
            profile = HydrogenProfile(star, model, profile_name)
            v_mod, r_mod = getvandr(profile)
            v_mag, r_mag = getvandr(profile_nophot)
            r_gauss_mod = r_mod .+ hotspot_gauss_model(v_mod, pars[i,6:end-1])
            plot!(plt, v_mod, r_gauss_mod, lc = :red, la = la)
            plot!(plt, v_mod, hotspot_gauss_model(v_mod, pars[i,6:end-1]) .+ 1, lc = :red, ls =:dash, la = la)
            plot!(plt, v_mod, r_mod, lc = :orange, ls = :dash, la = la)
            plot!(plt, v_mag, r_mag, lc = :blue, ls = :dash, la = la)
            # push!(plots, plt)
            push!(best_names, [model_name, profile_name])
            push!(best_pars, [i; pars[i,:]])
        end
    end
    plot!(plt, v_obs, r_obs, xlims = (-500, 600), ylims = (0.4, 1.2), lc = :black)

    for (name, pars) in zip(best_names, best_pars)
        println(name)
        println(pars)
    end
    plt
end

function putongrid(lgṀs, T_maxs, r_mis, Ws, angs, pars, names)
    n_models = length(pars[:,1])
    n_Ṁ = length(lgṀs)
    n_T = length(T_maxs)
    n_rmi = length(r_mis)
    n_W = length(Ws)
    n_ang = length(angs)
    n_pars = length(pars[1,:])
    gridded_pars = fill(0.0, (n_pars, n_Ṁ, n_T, n_rmi, n_W, n_ang))
    gridded_names = fill(["",""], n_Ṁ, n_T, n_rmi, n_W, n_ang)
    for i=1:n_models
        Ṁ, T_max, r_mi, W, ang = pars[i,1:5]
        lgṀ = log10(Ṁ)
        # println(pars[i, 1:5], " ", names[i])
        i_Ṁ = findall(x -> abs(x-lgṀ) < 1e-4, lgṀs)
        i_T = findall(x -> abs(x-T_max) < 1e-4, T_maxs)
        i_rmi = findall(x -> abs(x-r_mi) < 1e-4, r_mis)
        i_W = findall(x -> abs(x-W) < 1e-4, Ws)
        i_ang = findall(x -> abs(x-ang) < 1e-4, angs)
        try
            gridded_pars[:, i_Ṁ[1], i_T[1], i_rmi[1], i_W[1], i_ang[1]] = pars[i,:]
            gridded_names[i_Ṁ[1], i_T[1], i_rmi[1], i_W[1], i_ang[1]] = [names[i][1], names[i][2]]
        catch err
            println("Model $(names[i]) not on grid: $(pars[i, 1:5])")
            println("$i_Ṁ $i_T $i_rmi $i_W $i_ang")
        end
    end
    return gridded_pars, gridded_names
end



lgṀs = [-11:0.1:-9;]
T_maxs = [8e3:500:12000;]
Ws = [1.0:1:5;]
r_mis = [5.0:1:10;]
angs = [40:2.0:60;]

gridded_stat_pars, gridded_stat_names = putongrid(lgṀs, T_maxs, r_mis, Ws, angs, stat_pars, stat_names); ""
gridded_nonstat_pars, gridded_nonstat_names = putongrid(lgṀs, T_maxs, r_mis, Ws, angs, nonstat_pars, nonstat_names); ""


function plotδheatTM(gridded_pars, lgṀs, T_maxs, i_rm, i_W, i_ang, δ_cut)
    ang = angs[i_ang]
    r_mi = r_mis[i_rm]
    W = Ws[i_W]
    r_mo = r_mi + W
    title = L"i = %$(ang)\degree,\ r_\mathrm{mi} = %$(r_mi),\ r_\mathrm{mo} = %$(r_mo)"
    heatmap(T_maxs, lgṀs, log10.(gridded_pars[9,:,:,i_rm,i_W,i_ang]), title = title, clims = (-2, -1.2))
end

function plotδheat(gridded_pars, dims, xs, ys; clims = (0, 0.1))
    δ_i = size(gridded_pars)[1]
    griddedδ = selectdim(gridded_pars, 1, δ_i)
    hm_size = size(griddedδ, dims[1]), size(griddedδ, dims[2])
    n_1, n_2 = hm_size
    mindims = deleteat!([1:ndims(griddedδ);], dims)
    println(mindims)
    println(hm_size)
    hm = reshape(minimum(griddedδ, dims = mindims), hm_size)
    heatmap(ys, xs, hm, clims = clims)
end

function plotδfromgrid(gridded_pars, gridded_names, i_rm, i_W, i_ang, δ_cut)
    local names = vec(gridded_names[:,:,i_rm, i_W, i_ang])
    local pars = transpose(reshape(vec(gridded_pars[:,:,:,i_rm, i_W, i_ang]), (9, length(vec(gridded_pars[1,:,:,i_rm, i_W, i_ang])))))
    plotδ(pars, names, 0.02)
end

function addHb(star, model_name, profile_name; sinicorr = true)
    model = loadmodel(star, model_name)
    prof_Ha = HydrogenProfile(star, model, profile_name)
    prof_Ha_nos = HydrogenProfile(star, model, profile_name[1:5])
    prof_Ha = TTauUtils.addphotosphespecdoppler(prof_Ha_nos, 0.05, "spec/M_p5250g4.0z-5.00t1.0_a+0.40c0.00n0.00o+0.40r0.00s0.00_VIS.spec", sinicorr = sinicorr, Δλ = 1.98)
    sinicorr = occursin("sini", profile_name)
    i = parse(Float64, split(profile_name, "_")[2])
    prof_Hb_nos = HydrogenProfileDoppler(model, 4, 2, 1.0*i, 0.05, 0.05, 0.1, 200)
    prof_Hb = TTauUtils.addphotosphespecdoppler(prof_Hb_nos, 0.05, "spec/M_p5250g4.0z-5.00t1.0_a+0.40c0.00n0.00o+0.40r0.00s0.00_VIS.spec", sinicorr = sinicorr, Δλ = 1.4)
    Hb_name = "Hb"*profile_name[3:end]
    saveprofile(prof_Hb, Hb_name)
    v_Ha_mod, r_Ha_mod = getvandr(prof_Ha)
    v_Hb_mod, r_Hb_mod = getvandr(prof_Hb)
    v_Ha_nos, r_Ha_nos = getvandr(prof_Ha_nos)
    v_Hb_nos, r_Hb_nos = getvandr(prof_Hb_nos)
    v_Ha_obs, r_Ha_obs = readobservation("spec/RZPsc_16-11-2013_proc.dat")
    v_Hb_obs, r_Hb_obs = readobservation("spec/RZPsc_Hb_16-11-2013_proc.dat")
    Ha_plt = plot(v_Ha_obs, r_Ha_obs, legend = false)
    plot!(Ha_plt, v_Ha_mod, r_Ha_mod .+ 0.01)
    plot!(Ha_plt, v_Ha_nos, r_Ha_nos .+ 0.01)
    Hb_plt = plot(v_Hb_obs, r_Hb_obs, legend = false)
    plot!(Hb_plt, v_Hb_mod, r_Hb_mod .- 0.03)
    plot!(Hb_plt, v_Hb_nos, r_Hb_nos .- 0.03)
    plot(Ha_plt, Hb_plt, layout = @layout([A B]), size = (1400, 700), dpi = 600)
end

function plotgauss(pars, names, δ_cut; la = 0.1)
    
    plt = plot(clims = (-11,-9))
    # min_δ = min(min_stat_δ, min_nonstat_δ)
    best_names = []
    best_pars = []
    min_δ = minimum(pars[:,end])
    δ_range = δ_cut - min_δ
    alphacolor(δ) = (1.0 - (δ - min_δ)/δ_range)
    for i=1:length(names)
        model_name, profile_name = names[i]
        type = split(model_name, '_')[4]
        # if type == "stat"
        #     min_δ = min_stat_δ
        # elseif type == "nonstat"
        #     min_δ = min_nonstat_δ
        # end
        if pars[i,end] < δ_cut
            # println(names[i], " ", pars[i,5:end])
            # push!(plots, plt)
            scatter!(plt, [pars[i,6]], [pars[i,7]], marker_z=log10.(pars[i,1]), ma = (alphacolor(pars[i,end]))^(2), label = false)
            push!(best_names, [model_name, profile_name])
            push!(best_pars, [i; pars[i,:]])
        end
    end
    # best_pars = hcat(best_pars...)
    # scatter!(plt, best_pars[7,:], best_pars[8,:], marker_z=log10.(best_pars[2,:]))
    # show(best_pars)
    # for (name, pars) in zip(best_names, best_pars)
    #     println(name)
    #     println(pars)
    # end
    plt
end
# trunc = pars[δ .< 0.02, :]
# n = size(trunc)[1]
# for i=1:n
#     println(trunc[i,:])
# end