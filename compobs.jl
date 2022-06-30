using TTauUtils
using Dierckx
using Interpolations
using Statistics
using LsqFit
using LaTeXStrings
using SpecialFunctions
using Printf

function calcgridrelation(gridded_pars_1, gridded_pars_2)
    grid_size = size(gridded_pars_1)
    if grid_size != size(gridded_pars_2)
        throw(ErrorException("grid or parameters shapes don't match!"))
    end
    
    gridded_relation = zeros(grid_size)
end

function calcmodel(name, star, Ṁ, T_max, R_in, W; args...) 
    local_id = split(name, '_')[end]
    non_local = local_id == "nonlocal"
    alg_id, local_name = if non_local
        split(name, '_')[end-1], name[1:end-9]
    else
        local_id, name
    end
    model = if alg_id == "stat"
        StationarySolidMagnetosphereNHCool(local_name, star, R_in, R_in + W, Ṁ, T_max, 10; args...)
    elseif alg_id == "nonstat"
        NonStationarySolidMagnetosphereNHCool(local_name, star, R_in, R_in + W, Ṁ, T_max, 10; args...)
    end
    model = if non_local
        addnonlocal(model)
    else
        model
    end
    savemodel(model)
    return model
end

function checkmodel(name, star; star_dir = "stars")
    isfile("$star_dir/$(star.name)/$name/$name.dat")
end

function checkprofile(profile_name, model_name, star; star_dir = "stars")
    isfile("$star_dir/$(star.name)/$model_name/$profile_name.dat")
end

series_name = ["L", "H", "Pa", "Br", "Pf", "Hu"]

linename(u, l) = series_name[l]*('a' - 1 + u-l)

linename(line :: Pair{Int}) = linename(line[1], line[2])

# function calcprofile(star, model_name, profile_name, angle, line; args...)
#     does_model_exist = checkmodel(model_name, star)
#     model = if does_model_exist
#         loadmodel(star, model_name)
#     else
#         throw(ErrorException("model does not exist!"))
#     end
#     profile = HydrogenProfileDoppler(model, line[1], line[2], i, 0.1, 0.1, 0.1, 50)
#     saveprofile(profile, profile_name)
# end

function findmeanvgrad(pars, names; star_dir = "stars", n_θ = 100, n_r = 100)
    n_models = size(pars)[1]
    mean_vgrad = zeros(n_models)
    source_map = TTauUtils.HydrogenPopulations.SourcesMap(100)
    directions = source_map.directions
    solid_angles = source_map.solid_angles
    for i = 1:n_models
        println("$(i-1) from $n_models")
        Ṁ, T_max, R_in, W = pars[i, 1:4]
        lte_name = join(split(names[i][1], '_')[1:3], '_')*"_lte"
        model = if !isfile("$star_dir/$(star.name)/$lte_name/$lte_name.dat")
            model = TTauUtils.Models.SolidMagnetosphereNHCoolLTE(model_name, star, R_in, R_in + W, Ṁ, T_max, 10)
            savemodel(model)
            model
        else
            loadmodel(star, lte_name)
        end
        θ_min = asin(√(1/model.geometry.r_mo))
        θ_step = (π/2 - θ_min)/(n_θ+1)/2
        θ_halfstep = θ_step/2
        v_grad = 0.0
        V = 0.0
        kin = model.kinematics
        for j=1:n_θ
            θ = θ_min + θ_halfstep + θ_step*(j-1)
            r_min = max(1, model.geometry.r_mi*sin(θ)^2)
            r_max = model.geometry.r_mo*sin(θ)^2
            r_step = (r_max - r_min)/(n_r+1)
            r_halfstep = r_step/2
            for k = 1:n_r
                r = r_min + r_halfstep + r_step*(k-1)
                mean_n∇vn = 0.0
                for (direction, dΩ) in zip(directions, solid_angles)
                    α, β = direction
                    n∇vn = abs(TTauUtils.Models.directionalgradient(kin, r, θ, α, β))
                    mean_n∇vn += n∇vn*dΩ/4π
                end
                v_grad += mean_n∇vn*2π*sin(θ)*r^2*θ_step*r_step
                V += 2π*sin(θ)*r^2*θ_step*r_step
            end
        end
        mean_vgrad[i] = v_grad/V
    end
    return mean_vgrad
end

function findmeannh(pars, names; star_dir = "stars", n_θ = 100, n_r = 100)
    n_models = size(pars)[1]
    mean_nh = zeros(n_models)
    for i = 1:n_models
        Ṁ, T_max, R_in, W = pars[i, 1:4]
        lte_name = join(split(names[i][1], '_')[1:3], '_')*"_lte"
        model = if !isfile("$star_dir/$(star.name)/$lte_name/$lte_name.dat")
            model = TTauUtils.Models.SolidMagnetosphereNHCoolLTE(model_name, star, R_in, R_in + W, Ṁ, T_max, 10)
            savemodel(model)
            model
        else
            loadmodel(star, lte_name)
        end
        θ_min = asin(√(1/model.geometry.r_mo))
        θ_step = (π/2 - θ_min)/(n_θ+1)/2
        θ_halfstep = θ_step/2
        NH = 0.0
        V = 0.0
        for j=1:n_θ
            θ = θ_min + θ_halfstep + θ_step*(j-1)
            r_min = max(1, model.geometry.r_mi*sin(θ)^2)
            r_max = model.geometry.r_mo*sin(θ)^2
            r_step = (r_max - r_min)/(n_r+1)
            r_halfstep = r_step/2
            for k = 1:n_r
                r = r_min + r_halfstep + r_step*(k-1)
                nh = TTauUtils.Models.hydrogenconcentration(model, r, θ)
                NH += nh*2π*sin(θ)*r^2*θ_step*r_step
                V += 2π*sin(θ)*r^2*θ_step*r_step
            end
        end
        mean_nh[i] = NH/V
    end
    return mean_nh
end

function findminmaxnh(pars, names; star_dir = "stars")
    n_models = size(pars)[1]
    min_nh = zeros(n_models)
    max_nh = zeros(n_models)
    for i = 1:n_models
        Ṁ, T_max, R_in, W = pars[i, 1:4]
        lte_name = join(split(names[i][1], '_')[1:3], '_')*"_lte"
        model = if !isfile("$star_dir/$(star.name)/$lte_name/$lte_name.dat")
            model = TTauUtils.Models.SolidMagnetosphereNHCoolLTE(model_name, star, R_in, R_in + W, Ṁ, T_max, 10)
            savemodel(model)
            model
        else
            loadmodel(star, lte_name)
        end
        ts = [0:0.01:0.5;]
        R_out = R_in + W
        min_nh[i] = 1e100
        max_nh[i] = 0.0
        for t in ts
            nh = 10^model.lgnh_spl2d(R_out, t)
            if nh < min_nh[i]; min_nh[i] = nh; end
            if nh > max_nh[i]; max_nh[i] = nh; end 
        end
    end
    return min_nh, max_nh
end

function bestmodels(pars, names)
    min_δ = findmin(pars[:,end])[1]
    δ_thres = √2*min_δ
    best = pars[:, end] .< δ_thres
    println(sum(best), " ", length(names))
    return pars[best, :], names[best]
end

function computeltemodels(star, pars, names; star_dir = "stars", args...)
    n_models = size(pars)[1]
    for i = 1:n_models
        Ṁ, T_max, R_in, W = pars[i, 1:4]
        model_name = join(split(names[i][1], '_')[1:3], '_')*"_lte"
        if !isfile("$star_dir/$(star.name)/$model_name/$model_name.dat")
            model = TTauUtils.Models.SolidMagnetosphereNHCoolLTE(model_name, star, R_in, R_in + W, Ṁ, T_max, 10)
            savemodel(model)
        end
    end
end

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

hotspot_gauss_model(x, p) = @. p[1]*exp(-x^2/(2*p[2]^2)) .+ abs(p[3])

function readmodels(star :: TTauUtils.AbstractStar, obs_file, suffix; prof_suffix = "")
    # Assuming there is a grid
    # getting all file names
    star_name = star.name
    model_names = readdir("stars/$star_name")
    
    suffix_words_count = length(split(suffix, "_"))

    # deleting star file
    deleteat!(model_names, findall(name -> name == "$star_name.dat", model_names))

    # clearing nonstat files
    deleteat!(model_names, findall(name -> (join(split(name, '_')[end-suffix_words_count+1:end], '_') != suffix), model_names))
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
    fit_par = [0.0,0.0,0.01]
    names = Vector{String}[]
    for i = 1:n_models
        model_name = model_names[i]
        model = TTauUtils.Models.loadmodel(star, model_name)
        lg10Ṁ = model.Mdot; T_max = model.T_max; r_mi = model.geometry.r_mi; W = model.geometry.r_mo - r_mi
        profile_files = readdir("stars/$star_name/$model_name")
        deleteat!(profile_files, findall(name -> name == "$model_name.dat", profile_files))
        for k = 1:length(profile_files)
            # println("$model_name $star_name ")
            profile_name = profile_files[k][1:end-4]
            if prof_suffix != "" 
                if split(profile_name, "_")[end] != prof_suffix
                    continue
                end
            end
            profile = HydrogenProfile(star, model, profile_name)
            v_mod, r_mod = getvandr(profile)
            r_obs_mod = velocitiesobstomod(v_mod, r_mod, v_obs, r_obs)
            r_to_fit = r_obs_mod .- r_mod 
            # r_mod_obs = velocitiesmodtoobs(v_mod, r_mod, v_obs, r_obs)
            # r_to_fit = r_obs .- r_mod_obs 
            # fit_par = [2.0,50.0,0.00]
            # fit = curve_fit(hotspot_gauss_model, v_obs, r_to_fit, fit_par)
            # fit = curve_fit(hotspot_gauss_model, v_mod, r_to_fit, fit_par) # <- this
            # fit_par = coef(fit) # [0.0,0.0]
            # println(fit_par)
            # res = abs.(r_mod_obs .- r_obs .+ hotspot_gauss_model(v_obs, fit_par))
            res = (r_obs_mod .- r_mod .+ fit_par[3]) .^ 2#) .^ 2 #.- hotspot_gauss_model(v_mod, fit_par)) .^ 2
            # δ = sum(res[abs.(v_obs) .> 0])/length(abs.(v_obs) .> 0)
            δ = sqrt(sum(res[abs.(v_mod) .≥ 30])/length(res[abs.(v_mod) .≥ 30]))
            println("$model_name $star_name $profile_name $δ")
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
    # println(n_model)
    # println(parameters)
    return parameters, names
end

function boundpars(pars, names, bounds...)
    indeces = Int[]
    n_bounds = length(bounds)
    n_models = length(names)
    for i = 1:n_models
        inbounds = true
        cur_par = pars[i, :]
        for j = 1:n_bounds
            par_i = bounds[j][1]
            par_max = bounds[j][3]
            par_min = bounds[j][2]
            par_val = cur_par[par_i]
            if (par_val ≥ par_max) | (par_val ≤ par_min)
                # push!(indeces, i)
                inbounds = false
            end
        end
        if inbounds
            push!(indeces, i)
        end
    end
    return pars[indeces, :], names[indeces]
end

function findnearTmodels(T_max, pars, names, threshold; T_threshold = 10)
    n_models = length(pars[:,1])
    for i=1:n_models
        T_max_i = pars[i,2]
        δ_i = pars[i,end]
        if (abs(T_max_i - T_max) < T_threshold) & (δ_i < threshold)
            println("$i: $(names[i]) $(pars[i,:])")
        end
    end

end

function putongrid(pars :: Matrix{Float64}, names, grid...)
    n_pars = length(pars[1,:])
    n_grid = length(grid)
    println(n_grid)

    grid_sigfig = fill(8, n_grid)
    par_indeces = [1:n_grid;]
    log_axis = fill(false, n_grid)
    pars_length = zeros(Int, n_grid)
    par_axes = fill([], n_grid)
    
    for i = 1:n_grid
        par_grid = grid[i]
        if par_grid isa Array
            par_axes[i] = par_grid
            pars_length[i] = length(par_grid)
        elseif par_grid isa Tuple
            pars_length[i] = length(par_grid[1])
            par_axes[i] = par_grid[1]
            if length(par_grid) > 1
                options = par_grid[2]
                if haskey(options, :index)
                    par_indeces[i] = options.index
                end
                if haskey(options, :sigfig)
                    grid_sigfig[i] = options.sigfig
                end
                if haskey(options, :axis)
                    log_axis[i] = (options.axis == :log)
                end
            end        
        end
    end

    gridded_pars = zeros(n_pars, pars_length...)
    grid_array = zeros(pars_length...)

    for index in keys(grid_array)
        gridded_pars[end, index] = 1e15
        for i in 1:n_grid
            i_ind = index[i]
            par = if log_axis[i]
                10.0 .^ par_axes[i][i_ind]
            else
                par_axes[i][i_ind]
            end
            gridded_pars[i, index] = par
        end
    end

    gridded_names = fill(["", ""], (pars_length...))
    # println(grid_sigfig)
    # println(par_indeces)
    # println(pars_length)
    n_models = length(pars[:,1])
    for i=1:n_models
        cur_pars = pars[i, par_indeces]
        grid_indeces = zeros(Int, n_grid)
        on_grid = true
        for j = 1:n_grid
            par_axis = if log_axis[j]
                10.0 .^ par_axes[j]
            else
                par_axes[j]
            end
            par_value = cur_pars[j]
            exp10 = floor(log10(par_value))
            grid_index_found = findall(x -> abs((x-par_value)/10.0^exp10) < 10.0^(-grid_sigfig[j]+1), par_axis)
            index_found = !isempty(grid_index_found)
            if index_found
                grid_indeces[j] = grid_index_found[1]
            else
                on_grid = false
                break
            end
        end
        if on_grid
            if pars[i, end] < gridded_pars[end, grid_indeces...]
                gridded_pars[:, grid_indeces...] = pars[i,:]
                gridded_names[grid_indeces...] = [names[i][1], names[i][2]]
            end
        end
    end
    return gridded_pars, gridded_names
end

function correctnonstatgrid!(gridded_nonstat_pars, gridded_nonstat_names, gridded_stat_pars, gridded_stat_names)
    indeces = keys(gridded_stat_names)
    for index in indeces
        if gridded_nonstat_pars[end, index] > 1e10
            gridded_nonstat_pars[:, index] .= gridded_stat_pars[:, index]
            gridded_nonstat_names[index] .= gridded_stat_names[index]
        end
    end
end

function correctgridforcorotation!(gridded_pars, r_corr)
    r_mis = gridded_pars[3, 1, 1, :, 1, 1]
    Ws = gridded_pars[4, 1, 1, 1, :, 1]
    r_cor_rounded = r_corr # r_mis[findmin(abs.(r_mis .- r_corr))[2]]
    for i = 1:length(r_mis)
        for j = 1:length(Ws)
            r_mi = r_mis[i]; W = Ws[j]
            # println("$r_mi $W $(r_mi + W > r_cor_rounded) $r_corr $r_cor_rounded")
            if r_mi + W > r_cor_rounded
                gridded_pars[end, :, :, i, j, :] .= 1.0
            end
        end
    end
end

# trunc = pars[δ .< 0.02, :]
# n = size(trunc)[1]
# for i=1:n
#     println(trunc[i,:])
# end
dispersion(xs, x0) = sqrt(sum((xs .- x0) .^ 2)/length(xs))

function getmeananderrors(gridded_pars)
    pars = separatepars(flattenpars(gridded_pars, priority = [1:5;]))
    pars[1] = log10.(pars[1])
    getmeananderrors(gridded_pars, pars)
end

function getmeananderrors(gridded_pars, pars)
    n_p = length(pars)
    means = zeros(n_p+3)
    unierr = zeros(3, n_p)
    err = zeros(3, n_p+3)
    sum_weights = 0.0
    δs = gridded_pars[end, :,:,:,:,:]
    σ = minimum(δs)
    # println(size(δs))
    for index in keys(δs)
        δ = δs[index]
        if δ < √2*σ
            weight = 1/δ
            #weight = exp(-((δ-σ)/σ)^2)
            # δ = 1
            # println("lol")
            sum_weights += weight
            for par_i = 1:n_p
                i = index[par_i]
                cur_par = pars[par_i]
                means[par_i] = means[par_i] + cur_par[i] * weight
            end
            for par_i = n_p+1:n_p+3
                cur_par = gridded_pars[par_i, index[1],index[2],index[3],index[4],index[5]]
                means[par_i] = means[par_i] + cur_par * weight
            end
        end
    end
    means = means / sum_weights
    sum_weights_m = zeros(n_p+3)
    sum_weights_p = zeros(n_p+3)
    for index in keys(δs)
        δ = δs[index]
        if δ < √2*σ
            weight = 1/δ
            # weight = exp(-((δ-σ)/σ)^2)
            # δ = 1
            for par_i = 1:n_p
                cur_par = pars[par_i][index[par_i]]
                if cur_par < means[par_i]
                    err[1, par_i] += (pars[par_i][index[par_i]]-means[par_i])^2 * weight
                    sum_weights_m[par_i] += weight
                else
                    err[3, par_i] += (pars[par_i][index[par_i]]-means[par_i])^2 * weight
                    sum_weights_p[par_i] += weight
                end
                err[2, par_i] += (pars[par_i][index[par_i]]-means[par_i])^2 * weight
            end
            for par_i = n_p+1:n_p+3
                cur_par = gridded_pars[par_i, index[1],index[2],index[3],index[4],index[5]]
                err[2, par_i] += (cur_par-means[par_i])^2 * weight
            end
        end
    end
    err[2,:] = @. √(err[2,:]/sum_weights)
    err[1,:] = @. √(err[1,:]/sum_weights_m)
    err[3,:] = @. √(err[3,:]/sum_weights_p)

    unierr[2,:] = dispersion.(pars, means[1:n_p])

    printmeansanderrors(means, err[2,:], unierr[2,:])
    return means, err[2,:], unierr[2,:]
end

function printmeansanderrors(means, err, unierr)
    parnames = ["lg Ṁ", "T_max", "r_mi", "W", "i", "A", "v", "f"]
    n_p = length(parnames) - 3
    for i=1:n_p
        @printf("%6s = %8.2f ± %6.2f (±%6.2f) [%.2f]\n", parnames[i], means[i], err[i], unierr[i], unierr[i]/err[i])
    end

    # for i=n_p+1:n_p+3
    #     @printf("%6s = %9.3f ± %8.4f \n", parnames[i], means[i], err[i])
    # end
end

function addfluxconstant(f, pars, names)
    new_pars = deepcopy(pars)
    new_pars[:,8] .= f
    for i in 1:length(names)
        model_name, profile_name = names[i]
        if (model_name == "") | (profile_name == ""); continue; end
        model = loadmodel(star, model_name)
        profile = HydrogenProfile(star, model, profile_name)
        v_mod, r_mod = getvandr(profile)
        r_obs_mod = velocitiesobstomod(v_mod, r_mod, v_obs, r_obs)
        res = (r_obs_mod .- (r_mod .+ f)) .^ 2# .- hotspot_gauss_model(v_mod, fit_par)) .^ 2
        δ = sqrt(sum(res[abs.(v_mod) .≥ 0])/length(abs.(v_mod) .≥ 0))
        new_pars[i,9] = δ
    end
    return new_pars
end

function savepars(file, parss, namess)
    open(file*".dat", "w") do io
        model_name = ""
        n_profiles = length(namess)
        for i = 1:n_profiles
            pars = parss[i,:]
            names = namess[i]
            if names[1] != model_name
                model_name = names[1]
                println(io, "Model: $model_name")
                print(io, "model pars: ")
                @printf(io, "%8.6e ", pars[1])
                for par in pars[2:4]
                    @printf(io, "%8.6f ", par)
                end
                print(io, "\n")
            end
            @printf(io, "     profile %s δ = %.6f ", names[2], pars[end])
            for par in pars[5:8]
                @printf(io, "%8.6f ", par)
            end
            print(io, "\n")
        end
    end
end

function findmodels(pars, pars_value, pars_indeces, thresholds)
    n_pars = length(pars_value)
    n_models = length(pars[:, 1])
    indeces = Int[]
    for i = 1:n_models
        index_right = true
        for pid = 1:n_pars
            p = pars_indeces[pid]
            value = pars[i, p]
            if p == 1; value = log10(value); end
            if abs(value - pars_value[pid]) > thresholds[pid]
                index_right = false
                break
            end
        end
        if index_right
            push!(indeces, i)
        end
    end
    return indeces
end

function loadparameters(file, n_model_pars, n_prof_pars)
    n_pars = n_model_pars + n_prof_pars
    pars = Array{Float64,2}(undef, 0, n_pars+1)
    names = Vector{String}[]
    # lines = readlines(file)
    model_pars = zeros(n_model_pars)
    model_name = ""
    profile_name = ""
    profile_pars = zeros(n_prof_pars)
    δ = 1e2
    open(file, "r") do io
    line = "111"
    while line != ""
        line = readline(io)
        # println(line)
        words = split(line, " ", keepempty = false)
        # println(words)
        if length(words) == 0; continue; end
        if words[1] == "Model:"
            model_name = words[end]
            # println(model_name)
            continue
        elseif words[1] == "model"
            # println(words[end - n_model_pars+1:end])
            model_pars = parse.(Float64, words[end - n_model_pars+1:end])
            continue
        elseif words[1] == "profile"
            # println(words[end - n_prof_pars+1:end])
            profile_pars = parse.(Float64, words[end - n_prof_pars+1:end])
            profile_name = words[2]
            δ = parse(Float64, words[end - n_prof_pars])
            # println(model_pars, profile_pars, δ)
            pars = vcat(pars, vcat(model_pars, profile_pars, δ)')
            push!(names, [model_name, profile_name])
        else 
            continue
        end             
    end
    end
    return pars, names
end

function reshapeindeces(grid, priority)
    n_model_pars = length(priority)
    length_model_pars = size(grid)[2:n_model_pars+1]
    n_models = prod(length_model_pars)
    indeces = fill(CartesianIndex(Tuple(1:5)), n_models)
    for i = 1:n_models
        model_index = i-1
        index = [1:5;]
        for par_i in priority
            index[par_i] = model_index % length_model_pars[par_i] + 1
            model_index = model_index ÷ length_model_pars[par_i]
        end
        indeces[i] = CartesianIndex(Tuple(index))
        # println(index)
    end
    indeces
end

function flattenpars(gridded_pars; priority = [5; 1:4], n_pars = 9) 
    indeces = reshapeindeces(gridded_pars, priority)
    pars = gridded_pars[:, indeces]
    reshape(pars, (n_pars, length(pars) ÷ n_pars))'
end

function flattengrid(gridded_pars, gridded_names; priority = [5; 1:4], n_pars = 9)
    indeces = reshapeindeces(gridded_pars, priority)
    pars = gridded_pars[:, indeces]
    names = gridded_names[indeces]
    return reshape(pars, (n_pars, length(pars) ÷ n_pars))', names[:]
end

function separatepars(pars; pars_indeces :: UnitRange{Int} = 1:5)
    separated_pars = Vector{Float64}[]
    n_pars = length(pars_indeces)
    for k = 1:n_pars; push!(separated_pars, Float64[]); end
    cur_pars = zeros(n_pars)
    n_models = size(pars)[1]
    for i = 1:n_models   
        cur_pars = pars[i,pars_indeces]
        for k = 1:n_pars
            cur_par = cur_pars[k]
            found_par_vals = separated_pars[k]
            is_new = true
            for found_par_val in found_par_vals
                if abs(cur_par/found_par_val - 1) < 1e-6
                    is_new = false
                end
            end
            if is_new
                push!(separated_pars[k], cur_pars[k])
            end
        end        
    end
    for k = 1:n_pars
        sort!(separated_pars[k])
    end
    return separated_pars
end


