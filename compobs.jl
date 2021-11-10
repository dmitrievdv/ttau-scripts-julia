using TTauUtils
using Dierckx
using Interpolations

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

function readmodelgrid(star :: TTauUtils.AbstractStar, obs_file, deviationfunction)
    # Assuming there is a grid
    # getting all file names
    star_name = star.name
    model_names = readdir("stars/$star_name")
    
    # deleting star file
    deleteat!(model_names, findall(name -> name == "RZPsc.dat", model_names))

    # clearing nonstat files
    deleteat!(model_names, findall(name -> (split(name, '_')[end-1] == "nonstat"), model_names))
    n_models = length(model_names)
    
    # counting profiles
    n_profiles = 0
    for model_name in model_names
        n_profiles += length(readdir("stars/$star_name/$model_name")) - 1
    end

    #reading observation profile
    v_obs, r_obs = readobservation(obs_file)

    # reading models and parameters
    models = []
    n_models = length(model_names)
    parameters = zeros(n_profiles, 6) # six parameters: lg Ṁ, T_max, r_mi, W (r_mo = r_mi + W), i
    n_model = 0
    for i = 1:n_models
        model_name = model_names[i]
        model = TTauUtils.Models.loadmodel(star, model_name)
        lg10Ṁ = model.Mdot; T_max = model.T_max; r_mi = model.geometry.r_mi; W = model.geometry.r_mo - r_mi
        profile_files = readdir("stars/$star_name/$model_name")
        deleteat!(profile_files, findall(name -> name == "$model_name.dat", profile_files))
        for k = 1:length(profile_files)
            profile_name = profile_files[k][1:end-4]
            profile = HydrogenProfile(star, model, profile_name)
            v_mod, r_mod = getvandr(profile)
            δ = deviationfunction(v_obs, r_obs, v_mod, r_mod)
            i_ang = profile.orientation.i
            n_model += 1
            parameters[n_model, :] .= [lg10Ṁ, T_max, r_mi, W, i_ang/π*180, δ]
        end
        println(model_name)
    end
    return parameters
end

v_obs, r_obs = readobservation("observation.dat")
star = Star("RZPsc")
pars = readmodelgrid(star, "observation.dat", absorbtionmeanabs)
δ = pars[:,6]
trunc = pars[δ .< 0.019, :]
n = size(trunc)[1]
for i=1:n
    println(trunc[i,:])
end