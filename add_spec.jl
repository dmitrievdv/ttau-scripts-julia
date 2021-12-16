using TTauUtils

function addphotspectomodels(star :: TTauUtils.AbstractStar, suffix)
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
        n_profiles += length(readdir("stars/$star_name/$model_name")) - 1
    end

    # reading models and parameters
    models = []
    n_models = length(model_names)
    n_model = 0
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
            if length(split(profile_name, "_")) > 2
                continue
            end
            profile = HydrogenProfile(star, model, profile_name)
            prof_spec = TTauUtils.addphotosphespecdoppler(profile, 0.1, "spec/M_p5250g4.0z-0.25t1.0_a+0.10c0.00n0.00o+0.10r0.00s0.00_VIS.spec")
            saveprofile(prof_spec, "$(profile_name)_phot2crude")
        end
    end
end

addphotspectomodels(Star("RZPsc"), "stat_nonlocal")
addphotspectomodels(Star("RZPsc"), "nonstat_nonlocal")