using TTauUtils
using SMTPClient

function addphotspectomodels(star :: TTauUtils.AbstractStar, suffix, prof_suffix)
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
        n_prof = 0
        deleteat!(profile_files, findall(name -> name == "$model_name.dat", profile_files))
        println(model_name)
        for k = 1:length(profile_files)
            profile_name = profile_files[k][1:end-4]
            prof_with_spec = (length(split(profile_name, "_")) > 2)
            spec_already_added = isfile("stars/$star_name/$model_name/$(profile_name)_$(prof_suffix).dat")
            if !(prof_with_spec | spec_already_added)
                n_prof += 1
                profile = HydrogenProfile(star, model, profile_name)
                println(profile_name)
                prof_spec = TTauUtils.addphotosphespecdoppler(profile, 0.1, "spec/RZ_Psc_Ha_syn_unwid_corr.dat", progress_output = false)
                saveprofile(prof_spec, "$(profile_name)_$prof_suffix")
                prof_spec = TTauUtils.addphotosphespecdoppler(profile, 0.1, "spec/RZ_Psc_Ha_syn_unwid_corr.dat", sinicorr = true, progress_output = false)
                saveprofile(prof_spec, "$(profile_name)_$(prof_suffix)-sini")
            end
        end
    end
end

addphotspectomodels(Star("RZPsc"), "stat_nonlocal", "phot3crude")
addphotspectomodels(Star("RZPsc"), "nonstat_nonlocal", "phot3crude")

opt = SendOptions(
  isSSL = true,
  username = "my.calc.results@gmail.com",
  passwd = "vuN2jtpUXgApc6t")
#Provide the message body as RFC5322 within an IO
body = IOBuffer("Photospheric profile added")
url = "smtps://smtp.gmail.com:465"
rcpt = ["<dmitrievdv242@gmail.com>"]
from = "<my.calc.results@gmail.com>"
resp = send(url, rcpt, from, body, opt)