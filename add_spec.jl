using SMTPClient
using Distributed
@everywhere using TTauUtils

n_proc = 5

@everywhere function linename(u, l)
    line_name = ""
    if l == 2
        line_name *= 'H'
    end
    if l == 3
        line_name *= "Pa"
    end
    if l == 4
        line_name *= "Br"
    end
    line_name *= 'a' + u-1-l
end

phot_spec = TTauUtils.Profiles.readpolluxspecfile("spec/RZ_Psc_Ha_syn_unwid_corr.dat", 3, 2)

star_name = "RZPsc"
star = Star(star_name)

model_names = readdir("stars/$star_name")
deleteat!(model_names, findall(name -> name == "$star_name.dat", model_names))
deleteat!(model_names, findall(name -> split(name, '_')[end] == "lte", model_names))

n_models = length(model_names)
models = Array{TTauUtils.Models.HydrogenModel, 1}(undef, n_proc)
angs = [35:5.0:60;]


n_iters = n_models ÷ n_proc

u, l = 3, 2

@time for iter_id = 100:n_iters
    n_models_iter = n_proc
    println("$iter_id of $n_iters")
    proc_iter_start = iter_id*n_proc + 1
    proc_iter_end = (iter_id+1)*n_proc
    if proc_iter_end > n_models
        n_models_iter = n_models % n_proc + 1
        proc_iter_end = n_models
    end

    for model_id = proc_iter_start:proc_iter_end
        i = model_id - proc_iter_start + 1
        model = loadmodel(star, model_names[model_id])
        models[i] = model
    end
    # profs = Array{TTauUtils.Profiles.HydrogenProfile, 1}(undef, 0)
    profs = @distributed vcat for model_id = 1:n_models_iter
        model = models[model_id]
        proc_profs = Array{TTauUtils.Profiles.HydrogenProfile, 1}(undef, 0)
        for ang in angs
            int_ang = round(Int, ang)
            profile_name = "$(linename(u,l))_$int_ang"
            # print("stars/$star_name/$(model.name)/$profile_name"*"_phot_vsini.dat")
            if isfile("stars/$star_name/$(model.name)/$profile_name"*"_phot_vsini.dat")
                # println(" file found")
                mv("stars/$star_name/$(model.name)/$profile_name"*"_phot_vsini.dat", "stars/$star_name/$(model.name)/$profile_name"*"_phot-vsini.dat", force = true)
                continue
            end
            try 
                prof = HydrogenProfile(star, model, profile_name)
                prof_phot_vsini = TTauUtils.Profiles.addphotospherespec(prof, 0.1, phot_spec, progress_output = false, sinicorr = true)
                push!(proc_profs, prof_phot_vsini)
            catch
                #
            end
        end
        proc_profs
    end
    for prof in profs
        i_ang = round(Int, prof.orientation.i/π*180)
        profile_name = "$(linename(u,l))_$i_ang"
        saveprofile(prof, profile_name*"_phot-vsini")
        println(prof.model.name, " ", profile_name)
    end
end

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