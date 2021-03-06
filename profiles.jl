using TTauUtils

for star_name in ["hart94"]
# star_name = "coldhart94"
star = Star(star_name)
models = readdir("stars/$star_name")[2:end]#["hart94_85_8500_4-6_stat_nonlocal", "hart94_85_8500_4-6_nonstat_nonlocal"]
models = ["hart94_90_8000_2-3_stat_nonlocal", "hart94_90_8000_2-3_nonstat_nonlocal"]

transitions = [(3, 2)]
i_angles = [0:30:90;]

log = open("profiles.log", "w")

function linename(u, l)
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

# models = [models[1], models[2]]
# i_angles = [10, 20]
forced = true
for model_name in models
    for i_ang in i_angles
        for (u,l) in transitions
            try 
                line_name = linename(u,l)
                println("$model_name, i = $i_ang, line $line_name")
                profile_file = "stars/$star_name/$model_name/$(line_name)_$i_ang.dat"
                if isfile(profile_file) && !forced
                    println("Profile exist already.")
                    continue
                end    
                model = SolidMagnetosphere(star, model_name)
                # println(split(model_name, '_')[end])
                if split(model_name, '_')[end] == "nonlocal"
                    @time prof = HydrogenProfileDoppler(model, u, l, i_ang, 0.05, 0.1, 0.1, 200, progress_output = false)
                    saveprofile(prof, linename(u,l)*"_$i_ang")
                end
            catch e
                println(log, model_name)
                bt = catch_backtrace()
                msg = sprint(showerror, e, bt)
                println(log, msg)
                continue
            end
        end
    end
end

close(log)
end