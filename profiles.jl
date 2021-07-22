using TTauUtils

star_name = "hart94"
star = Star(star_name)
models = readdir("stars/$star_name")[2:end]

transitions = [(3, 2), (4, 2)]
i_angles = [10:10:80;]

log = open("profiles.log", "w")

function linename(u, l)
    line_name = ""
    if l == 2
        line_name *= 'H'
    end
    line_name *= 'a' + u-1-l
end

# models = [models[1], models[2]]
# i_angles = [10, 20]
for model_name in models
    for i_ang in i_angles
        for (u,l) in transitions
            try 
                line_name = linename(u,l)
                println("$model_name, i = $i_ang, line $line_name")
                model = SolidMagnetosphere(star, model_name)
                prof = HydrogenProfileDoppler(model, u, l, i_ang, 0.05, 0.1, 0.1, 200, progress_output = false)
                saveprofile(prof, linename(u,l)*"_$i_ang")
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
