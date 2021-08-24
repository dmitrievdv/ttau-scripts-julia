using TTauUtils

star_name = "RZPsc"
star = Star(star_name)
savestar(star)

int_10log10_Ṁs = [90:5:110;]
Ṁ_strings = string.(int_10log10_Ṁs)
Ṁs = @. 10^(-int_10log10_Ṁs/10)

# Ṁs = [10^(-7.6), 10^(-7.7), 10^(-7.8), 10^(-7.9), 10^(-8.0), 10^(-8.1), 10^(-8.2), 10^(-8.3), 10^(-8.4)]
# Ṁ_strings = ["76", "77", "78", "79", "80", "81", "82", "83", "84"]

int_T_maxs = [7500, 8000, 8500, 9000, 9500, 10000]
T_max_strings = string.(int_T_maxs)
T_maxs = 1.0*int_T_maxs

mags = [(4.0, 6.0), (4.0, 8.0), (6.0, 8.0), (6.0, 10.0), (8.0, 10.0)]
mag_strings = ["4-6", "4-8", "6-8", "6-10", "8-10"]

log = open("models.log", "w")


for (m, Ṁ) in enumerate(Ṁs)
    Ṁ_str = Ṁ_strings[m]
    for (T_max, T_str) in zip(T_maxs, T_max_strings)
        for ((r_mi, r_mo), mag_str) in zip(mags, mag_strings)
            name = star.name*'_'*Ṁ_str*'_'*T_str*'_'*mag_str
            println(name)
            
            local mag
            local model_name

            try 
                model_name = name*"_stat"
                model_file = "stars/$star_name/$model_name/$model_name.dat"
                if !isfile(model_file)
                    mag = StationarySolidMagnetosphereNHCool(model_name, star, r_mi, r_mo, Ṁ, T_max, 10)
                    savemodel(mag)
                else 
                    println("model exist")
                end
            catch e
                println(log, "Critical failure! $model_name")
                bt = catch_backtrace()
                msg = sprint(showerror, e, bt)
                println(log, msg)
            end

            try
                model_name = name*"_stat_nonlocal"
                model_file = "stars/$star_name/$model_name/$model_name.dat"
                if !isfile(model_file)
                    nonloc_mag = addnonlocal(mag)
                    savemodel(nonloc_mag)
                end
            catch e
                println(log, "Critical failure! $model_name")
                bt = catch_backtrace()
                msg = sprint(showerror, e, bt)
                println(log, msg)
            end

            try
                model_name = name*"_nonstat"
                model_file = "stars/$star_name/$model_name/$model_name.dat"
                if !isfile(model_file)
                    mag = NonStationarySolidMagnetosphereNHCool(model_name, star, r_mi, r_mo, Ṁ, T_max, 10)
                end
                savemodel(mag)
            catch e
                println(log, "Critical failure! $model_name")
                bt = catch_backtrace()
                msg = sprint(showerror, e, bt)
                println(log, msg)
                continue
            end

            try
                model_name = name*"_nonstat_nonlocal"
                model_file = "stars/$star_name/$model_name/$model_name.dat"
                if !isfile(model_file)
                    nonloc_mag = addnonlocal(mag)
                end
                savemodel(nonloc_mag)
            catch e
                println(log, "Critical failure! $model_name")
                bt = catch_backtrace()
                msg = sprint(showerror, e, bt)
                println(log, msg)
            end
        end
    end
end

close(log)
