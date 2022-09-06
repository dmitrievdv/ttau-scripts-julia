using TTauUtils

dec_char = ["d", "c"]

function rstring(r; max_d = 2)
    r_int = max_d
    while (round(Int, r*10^(r_int*1.0)) % 10) == 0
        r_int -= 1
    end
    r_int
    whole = floor(Int, r)
    dec_str = if r_int <= 0
        string(whole)
    else 
        part = floor(Int, r*10^r_int) - whole*10^r_int
        zeros = r_int - floor(Int, log10(part)) - 1
        dec_char[r_int]*string(whole)*("0"^zeros)*string(part)
    end
    return dec_str
end

function modelname(Ṁ, T_max, r_mi, r_mo)
    mag_string = "$(rstring(r_mi))-$(rstring(r_mo))"
    T_string = string(round(Int, T_max))
    lg_Ṁ = log10(Ṁ)
    Ṁ_string = string(round(Int, -lg_Ṁ*10))
    return "$(Ṁ_string)_$(T_string)_$(mag_string)"
end

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

algsuffix = Dict("HartNHCoolStat" => "stat_nonlocal", 
                 "HartNHCoolNonStat" => "nonstat_nonlocal", 
                 "HartNHCoolLTE" => "lte")

r_mis = [2.0:1:10.0;]
Ws = [1:0.2:4;]
T_maxs = [7000:1000:15000;]
lgṀs = [-11:0.2:-8.4;]
angs = [35:5:60;]


alg = "HartNHCoolStat"
suffix = algsuffix[alg]

prof_suffix = "phot3crude"
u = 3; l = 2

line = linename(u, l)

star = Star("RZPsc")

star_dir = "stars/$(star.name)"

log = open("checmodels.log", "w")
problem = open("checmodels.dat", "w")

for lgṀ in lgṀs, T_max in T_maxs, r_mi in r_mis, W in Ws
    r_mo = r_mi + W
    if r_mo > TTauUtils.corotationradius(star); continue; end
    model_name = modelname(10^lgṀ, T_max, r_mi, r_mo) * "_" * suffix
    model_dir = "$star_dir/$model_name"
    model_file_exist = isfile("$model_dir/$model_name.dat")
    if !model_file_exist
        if alg != "HartNHCoolNonStat"
            println(log, "Model $model_name does not exist!")
            println(problem, model_name)
        end
        continue
    end
    model = TTauUtils.loadmodel(star, model_name)
    parameters_ok = ((abs(log10(model.Mdot) - lgṀ) < 5e-3) & 
                     (abs(model.T_max - T_max) < 1) & 
                     (abs(model.geometry.r_mi - r_mi) < 1e-2) &
                     (abs(model.geometry.r_mo - r_mo) < 1e-2))
    if !parameters_ok
        println(log, "Model $model_name has wrong parameters!")
        println(problem, model_name)
        continue
    end
    if model.alg != alg
        println(log, "Wrong algorithm in $(model_name)!")
        println(problem, model_name)
        continue
    end
    for ang in angs
        main_profile_name = line * "_" * string(round(Int, ang))
        main_profile_file_exist = isfile("$model_dir/$main_profile_name.dat")
        if !main_profile_file_exist
            println(log, "Profile $main_profile_name is not computed for model $(model_name)!")
            println(problem, model_name)
            continue
        end
        profile = HydrogenProfile(star, model, main_profile_name)
        prof_ang = profile.orientation.i/π*180
        parameters_ok = (u == profile.upper_level) & (l == profile.lower_level) & (abs(prof_ang - ang) < 0.1)
        if !parameters_ok
            println(log, "Profile $main_profile_name has wrong parameters (in $(model_name))!")
            println(problem, model_name)
            continue
        end
        phot_profile_name = main_profile_name * "_" * prof_suffix
        phot_profile_file_exist = isfile("$model_dir/$phot_profile_name.dat")
        if !phot_profile_file_exist
            println(log, "Profile $main_profile_name is not corrected for photosphere profile for model $(model_name)!")
            println(problem, model_name)
            continue
        end
    end
end

close(log)
close(problem)
