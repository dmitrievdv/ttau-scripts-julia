using Distributed
using SMTPClient

n_proc = 6
addprocs(n_proc)

@everywhere using TTauUtils

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

@everywhere function modelname(Ṁ, T_max, r_mi, r_mo)
    mag_string = "$(rstring(r_mi))-$(rstring(r_mo))"
    T_string = string(round(Int, T_max))
    lg_Ṁ = log10(Ṁ)
    Ṁ_string = string(round(Int, -lg_Ṁ*10))
    return "$(Ṁ_string)_$(T_string)_$(mag_string)"
end

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

star_name = "SCrA_SE"
star = Star(star_name)

# PDS 70
# r_mis = [2.0:0.5:6.0;]
# mag_widths = [0.1:0.1:2;]
# T_maxs = [10000:1000:15000;]
# lg10_Ṁs = [-12:0.2:-9;]
# i_angs = [35:5:65;]

# RZ Psc
# r_mis = [2.0:1:10.0;]
# mag_widths = [1:0.2:4;]
# T_maxs = [7000,8000]#:1000:15000;]
# lg10_Ṁs = [-10:0.2:-8.4;]
# i_angs = [35:5:60;]

# SCrA_SE
r_mis = [2.0:1:6.0;]
mag_widths = [0.5:0.5:4;]
T_maxs = [6000:1000:10000;]#:1000:15000;]
lg10_Ṁs = [-8:0.2:-6;]
i_angs = [40:5:80;]

u = 5; l = 2

parameters = Vector{Float64}[]

stats_exist = Bool[]
nonstats_exist = Bool[]

stats_ignore = Bool[]
nonstats_ignore = Bool[]

model_names = String[]

stats_angles = Vector{Float64}[]
nonstats_angles = Vector{Float64}[]

n_models = 0

for r_mi in r_mis, mag_width in mag_widths, T_max in T_maxs, lg10_Ṁ in lg10_Ṁs
    r_mo = r_mi + mag_width
    if r_mo > ceil(TTauUtils.Stars.corotationradius(star)); continue; end
    mag_string = "$(round(Int, r_mi))-$(round(Int, r_mo))"
    T_string = string(round(Int, T_max))
    Ṁ_string = string(round(Int, -lg10_Ṁ*10))
    Ṁ = 10.0^lg10_Ṁ
    model_name = modelname(Ṁ, T_max, r_mi, r_mo)
    rm("stars/$star_name/$(model_name)_lte", force = true, recursive = true)
    # stat_exist = false; rm("stars/$star_name/$(model_name)_stat_nonlocal", force = true, recursive = true) # isfile("stars/$star_name/$(model_name)_stat_nonlocal/$(model_name)_stat_nonlocal.dat")
    stat_exist = isfile("stars/$star_name/$(model_name)_stat_nonlocal/$(model_name)_stat_nonlocal.dat")
    # nonstat_exist = false; rm("stars/$star_name/$(model_name)_nonstat_nonlocal", force = true, recursive = true) # isfile("stars/$star_name/$(model_name)_nonstat_nonlocal/$(model_name)_nonstat_nonlocal.dat")
    nonstat_exist = isfile("stars/$star_name/$(model_name)_nonstat_nonlocal/$(model_name)_nonstat_nonlocal.dat")
    stat_angles = Float64[]
    nonstat_angles = Float64[]
    stat_ignore = true; nonstat_ignore = true
    for i_ang in i_angs
        profile_name = "$(linename(u,l))_$i_ang"
        stat_profile_exist = isfile("stars/$star_name/$(model_name)_stat_nonlocal/$profile_name.dat")
        nonstat_profile_exist = isfile("stars/$star_name/$(model_name)_nonstat_nonlocal/$profile_name.dat")
        if !stat_profile_exist
            stat_ignore = false
            push!(stat_angles, i_ang)
        end
        if !nonstat_profile_exist
            nonstat_ignore = false
            push!(nonstat_angles, i_ang)
        end
    end

    

    if (!isempty(stat_angles) | !isempty(nonstat_angles)) & !(stat_ignore & nonstat_ignore)
        push!(parameters, [Ṁ, T_max, r_mi, r_mo])
        push!(model_names, model_name)
        push!(stats_exist, stat_exist)
        push!(stats_ignore, stat_ignore)
        push!(stats_angles, stat_angles)
        push!(nonstats_exist, nonstat_exist)
        push!(nonstats_ignore, nonstat_ignore)
        push!(nonstats_angles, nonstat_angles)
        global n_models += 1
    end
end

function findprocmodels(n_proc, n_models)
    n_proc_models = n_proc
    for i = n_proc:n_models
        if n_models % i == 0
            n_proc_models = i
            break
        end
    end
    n_proc_models
end

n_proc_models = findprocmodels(n_proc, n_models)
n_iters = n_models ÷ n_proc
println(n_iters, " ", n_models)
# println(stats_angles)
# println(nonstats_angles)

@time for iter_id = 0:n_iters
    proc_iter_start = iter_id*n_proc + 1
    proc_iter_end = (iter_id+1)*n_proc
    if proc_iter_end > n_models
        proc_iter_end = n_models
    end
#     println("$proc_iter_start $proc_iter_end")
    profs = @distributed vcat for model_id = proc_iter_start:proc_iter_end
        proc_profs = []
#         model_id = n_proc_models*(iter_id-1) + proc_iter_id
        model_name = model_names[model_id]
        Ṁ, T_max, r_mi, r_mo = parameters[model_id]
        stat_exist = stats_exist[model_id]; nonstat_exist = nonstats_exist[model_id]
        stat_ignore = stats_ignore[model_id]; nonstat_ignore = nonstats_ignore[model_id]
        stat_angles = stats_angles[model_id]
        nonstat_angles = nonstats_angles[model_id]

        stat_ok = true
        mag_stat = if stat_exist
            TTauUtils.Models.loadmodel(star, "$(model_name)_stat_nonlocal")
        else
            try 
                local_mag = StationarySolidMagnetosphereNHCool("$(model_name)_stat", star, r_mi, r_mo, Ṁ, T_max, 10, n_t = 20, progress_output = false)
                addnonlocal(local_mag, progress_output = false)
            catch 
                stat_ok = false
            end
        end

        nonstat_ok = false
        # mag_nonstat = if nonstat_exist
        #     TTauUtils.Models.loadmodel(star, "$(model_name)_nonstat_nonlocal")
        # else
        #     try 
        #         local_mag = NonStationarySolidMagnetosphereNHCool("$(model_name)_nonstat", star, r_mi, r_mo, Ṁ, T_max, 10, n_t = 20, progress_output = false)
        #         addnonlocal(local_mag, progress_output = false)
        #     catch 
        #         nonstat_ok = false
        #     end
        # end

        if stat_ok
            for i_ang in stat_angles
#                 println(profile_name)
                prof = HydrogenProfileDoppler(mag_stat, u, l, i_ang, 0.05, 0.05, 0.05, 50, progress_output = false, blue_v_max = 300, red_v_max = 400)
                push!(proc_profs, prof)
#                 println(prof.orientation)
            end
        end

        if nonstat_ok
            for i_ang in nonstat_angles
#                 profile_name = "$(linename(u,l))_$i_ang"
#                 println(profile_name)
                prof = HydrogenProfileDoppler(mag_nonstat, u, l, i_ang, 0.05, 0.05, 0.05, 50, progress_output = false, blue_v_max = 300, red_v_max = 400)
                push!(proc_profs, prof)
#                 println(prof.orientation)
            end
        end


#         println("$model_id $iter_id $model_name")
        proc_profs
    end
    println(iter_id)
    for prof in profs
        i_ang = round(Int, prof.orientation.i/π*180)
#         u = prof.upper_level
#         l = prof.lower_level
        profile_name = "$(linename(u,l))_$i_ang"
        saveprofile(prof, profile_name)
        println(profile_name, " ", prof.model.name)
    end
    println("$iter_id from $n_iters")
end

opt = SendOptions(
  isSSL = true,
  username = "my.calc.results@gmail.com",
  passwd = "vuN2jtpUXgApc6t")
#Provide the message body as RFC5322 within an IO
body = IOBuffer("Models computed")
url = "smtps://smtp.gmail.com:465"
rcpt = ["<dmitrievdv242@gmail.com>"]
from = "<my.calc.results@gmail.com>"
resp = send(url, rcpt, from, body, opt)
