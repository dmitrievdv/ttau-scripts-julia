using Distributed

n_proc = 6
addprocs(n_proc)

@everywhere using TTauUtils

@everywhere function modelname(Ṁ, T_max, r_mi, r_mo)
    mag_string = "$(round(Int, r_mi))-$(round(Int, r_mo))"
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

star_name = "RZPsc"
star = Star(star_name)

r_mis = [6.0, 7.0, 8.0]
mag_widths = [1.0, 2.0, 3.0]
T_maxs = [8000:500:12000;]
lg10_Ṁs = [-11:0.1:-8.5;]

i_angs = [40:2:50;]
u = 3; l =2

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
    mag_string = "$(round(Int, r_mi))-$(round(Int, r_mo))"
    T_string = string(round(Int, T_max))
    Ṁ_string = string(round(Int, -lg10_Ṁ*10))
    Ṁ = 10.0^lg10_Ṁ
    model_name = modelname(Ṁ, T_max, r_mi, r_mo)
    stat_exist = isfile("stars/$star_name/$(model_name)_stat_nonlocal/$(model_name)_stat_nonlocal.dat")
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

    

    if (!isempty(stat_angles) | !isempty(nonstat_angles)) & !stat_exist
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
n_iters = n_models ÷ n_proc_models
println(n_iters, " ", n_proc_models, " ", n_models)
# println(stats_angles)
# println(nonstats_angles)

@time for iter_id = 1:n_iters
    profs = @distributed vcat for proc_iter_id = 1:n_proc_models
        proc_profs = []
        model_id = n_proc_models*(iter_id-1) + proc_iter_id
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
                local_mag = StationarySolidMagnetosphereNHCool("$(model_name)_stat", star, r_mi, r_mo, Ṁ, T_max, 10, n_t = 40, progress_output = false)
                addnonlocal(local_mag, progress_output = false)
            catch 
                stat_ok = false
            end
        end

        nonstat_ok = true
        mag_nonstat = if nonstat_exist
            TTauUtils.Models.loadmodel(star, "$(model_name)_stat_nonlocal")
        else
            try 
                local_mag = StationarySolidMagnetosphereNHCool("$(model_name)_nonstat", star, r_mi, r_mo, Ṁ, T_max, 10, n_t = 40, progress_output = false)
                addnonlocal(local_mag, progress_output = false)
            catch 
                nonstat_ok = false
            end
        end

        if stat_ok
            for i_ang in stat_angles
                # println(profile_name)
                prof = HydrogenProfileDoppler(mag_stat, u, l, i_ang, 0.05, 0.1, 0.1, 200, progress_output = false)
                push!(proc_profs, prof)
                # println(prof.orientation)
            end
        end

        if nonstat_ok
            for i_ang in nonstat_angles
                # profile_name = "$(linename(u,l))_$i_ang"
                # println(profile_name)
                prof = HydrogenProfileDoppler(mag_nonstat, u, l, i_ang, 0.05, 0.1, 0.1, 200, progress_output = false)
                push!(proc_profs, prof)
                # println(prof.orientation)
            end
        end


        # println("$model_id $iter_id $model_name")
        proc_profs
    end
    println(iter_id)
    for prof in profs
        i_ang = round(Int, prof.orientation.i/π*180)
        u = prof.upper_level
        l = prof.lower_level
        # profile_name = "$(linename(u,l))_$i_ang"
        saveprofile(prof, profile_name)
        # println(profile_name)
    end
    println("$iter_id from $n_iters")
end

# @time for r_mi in r_mis
#     for mag_width in mag_widths
#         r_mo = r_mi + mag_width
#         mag_string = "$(round(Int, r_mi))-$(round(Int, r_mo))"
#         for T_max in T_maxs
#             T_string = string(round(Int, T_max))
#             profs = @distributed vcat for lg10_Ṁ in lg10_Ṁs
#                 proc_profs = []
#                 Ṁ_string = string(round(Int, -lg10_Ṁ*10))
#                 Ṁ = 10.0^lg10_Ṁ
#                 model_name = "$(Ṁ_string)_$(T_string)_$(mag_string)"
#                 stat_ok = true
#                 nonstat_ok = true
#                 mag_stat = if !isfile("stars/$star_name/$(model_name)_stat_nonlocal/$(model_name)_stat_nonlocal.dat")
#                     try 
#                         local_mag = StationarySolidMagnetosphereNHCool("$(model_name)_stat", star, r_mi, r_mo, Ṁ, T_max, 10, n_t = 40, progress_output = false)
#                         addnonlocal(local_mag, progress_output = false)
#                     catch 
#                         stat_ok = false
#                     end
#                 else
#                     TTauUtils.Models.loadmodel(star, "$(model_name)_stat_nonlocal")
#                 end
#                 println(mag_stat.Mdot)
#                 mag_nonstat = if !isfile("stars/$star_name/$(model_name)_nonstat_nonlocal/$(model_name)_nonstat_nonlocal.dat")
#                     try 
#                         local_mag = NonStationarySolidMagnetosphereNHCool("$(model_name)_nonstat", star, r_mi, r_mo, Ṁ, T_max, 10, n_t = 40, progress_output = false)
#                         addnonlocal(local_mag, progress_output = false)
#                     catch 
#                         nonstat_ok = false
#                     end
#                 else
#                     nonstat_ok = true
#                     TTauUtils.Models.loadmodel(star, "$(model_name)_nonstat_nonlocal")
#                 end
#                 println(mag_nonstat.Mdot)
#                 if stat_ok
#                     for i_ang in i_angs
#                         profile_name = "$(linename(u,l))_$i_ang"
#                         if !isfile("stars/$star_name/$(model_name)_stat_nonlocal/$profile_name.dat")
#                             prof = HydrogenProfileDoppler(mag_stat, u, l, i_ang, 0.05, 0.1, 0.1, 200, progress_output = false)
#                             push!(proc_profs, prof)
#                             println(prof.orientation)
#                         end
#                     end
#                 end
#                 if nonstat_ok
#                     for i_ang in i_angs
#                         if !isfile("stars/$star_name/$(model_name)_nonstat_nonlocal/$profile_name.dat")
#                             prof = HydrogenProfileDoppler(mag_nonstat, u, l, i_ang, 0.05, 0.1, 0.1, 200, progress_output = false)
#                             push!(proc_profs, prof)
#                             println(prof.orientation)
#                         end
#                     end
#                 end
#                 proc_profs
#             end
#             for prof in profs
#                 i_ang = round(Int, prof.orientation.i/π*180)
#                 u = prof.upper_level
#                 l = prof.lower_level
#                 saveprofile(prof, "$(linename(u,l))_$i_ang")
#             end
#         end
#     end
# end
