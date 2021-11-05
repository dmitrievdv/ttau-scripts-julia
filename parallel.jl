using Distributed

addprocs(6)

@everywhere using TTauUtils

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
T_maxs = [8000:1000:12000;]
lg10_Ṁs = [-11:0.5:-8.5;]

i_angs = [40:2:50;]
u = 3; l =2

for r_mi in r_mis
    for mag_width in mag_widths
        r_mo = r_mi + mag_width
        mag_string = "$(floor(Int, r_mi))-$(floor(Int, r_mo))"
        for T_max in T_maxs
            T_string = string(floor(Int, T_max))
            profs = @distributed vcat for lg10_Ṁ in lg10_Ṁs
                proc_profs = []
                Ṁ_string = string(floor(Int, -lg10_Ṁ*10))
                Ṁ = 10.0^lg10_Ṁ
                model_name = "$(Ṁ_string)_$(T_string)_$(mag_string)"
                stat_ok = true
                nonstat_ok = true
                mag_stat = try 
                    local_mag = StationarySolidMagnetosphereNHCool("$(model_name)_stat", star, r_mi, r_mo, Ṁ, T_max, 10, n_t = 40, progress_output = false)
                    nonlocal_mag = addnonlocal(local_mag, progress_output = false)
                    nonlocal_mag
                catch 
                    stat_ok = false
                end
                mag_nonstat = try 
                    local_mag = NonStationarySolidMagnetosphereNHCool("$(model_name)_nonstat", star, r_mi, r_mo, Ṁ, T_max, 10, n_t = 40, progress_output = false)
                    nonlocal_mag = addnonlocal(local_mag, progress_output = false)
                    nonlocal_mag
                catch 
                    nonstat_ok = false
                end
                if stat_ok
                    for i_ang in i_angs
                        prof = HydrogenProfileDoppler(mag_stat, u, l, i_ang, 0.05, 0.1, 0.1, 200, progress_output = false)
                        push!(proc_profs, prof)
                    
                    end
                end
                if nonstat_ok
                    for i_ang in i_angs
                        prof = HydrogenProfileDoppler(mag_nonstat, u, l, i_ang, 0.05, 0.1, 0.1, 200, progress_output = false)
                        push!(proc_profs, prof)
                    end
                end
                proc_profs
            end
            for prof in profs
                i_ang = floor(Int, prof.orientation.i/π*180)
                u = prof.upper_level
                l = prof.lower_level
                saveprofile(prof, "$(linename(u,l))_$i_ang")
            end
        end
    end
end