using TTauUtils
using Plots
using Plots.Measures
using LaTeXStrings
using Printf
using Dierckx


function getprofile(prof :: HydrogenProfile)
    r_ν = prof.profile
    ν_0 = TTauUtils.HydrogenPopulations.linefrequency(prof.upper_level, prof.lower_level)
    vs = @. (ν_0 - prof.frequencies)/ν_0*3e5 
    return vs, r_ν
end

function eqwidth(prof :: HydrogenProfile)
    vs, r_ν = getprofile(prof)
    n = prof.frequency_points
    dv = zeros(n)
    dv[1] = vs[2] - vs[1]
    for i=2:n
        dv[i] = 2*(vs[i] - vs[i-1]) - dv[i-1]
    end
    fs = @. (r_ν - 1)*dv
    return maximum(r_ν .- 1)
end

function maxintensity(prof :: HydrogenProfile)
    vs, r_ν = getprofile(prof)
    return maximum(r_ν) - 1
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

function plotall(star_name)

    star = Star(star_name)
    lines = ["Ha", "Hb"]; n_lines = length(lines)
    i_angs = [10:10:80;]; n_angs = length(i_angs)
    # model = SolidMagnetosphere(star, "hart94_70_7500_2-3_stat_nonlocal")

    models = readdir("stars/$star_name")[2:end]
    for model_name in models
        model = SolidMagnetosphere(star, model_name)
    

        layout = @layout grid(n_lines, n_angs)
        plts = []
        profs = []
        y_maxs = fill(1.0, n_lines); y_mins = fill(1.0, n_lines)
        mult = 1.1
        for (l, line) in enumerate(lines)
            for i_ang in i_angs
                prof = HydrogenProfile(star, model, "$(line)_$i_ang")
                prof_max = maximum(prof.profile); prof_min = minimum(prof.profile)
                if prof_max*mult > y_maxs[l]
                    y_maxs[l] = prof_max*mult
                end
                if prof_min/mult < y_mins[l]
                    y_mins[l] = prof_min/mult
                end
                push!(profs, prof)
            end
        end

        for l=1:n_lines
            y_min = y_mins[l]; y_max = y_maxs[l]
            line = lines[l]
            for i=1:n_angs
                i_ang = i_angs[i]
                k = (l-1)*n_angs + i
                vs, r_ν = getprofile(profs[k])
                plt = plot(vs, r_ν, label = false, title = "$(line)_$i_ang\nmin = $(minimum(r_ν)), max = $(maximum(r_ν))", ylims = (y_min, y_max))
                push!(plts, plt)
            end
        end
        plt = plot(plts..., layout = layout, size = (400*length(i_angs), 400*length(lines)), margin = 5mm)
        savefig(plt, "stars/$(star.name)/plots/$(model.name).pdf")
    end
end


begin
    star_name = "hart94"
    star = Star(star_name)
    model_firstname = "$(star_name)_93_10000_2-3"
    model_name_1 = "$(model_firstname)_stat_nonlocal"
    model_name_2 = "$(model_firstname)_nonstat_nonlocal"
    model_1 = SolidMagnetosphere(star, model_name_1)
    model_2 = model_1
    nonstat_undef = false
    try
        model_2 = SolidMagnetosphere(star, model_name_2)
    catch
        model_2 = model_1
        nonstat_undef = true
    end

    u, l = 5, 2

    line_name = linename(u, l)
    i_angs = [80;]
    layout = @layout grid(1,1)
    plts = []

    for i_ang in i_angs
        prof_1 = HydrogenProfile(star, model_1, "$(line_name)_$i_ang")
        
        if nonstat_undef
            plt = plot(getprofile(prof_1)..., title = "i = $i_ang", label = false, lc = :black, xlabel = L"v,\ \mathrm{km/s}", ylabel = L"r_\nu")
        else
            plt = plot(getprofile(prof_1)..., title = "i = $i_ang", label = false, lc = :red)
            prof_2 = HydrogenProfile(star, model_2, "$(line_name)_$i_ang")
            plot!(plt, getprofile(prof_2)..., label = false, lc = :blue, xlabel = L"v,\ \mathrm{km/s}", ylabel = L"r_\nu")
        end
        push!(plts, plt)
    end
    prof_plot = plot(plts..., layout = layout, size = (800, 800), margin = 5mm)
end

begin 
    star_name = "coldhart94"
    star = Star(star_name)
    model_firstname = "$(star_name)_84_9000_4-6"
    model_name_1 = "$(model_firstname)_stat_nonlocal"
    model_name_2 = "$(model_firstname)_nonstat_nonlocal"
    model_1 = SolidMagnetosphere(star, model_name_1)
    nonstat_undef = false
    try
        model_2 = SolidMagnetosphere(star, model_name_2)
    catch
        model_2 = model_1
        nonstat_undef = true
    end
    
    u, l = 7, 4

    r_ms = model_1.r_m_grid
    t_points = model_1.t_grid
    t_lines = collect(range(0, 1, length = 100))
    
    T_e_plot = plot()
    n_h_plot = plot()
    f_plot = plot()
    S_plot = plot()
    n_u_plot = plot()
    n_l_plot = plot()

    for r_m in r_ms
        θs = @. asin(√(1/r_m)) + t_points * (π/2 - asin(√(1/r_m)))
        r_points = @. r_m*sin(θs)^2
        θs = @. asin(√(1/r_m)) + t_lines * (π/2 - asin(√(1/r_m)))
        r_lines = @. r_m*sin(θs)^2
        n_e_1_lines= @. 10^model_1.lgne_spl2d.(r_m, t_lines)
        n_h_1_lines= @. 10^model_1.lgnh_spl2d.(r_m, t_lines)

        T_e_1_lines= @. model_1.Te_spl2d.(r_m, t_lines)
        T_e_2_lines= @. model_2.Te_spl2d.(r_m, t_lines)
        
        n_e_2_lines= @. 10^model_2.lgne_spl2d.(r_m, t_lines)
        n_h_2_lines= @. 10^model_2.lgnh_spl2d.(r_m, t_lines)
        
        n_u_1_lines= @. 10^model_1.lgni_spl2d[u].(r_m, t_lines)
        n_l_1_lines= @. 10^model_1.lgni_spl2d[l].(r_m, t_lines)

        n_u_2_lines= @. 10^model_2.lgni_spl2d[u].(r_m, t_lines)
        n_l_2_lines= @. 10^model_2.lgni_spl2d[l].(r_m, t_lines)

        n_e_1_points= @. 10^model_1.lgne_spl2d.(r_m, t_points)
        n_h_1_points= @. 10^model_1.lgnh_spl2d.(r_m, t_points)

        T_e_1_points= @. model_1.Te_spl2d.(r_m, t_points)
        T_e_2_points= @. model_2.Te_spl2d.(r_m, t_points)
        
        n_e_2_points= @. 10^model_2.lgne_spl2d.(r_m, t_points)
        n_h_2_points= @. 10^model_2.lgnh_spl2d.(r_m, t_points)
        
        n_u_1_points= @. 10^model_1.lgni_spl2d[u].(r_m, t_points)
        n_l_1_points= @. 10^model_1.lgni_spl2d[l].(r_m, t_points)

        n_u_2_points= @. 10^model_2.lgni_spl2d[u].(r_m, t_points)
        n_l_2_points= @. 10^model_2.lgni_spl2d[l].(r_m, t_points)

        


        # plot!(plt, rs, (@. 1/(n_l_1s/n_u_1s*u^2/l^2 - 1)), lc = :red)
        # plot!(plt, rs, (@. 1/(n_l_2s/n_u_2s*u^2/l^2 - 1)), lc = :blue)
    
        plot!(T_e_plot, r_lines, T_e_1_lines, lc = :black, legend = false, ylabel = L"T_e,\ \mathrm{K}", xlabel = L"r,\ \mathrm{R}_\star")
        scatter!(T_e_plot, r_points, T_e_1_points, mc = :black, ms = 2)

        plot!(n_h_plot, r_lines, (@.  log10(n_h_1_lines) ), lc = :black, legend = false, ylabel = L"\lg\ n_H,\ \mathrm{cm}^{-3}", xlabel = L"r,\ \mathrm{R}_\star")
        scatter!(n_h_plot, r_points, (@. log10(n_h_1_points) ), mc = :black, ms = 2)

        if nonstat_undef
            plot!(f_plot, r_lines, (@. n_e_1_lines/n_h_1_lines), lc = :black, ylabel = L"n_e/n_H", xlabel = L"r,\ \mathrm{R}_\star")
            scatter!(f_plot, r_points, (@. n_e_1_points/n_h_1points), mc = :black, ms = 2, legend = false)

            plot!(S_plot, r_lines, log10.(@. 1/(n_l_1_lines/n_u_1_lines*u^2/l^2 - 1)), lc = :black, ylabel = L"S_{%$u%$l},\ \mathrm{erg}\cdot\mathrm{cm}^{-2}\cdot\mathrm{s}^{-1}")
            scatter!(S_plot, r_points, log10.(@. 1/(n_l_1_points/n_u_1_points*u^2/l^2 - 1)), mc = :black, ms = 2, legend = false, xlabel = L"r,\ \mathrm{R}_\star")

            plot!(n_l_plot, r_lines, @.( log10(n_l_1_lines) ), lc = :black, ylabel = L"\lg\ n_{%$l},\ \mathrm{cm}^{-3}")
            scatter!(n_l_plot, r_points, @.( log10(n_l_1_points) ), mc = :black, ms = 2, legend = false, xlabel = L"r,\ \mathrm{R}_\star")

            plot!(n_u_plot, r_lines, @.( log10(n_u_1_lines) ), lc = :black, ylabel = L"\lg\ n_{%$u},\ \mathrm{cm}^{-3}")
            scatter!(n_u_plot, r_points, @.( log10(n_u_1_points) ), mc = :black, ms = 2, legend = false, xlabel = L"r,\ \mathrm{R}_\star")
        else
            plot!(f_plot, r_lines, (@. n_e_1_lines/n_h_1_lines), lc = :red, ylabel = L"n_e/n_H")
            plot!(f_plot, r_lines, (@. n_e_2_lines/n_h_2_lines), lc = :blue, xlabel = L"r,\ \mathrm{R}_\star", legend = false)
            scatter!(f_plot, r_points, (@. n_e_1_points/n_h_1_points), mc = :red, ms = 2, legend = false)
            scatter!(f_plot, r_points, (@. n_e_2_points/n_h_2_points), mc = :blue, ms = 2)

            plot!(S_plot, r_lines, log10.(@. 1/(n_l_1_lines/n_u_1_lines*u^2/l^2 - 1)), lc = :red, ylabel = L"S_{%$u%$l},\ \mathrm{erg}\cdot\mathrm{cm}^{-2}\cdot\mathrm{s}^{-1}")
            plot!(S_plot, r_lines, log10.(@. 1/(n_l_2_lines/n_u_2_lines*u^2/l^2 - 1)), lc = :blue, xlabel = L"r,\ \mathrm{R}_\star", legend = false)
            scatter!(S_plot, r_points, log10.(@. 1/(n_l_1_points/n_u_1_points*u^2/l^2 - 1)), mc = :red, ms = 2, legend = false)
            scatter!(S_plot, r_points, log10.(@. 1/(n_l_2_points/n_u_2_points*u^2/l^2 - 1)), mc = :blue, ms = 2)

            plot!(n_l_plot, r_lines, @.( log10(n_l_1_lines) ), lc = :red, ylabel = L"\lg\ n_{%$l},\ \mathrm{cm}^{-3}")
            plot!(n_l_plot, r_lines, @.( log10(n_l_2_lines) ), lc = :blue, xlabel = L"r,\ \mathrm{R}_\star", legend = false)
            scatter!(n_l_plot, r_points, @.( log10(n_l_1_points) ), mc = :red, ms = 2, legend = false)
            scatter!(n_l_plot, r_points, @.( log10(n_l_2_points) ), mc = :blue, ms = 2)

            plot!(n_u_plot, r_lines, @.( log10(n_u_1_lines) ), lc = :red, ylabel = L"\lg\ n_{%$u},\ \mathrm{cm}^{-3}")
            plot!(n_u_plot, r_lines, @.( log10(n_u_2_lines) ), lc = :blue, xlabel = L"r,\ \mathrm{R}_\star", legend = false)
            scatter!(n_u_plot, r_points, @.( log10(n_u_1_points) ), mc = :red, ms = 2, legend = false)
            scatter!(n_u_plot, r_points, @.( log10(n_u_2_points) ), mc = :blue, ms = 2)
        end
    end
    plt = plot(f_plot, S_plot, n_u_plot, n_l_plot, layout = @layout([A B; C D]), size = (1000,1000), dpi = 300, margin = 5mm)

    # line_name = linename(u, l)
    # i_ang = 40
    # prof_1 = HydrogenProfile(star, model_1, "$(line_name)_$i_ang")

    # obs_v, obs_r = open("observation.dat", "r") do io
    #     obs_v = Float64[]; obs_r = Float64[]
    #     for line in readlines(io)
    #         v, r = parse.(Float64, split(line))
    #         push!(obs_v, v); push!(obs_r, r)
    #     end
    #     obs_v, obs_r 
    # end

    # if nonstat_undef
    #     prof_plt = plot(getprofile(prof_1)..., title = "i = $i_ang", label = false, lc = :black, xlabel = L"v,\ \mathrm{km/s}", ylabel = L"r_\nu")
    # else
    #     prof_plt = plot(getprofile(prof_1)..., title = L"H_\alpha,\ i = %$i_ang", label = false, lc = :red)
    #     prof_2 = HydrogenProfile(star, model_2, "$(line_name)_$i_ang")
    #     plot!(prof_plt, getprofile(prof_2)..., label = false, lc = :blue, xlabel = L"v,\ \mathrm{km/s}", ylabel = L"r_\nu")
    #     plot!(prof_plt, obs_v, obs_r, label = false, lc = :black, ls = :dash)
    # end
    # plt = plot(prof_plt, size = (1000,600), dpi = 300, margin = 5mm)plt = plot(f_plot, S_plot, layout = @layout([A B]), size = (1000,600), dpi = 300, margin = 5mm)

    # line_name = linename(u, l)
    # i_ang = 40
    # prof_1 = HydrogenProfile(star, model_1, "$(line_name)_$i_ang")

    # obs_v, obs_r = open("observation.dat", "r") do io
    #     obs_v = Float64[]; obs_r = Float64[]
    #     for line in readlines(io)
    #         v, r = parse.(Float64, split(line))
    #         push!(obs_v, v); push!(obs_r, r)
    #     end
    #     obs_v, obs_r 
    # end

    # if nonstat_undef
    #     prof_plt = plot(getprofile(prof_1)..., title = "i = $i_ang", label = false, lc = :black, xlabel = L"v,\ \mathrm{km/s}", ylabel = L"r_\nu")
    # else
    #     prof_plt = plot(getprofile(prof_1)..., title = L"H_\alpha,\ i = %$i_ang", label = false, lc = :red)
    #     prof_2 = HydrogenProfile(star, model_2, "$(line_name)_$i_ang")
    #     plot!(prof_plt, getprofile(prof_2)..., label = false, lc = :blue, xlabel = L"v,\ \mathrm{km/s}", ylabel = L"r_\nu")
    #     plot!(prof_plt, obs_v, obs_r, label = false, lc = :black, ls = :dash)
    # end
    # plt = plot(prof_plt, size = (1000,600), dpi = 300, margin = 5mm)
end

# begin
#     star_name = "hart94"
#     star = Star(star_name)
#     models_all = readdir("stars/$star_name")[2:end]
#     T_e_str = "9500"
#     mag_str = "2-3"
    
#     models = String[]
#     regex = Regex("$(star_name)_[0-9]{2,2}_$(T_e_str)_$(mag_str)")
#     for model_name in models_all
#         occur = occursin(regex, model_name)
#         if occur
#             m = match(regex, model_name)
#             if !isempty(m.match) && (isempty(models) || m.match != models[end])
#                 push!(models, m.match)
#                 println(m.match)
#             end
#         end
#     end

#     regex = Regex("$(star_name)_[0-9]{3,}_$(T_e_str)_$(mag_str)")
#     for model_name in models_all
#         occur = occursin(regex, model_name)
#         if occur
#             m = match(regex, model_name)
#             if !isempty(m.match) && (isempty(models) || m.match != models[end])
#                 push!(models, m.match)
#                 println(m.match)
#             end
#         end
#     end

#     u, l = 4, 2

    

#     anim = @animate for model_name in models
#         model_name_1 = "$(model_name)_stat_nonlocal"
#         model_name_2 = "$(model_name)_nonstat_nonlocal"
#         model_1 = SolidMagnetosphere(star, model_name_1)
#         model_2 = model_1
#         nonstat_undef = false
#         try
#             model_2 = SolidMagnetosphere(star, model_name_2)
#         catch
#             model_2 = model_1
#             nonstat_undef = true
#         end
#         u, l = 4, 2
#         line_name = linename(u, l)
#         i_angs = [10:10:80;]
#         layout = @layout grid(2,4)
#         plts = []

        

#         for i_ang in i_angs
#             prof_1 = HydrogenProfile(star, model_1, "$(line_name)_$i_ang")
        
#             if nonstat_undef
#                 plt = plot(getprofile(prof_1)..., title = "i = $i_ang", label = false, lc = :black)
#             else
#                 plt = plot(getprofile(prof_1)..., title = "i = $i_ang", label = false, lc = :red)
#                 prof_2 = HydrogenProfile(star, model_2, "$(line_name)_$i_ang")
#                 plot!(plt, getprofile(prof_2)..., label = false, lc = :blue)
#             end
#             push!(plts, plt)
#         end
#         title = plot(title = "$model_name, $line_name", grid = false, axis = ([], false), showaxis = false, bottom_margin = -50Plots.px)
#         prof_plot = plot(title, plts..., layout = @layout([A{0.02h}; [A B C D; E F G H]]), size = (1600, 900), dpi = 300)
#         # prof_plot = plot(plts..., layout = layout, size = (1600, 800))
#     end
#     gif(anim, "$(line_name)_anim.gif", fps = 1)
# end

function calceqwidths(models, u, l, i; suffix = "")
    Ṁs = Float64[]
    stat_eqwidths = Float64[]
    nonstat_eqwidths = Float64[]

    for model_name in models
        if suffix != ""
            model_name = model_name*"_"*suffix
        end
        model_name_1 = "$(model_name)_stat_nonlocal"
        model_name_2 = "$(model_name)_nonstat_nonlocal"
        model_1 = SolidMagnetosphere(star, model_name_1)
        model_2 = model_1
        nonstat_undef = false
        try
            model_2 = SolidMagnetosphere(star, model_name_2)
        catch
            model_2 = model_1
            nonstat_undef = true
        end
        line_name = linename(u, l)

        Ṁ = model_1.Mdot
        prof_1 = try
            prof_1 = HydrogenProfile(star, model_1, "$(line_name)_$i")
        catch
            continue
        end
        
        push!(Ṁs, Ṁ)
        push!(stat_eqwidths, eqwidth(prof_1))
        if nonstat_undef
            push!(nonstat_eqwidths, eqwidth(prof_1))
        else
            prof_2 = HydrogenProfile(star, model_2, "$(line_name)_$i")
            push!(nonstat_eqwidths, eqwidth(prof_2))
        end
    end
    return Ṁs, stat_eqwidths, nonstat_eqwidths
end

begin
    star_name = "hart94"
    star = Star(star_name)
    models_all = readdir("stars/$star_name")[2:end]
    int_T_maxs = [7000:1000:11000;]
    # int_T_maxs = [7000,8000,9000,10000,11000]
    T_max_strings = string.(int_T_maxs)
    plt = plot()
    mag_str = "2-3"
    line_intencities = [0.1, 0.15, 0.2, 0.3, 0.5, 1.0, 2.0, 5.0]
    line_intencities_strs = ["0.1", "0.15", "0.2", "0.3", "0.5", "1.0", "2.0", "5.0"]
    n_int = length(line_intencities)
    n_T = length(T_max_strings)
    line_intencities_rels = zeros(n_int, n_T)
    line_intencities_Ṁs = zeros(n_int, n_T)
    annotations_pos = zeros(n_T)
    annotations_pos[int_T_maxs .< 9000] .= -0.05
    maximum_lgṀ = -100.0
    for (j,T_e_str) in enumerate(T_max_strings)

        models = String[]
        regex = Regex("$(star_name)_[0-9]{2,2}_$(T_e_str)_$(mag_str)")
        for model_name in models_all
            occur = occursin(regex, model_name)
            if occur
                m = match(regex, model_name)
                if !isempty(m.match) && (isempty(models) || m.match != models[end])
                    push!(models, m.match)
                    println(m.match)
                end
            end
        end

        regex = Regex("$(star_name)_[0-9]{3,}_$(T_e_str)_$(mag_str)")
        for model_name in models_all
            occur = occursin(regex, model_name)
            if occur
                m = match(regex, model_name)
                if !isempty(m.match) && (isempty(models) || m.match != models[end])
                    push!(models, m.match)
                    println(m.match)
                end
            end
        end

        u, l = 3, 2
        i = 80

        Ṁs, stat_eqwidths, nonstat_eqwidths = calceqwidths(models, u, l, i)
        rel = abs.(stat_eqwidths ./ nonstat_eqwidths)
        println(log10.(Ṁs))
        maximum_lgṀ = max(maximum_lgṀ, maximum(log10.(Ṁs)))
        # Ṁs, v0_stat_eqwidths, v0_nonstat_eqwidths = calceqwidths(models, u, l, i, suffix = "v0")
        for k = 1:n_int
            intencity = line_intencities[k]
            for l = 1:length(stat_eqwidths)-1
                if intencity > stat_eqwidths[l+1]
                    int_step = stat_eqwidths[l+1] - stat_eqwidths[l]
                    Ṁs_step = -0.1
                    rels_step = rel[l+1] - rel[l]
                    lgṀ = log10(Ṁs[l]) + Ṁs_step/int_step*(intencity - stat_eqwidths[l])
                    line_intencities_rels[k, j] = rel[l] + rels_step/int_step*(intencity - stat_eqwidths[l])
                    line_intencities_Ṁs[k,j] = 10^lgṀ
                    # println(int_step, " ", log10(Ṁs[l]), " ", rels_step)
                    break
                end
            end
        end

        # line_intencities_Ṁs = 10 .^ line_intencities_Ṁs
        
        # rel_rel = rel[nonstat_eqwidths .> 0.05 ]
        # v0_rel = abs.(v0_stat_eqwidths ./ v0_nonstat_eqwidths)
        rel_rel = []
        strs = []
        last = 0
        println(length(rel_rel))
        for k = 1:length(rel)
        
            eqw = nonstat_eqwidths[k]
            str = @sprintf("%4.2f", eqw)
            if (stat_eqwidths[k] > 0.1)
                last += 1
                push!(rel_rel, stat_eqwidths[k] / nonstat_eqwidths[k])
                push!(strs, "")
            else
                push!(rel_rel, NaN)
                # push!(strs, "")
            end
        end
        # push!(rel_rel, NaN)
        # push!(Ṁs, NaN)
        push!(strs, T_e_str)
        rel_rel[last+1] = line_intencities_rels[1,j]
        Ṁs[last+1] = line_intencities_Ṁs[1,j]
        if int_T_maxs[j] % 1000 == 0
        plot!(plt, Ṁs[1:last+1], log10.(rel_rel[1:last+1]), lc = :black, xaxis = :log, ylims = (-0.1,1), yticks = (log10.(1:10), string.(1:10)),
                    label = false, xminorgrid = true)
        
        scatter!(plt, Ṁs[1:last+1], log10.(rel_rel[1:last+1]), mc = :black, ms = 3, xaxis = :log,
                    label = false, xminorgrid = true)
        
        annotate!(plt, Ṁs[last+1], log10(rel_rel[last+1]) + annotations_pos[j], text(T_e_str, :bottom, 7))
        plot!(xlabel = L"\lg \dot{M}", ylabel = L"\mathrm{max}\{F^\mathrm{stat}_{\mathrm{H}_\alpha}\}/\mathrm{max}\{F^\mathrm{adv}_{\mathrm{H}_\alpha}\}")
        end
    end

    for k = 1:n_int
        intencity = line_intencities[k]
        for j=1:n_T-1
            lgṀ = log10(line_intencities_Ṁs[k, j])
            if lgṀ > maximum_lgṀ
                rel_step = line_intencities_rels[k, j+1] - line_intencities_rels[k, j]
                Ṁs_step = log10(line_intencities_Ṁs[k, j+1]) - lgṀ
                line_intencity_rel = line_intencities_rels[k, j] + rel_step/Ṁs_step*(maximum_lgṀ - lgṀ)
                println(line_intencity_rel)
                line_intencities_rels[k, 1:j] .= line_intencity_rel
                line_intencities_Ṁs[k,1:j] .= 10^maximum_lgṀ
            end
        end
    end

    for k = 2:n_int
        rels = line_intencities_rels[k,:]
        Ṁs = line_intencities_Ṁs[k,:]
        str = line_intencities_strs[k]
        plot!(plt, Ṁs, log10.(rels), ls = :dot, lc = :black, label = false)
        annotate!(plt, Ṁs[end]/1.1, log10.(rels[end])-0.01, text(str, :top, "Italics", 8))
    end

    plt
    # plot!(plt, Ṁs, v0_rel, label = L"v_\mathrm{start} = 0", lc = :black, xaxis = :log, ylims = (1,5))
    # strs = []
    # for (k,eqw) in enumerate(v0_nonstat_eqwidths)
    #     str = @sprintf("%3.1f", eqw)
    #     if (abs(v0_rel[k] - rel[k]) > 0.1) & (v0_nonstat_eqwidths[k] > 0.05)
    #         push!(strs, str)
    #     else
    #         push!(strs, "")
    #     end
    # end
    # scatter!(plt, Ṁs, v0_rel, label = false, ms = 3, mc = :black, series_annotations = text.(strs, :bottom, 7))
    # vline!(plt, [10^(-9.5), 10^(-9.7)], label = false, c = :black, ls = :dot, xlabel = L"\dot{M}")
    # annotate!(plt, 10^(-9.5), 3, text("Radiative transitions become\n more important than collisional", :bottom, 8, rotation = 90))
    # annotate!(plt, 10^(-9.7), 3, text("Line becomes optically thin", :bottom, 8, rotation = 90), ylabel = L"(\mathrm{max}(r^\mathrm{stat}_\nu) - 1)/(\mathrm{max}(r^\mathrm{nonstat}_\nu) - 1)")
    # plot!(plt, Ṁs, log10.(abs.(nonstat_eqwidths)), label = "nonstat", lc = :blue)
    # scatter!(plt, Ṁs, log10.(abs.(nonstat_eqwidths)), label = false, ms = 2, mc = :blue)
    # savefig(plt, "rel_v0_hart94_9500.pdf")
end

begin
for star_name in ["coldhart94"]
    # star_name = "hart94"
    star = Star(star_name)
    models_all = readdir("stars/$star_name")[2:end]
    int_T_maxs = [8000:1000:12000;]
    # int_T_maxs = [7000,8000,9000,10000,11000]
    T_max_strings = string.(int_T_maxs)
    plt = plot()
    
    diag_x = [-10.5, -8.6]
    ylims = [-11, -8.5]
    plt2 = plot(xlims = 10 .^ diag_x, ylims = 10 .^ ylims, xlabel = L"\dot{M}_\mathrm{adv},\ \mathrm{M_\odot/yr}", ylabel = L"\dot{M}_\mathrm{stat},\ \mathrm{M_\odot/yr}")
    plot!(plt2, 10 .^ diag_x, 10 .^ diag_x, ls = :dash, label = false, legend = :topleft, lc = :gray)
    plot!(plt2, 10 .^ diag_x, 10 .^ (diag_x .- log10(2)), ls = :dot, lc = :gray, label = false, xaxis = :log)
    plot!(plt2, 10 .^ diag_x, 10 .^ ( diag_x .- log10(4)), ls = :dot, lc = :gray, label = false, yaxis = :log)
    mag_str = "4-6"
    r_mi, r_mo = parse.(Int, split(mag_str, '-'))
    plot!(plt2, title = L"T_\star = %$(Int(star.T))\ \mathrm{K};\ r_\mathrm{mi} = %$r_mi;\ r_\mathrm{mo} = %$r_mo")
    
    for (j,T_e_str) in enumerate(T_max_strings)

        models = String[]
        regex = Regex("$(star_name)_[0-9]{2,2}_$(T_e_str)_$(mag_str)")
        for model_name in models_all
            occur = occursin(regex, model_name)
            if occur
                m = match(regex, model_name)
                if !isempty(m.match) && (isempty(models) || m.match != models[end])
                    push!(models, m.match)
                    println(m.match)
                end
            end
        end

        regex = Regex("$(star_name)_[0-9]{3,}_$(T_e_str)_$(mag_str)")
        for model_name in models_all
            occur = occursin(regex, model_name)
            if occur
                m = match(regex, model_name)
                if !isempty(m.match) && (isempty(models) || m.match != models[end])
                    push!(models, m.match)
                    println(m.match)
                end
            end
        end

        u, l = 3, 2
        is = [40]
        line_styles = [:solid]
        for k in 1:length(is)
            i = is[k]
            line_style = line_styles[k]
            Ṁs, stat_eqwidths, nonstat_eqwidths = calceqwidths(models, u, l, i)
            stat_lgṀs = log10.(Ṁs[stat_eqwidths .> 0.05])
            lgstat = log10.(stat_eqwidths[stat_eqwidths .> 0.05])

            nonstat_lgṀs = log10.(Ṁs[nonstat_eqwidths .> 0.05])
            lgnonstat = log10.(nonstat_eqwidths[nonstat_eqwidths .> 0.05])

            lgstat_sorted_id = sortperm(lgstat)
            lgnonstat_sorted_id = sortperm(lgnonstat)

            lgṀ_stat_spl = Spline1D(lgstat[lgstat_sorted_id], stat_lgṀs[lgstat_sorted_id], k = 1)
            lgṀ_nonstat_spl = Spline1D(lgnonstat[lgnonstat_sorted_id], nonstat_lgṀs[lgnonstat_sorted_id], k = 1)

            lgint = [-1:0.1:1;]
            xticks = [0.2:0.2:1; 2:2:10]
            # yticks = [1:]
        
            plt_x = lgint
            plt_y = 10 .^ (lgṀ_nonstat_spl.(lgint) .- lgṀ_stat_spl.(lgint))

            max_y, max_id = findmax(plt_y)
            max_x = plt_x[max_id]
        # max_y = plt_y[max_id]

            plot!(plt, lgint, 10 .^ (lgṀ_nonstat_spl.(lgint) .- lgṀ_stat_spl.(lgint)), xticks = (log10.(xticks), string.(xticks)), label = T_e_str)
            annotate!(plt, max_x, max_y, text(T_e_str, :bottom, "Italics", 8))

            plt_x = lgṀ_nonstat_spl.(lgint)
            plt_y = lgṀ_stat_spl.(lgint)
            ids = abs.(plt_x - plt_y) .> 0.01
            max_x, max_id = findmax(abs.(plt_x - plt_y))
            max_x = plt_x[max_id]
            max_y = plt_y[max_id]-0.05

            plot!(plt2, 10 .^ plt_x[ids], 10 .^ plt_y[ids],label = T_e_str, xminorgrid = :true, yminorgrid = :true, lc = :black, ls = line_style)
            if k == 1
                annotate!(plt2, 10^max_x, 10^max_y, text(T_e_str, :left, "Italics", 8), legend = false)
            end
        end
        # lgṀ_str = L"10^{-4}"
        # annotate!(plt, max_x, max_y, text(lgṀ_nonstat_spl(max_x), :top, "Italics", 8))
    end
    plt
    plt2

    # plot!(plt, Ṁs, v0_rel, label = L"v_\mathrm{start} = 0", lc = :black, xaxis = :log, ylims = (1,5))
    # strs = []
    # for (k,eqw) in enumerate(v0_nonstat_eqwidths)
    #     str = @sprintf("%3.1f", eqw)
    #     if (abs(v0_rel[k] - rel[k]) > 0.1) & (v0_nonstat_eqwidths[k] > 0.05)
    #         push!(strs, str)
    #     else
    #         push!(strs, "")
    #     end
    # end
    # scatter!(plt, Ṁs, v0_rel, label = false, ms = 3, mc = :black, series_annotations = text.(strs, :bottom, 7))
    # vline!(plt, [10^(-9.5), 10^(-9.7)], label = false, c = :black, ls = :dot, xlabel = L"\dot{M}")
    # annotate!(plt, 10^(-9.5), 3, text("Radiative transitions become\n more important than collisional", :bottom, 8, rotation = 90))
    # annotate!(plt, 10^(-9.7), 3, text("Line becomes optically thin", :bottom, 8, rotation = 90), ylabel = L"(\mathrm{max}(r^\mathrm{stat}_\nu) - 1)/(\mathrm{max}(r^\mathrm{nonstat}_\nu) - 1)")
    # plot!(plt, Ṁs, log10.(abs.(nonstat_eqwidths)), label = "nonstat", lc = :blue)
    # scatter!(plt, Ṁs, log10.(abs.(nonstat_eqwidths)), label = false, ms = 2, mc = :blue)
    # savefig(plt, "rel_v0_hart94_9500.pdf")
    savefig(plt2, "accrchange/$(star_name)_$(mag_str).pdf")
end

end

function calceqwidths(models, u_1, l_1, u_2, l_2, i; suffix = "")
    Ṁs = Float64[]
    eqwidths_1 = Float64[]
    eqwidths_2 = Float64[]

    for model_name in models
        if suffix != ""
            model_name = model_name*"_"*suffix
        end
        model_name = "$(model_name)_nonstat_nonlocal"
        try
            model = SolidMagnetosphere(star, model_name)
        
            line_name_1 = linename(u_1, l_1)
            line_name_2 = linename(u_2, l_2)

            Ṁ = model.Mdot
            prof_1 = try
                prof_1 = HydrogenProfile(star, model, "$(line_name_1)_$i")
            catch
                continue
            end

            prof_2 = try
                prof_2 = HydrogenProfile(star, model, "$(line_name_2)_$i")
            catch
                continue
            end
        
            push!(Ṁs, Ṁ)
            push!(eqwidths_1, eqwidth(prof_1))
            push!(eqwidths_2, eqwidth(prof_2))
        catch
            continue
        end
    end
    return Ṁs, eqwidths_1, eqwidths_2
end

begin
    star_name = "coldhart94"
    star = Star(star_name)
    models_all = readdir("stars/$star_name")[2:end]
    int_T_maxs = [7000:1000:11000;]
    T_max_strings = string.(int_T_maxs)
    plt = plot()
    mag_str = "2-3"
    line_intencities = [0.1, 0.15, 0.2, 0.3, 0.5, 1.0, 2.0, 5.0]
    line_intencities_strs = ["0.1", "0.15", "0.2", "0.3", "0.5", "1.0", "2.0", "5.0"]
    n_int = length(line_intencities)
    n_T = length(T_max_strings)
    line_intencities_rels = zeros(n_int, n_T)
    line_intencities_Ṁs = zeros(n_int, n_T)
    annotations_pos = zeros(n_T)
    maximum_lgṀ = -100.0
    for (j,T_e_str) in enumerate(T_max_strings)

        models = String[]
        regex = Regex("$(star_name)_[0-9]{2,2}_$(T_e_str)_$(mag_str)")
        for model_name in models_all
            occur = occursin(regex, model_name)
            if occur
                m = match(regex, model_name)
                if !isempty(m.match) && (isempty(models) || m.match != models[end])
                    push!(models, m.match)
                    println(m.match)
                end
            end
        end

        regex = Regex("$(star_name)_[0-9]{3,}_$(T_e_str)_$(mag_str)")
        for model_name in models_all
            occur = occursin(regex, model_name)
            if occur
                m = match(regex, model_name)
                if !isempty(m.match) && (isempty(models) || m.match != models[end])
                    push!(models, m.match)
                    println(m.match)
                end
            end
        end

        u_1, l_1 = 5, 3
        u_2, l_2 = 7, 4
        ν_1 = TTauUtils.HydrogenPopulations.linefrequency(u_1, l_1)
        ν_2 = TTauUtils.HydrogenPopulations.linefrequency(u_2, l_2)
        cont_1 = TTauUtils.Stars.starcontinuum(star, ν_1)
        cont_2 = TTauUtils.Stars.starcontinuum(star, ν_2)
        i = 40


        Ṁs, eqwidths_1, eqwidths_2 = calceqwidths(models, u_1, l_1, u_2, l_2, i)
        rel = abs.(eqwidths_1 ./ eqwidths_2) .* cont_1/cont_2
        println(log10.(Ṁs))
        maximum_lgṀ = max(maximum_lgṀ, maximum(log10.(Ṁs)))
        # Ṁs, v0_eqwidths_1, v0_eqwidths_2 = calceqwidths(models, u, l, i, suffix = "v0")
        for k = 1:n_int
            intencity = line_intencities[k]
            for l = 1:length(eqwidths_1)-1
                if intencity > eqwidths_1[l+1]
                    int_step = eqwidths_1[l+1] - eqwidths_1[l]
                    Ṁs_step = -0.1
                    rels_step = rel[l+1] - rel[l]
                    lgṀ = log10(Ṁs[l]) + Ṁs_step/int_step*(intencity - eqwidths_1[l])
                    line_intencities_rels[k, j] = rel[l] + rels_step/int_step*(intencity - eqwidths_1[l])
                    line_intencities_Ṁs[k,j] = 10^lgṀ
                    # println(int_step, " ", log10(Ṁs[l]), " ", rels_step)
                    break
                end
            end
        end

        # line_intencities_Ṁs = 10 .^ line_intencities_Ṁs
        
        # rel_rel = rel[eqwidths_2 .> 0.05 ]
        # v0_rel = abs.(v0_eqwidths_1 ./ v0_eqwidths_2)
        rel_rel = []
        strs = []
        last = 0
        println(length(rel_rel))
        for k = 1:length(rel)
        
            eqw = eqwidths_2[k]
            str = @sprintf("%4.2f", eqw)
            if (eqwidths_1[k] > 0.1)
                last += 1
                push!(rel_rel, eqwidths_1[k] / eqwidths_2[k] * cont_1/cont_2)
                push!(strs, "")
            else
                push!(rel_rel, NaN)
                # push!(strs, "")
            end
        end
        # push!(rel_rel, NaN)
        # push!(Ṁs, NaN)
        push!(strs, T_e_str)
        rel_rel[last+1] = line_intencities_rels[1,j]
        Ṁs[last+1] = line_intencities_Ṁs[1,j]
        plot!(plt, Ṁs[1:last+1], log10.(rel_rel[1:last+1]), lc = :black, xaxis = :log, ylims = (-0.1,1), yticks = (log10.(1:10), string.(1:10)),
                    label = false, xminorgrid = true)
        
        scatter!(plt, Ṁs[1:last+1], log10.(rel_rel[1:last+1]), mc = :black, ms = 3, xaxis = :log,
                    label = false, xminorgrid = true)

        annotate!(plt, Ṁs[last+1], log10(rel_rel[last+1]) + annotations_pos[j], text(T_e_str, :bottom, 7))
        plot!(xlabel = L"\lg \dot{M}", ylabel = L"\mathrm{max}\{F^\mathrm{stat}_{\mathrm{H}_\beta}\}/\mathrm{max}\{F^\mathrm{stat}_{\mathrm{Br}_\gamma}\}")
    end

    for k = 1:n_int
        intencity = line_intencities[k]
        for j=1:n_T-1
            lgṀ = log10(line_intencities_Ṁs[k, j])
            if lgṀ > maximum_lgṀ
                rel_step = line_intencities_rels[k, j+1] - line_intencities_rels[k, j]
                Ṁs_step = log10(line_intencities_Ṁs[k, j+1]) - lgṀ
                line_intencity_rel = line_intencities_rels[k, j] + rel_step/Ṁs_step*(maximum_lgṀ - lgṀ)
                println(line_intencity_rel)
                line_intencities_rels[k, 1:j] .= line_intencity_rel
                line_intencities_Ṁs[k,1:j] .= 10^maximum_lgṀ
            end
        end
    end

    for k = 2:n_int
        rels = line_intencities_rels[k,:]
        Ṁs = line_intencities_Ṁs[k,:]
        str = line_intencities_strs[k]
        if isempty(rels[rels .<= 0.0])
            plot!(plt, Ṁs, log10.(rels), ls = :dot, lc = :black, label = false)
            annotate!(plt, Ṁs[1]*1.05, log10.(rels[1])+0.02, text(str, :left, "Italics", 8))
        end
    end

    plt
    # plot!(plt, Ṁs, v0_rel, label = L"v_\mathrm{start} = 0", lc = :black, xaxis = :log, ylims = (1,5))
    # strs = []
    # for (k,eqw) in enumerate(v0_eqwidths_2)
    #     str = @sprintf("%3.1f", eqw)
    #     if (abs(v0_rel[k] - rel[k]) > 0.1) & (v0_eqwidths_2[k] > 0.05)
    #         push!(strs, str)
    #     else
    #         push!(strs, "")
    #     end
    # end
    # scatter!(plt, Ṁs, v0_rel, label = false, ms = 3, mc = :black, series_annotations = text.(strs, :bottom, 7))
    # vline!(plt, [10^(-9.5), 10^(-9.7)], label = false, c = :black, ls = :dot, xlabel = L"\dot{M}")
    # annotate!(plt, 10^(-9.5), 3, text("Radiative transitions become\n more important than collisional", :bottom, 8, rotation = 90))
    # annotate!(plt, 10^(-9.7), 3, text("Line becomes optically thin", :bottom, 8, rotation = 90), ylabel = L"(\mathrm{max}(r^\mathrm{stat}_\nu) - 1)/(\mathrm{max}(r^\mathrm{nonstat}_\nu) - 1)")
    # plot!(plt, Ṁs, log10.(abs.(eqwidths_2)), label = "nonstat", lc = :blue)
    # scatter!(plt, Ṁs, log10.(abs.(eqwidths_2)), label = false, ms = 2, mc = :blue)
    # savefig(plt, "rel_v0_hart94_9500.pdf")
end

function stateplot(stat_model, nonstat_model, u, l)
    model_1 = stat_model
    model_2 = nonstat_model
    
    r_ms = model_1.r_m_grid
    t_points = model_1.t_grid
    t_lines = collect(range(0, 1, length = 100))
    
    T_e_plot = plot()
    n_h_plot = plot()
    f_plot = plot()
    S_plot = plot()
    n_u_plot = plot()
    n_l_plot = plot()

    for r_m in r_ms
        θs = @. asin(√(1/r_m)) + t_points * (π/2 - asin(√(1/r_m)))
        r_points = @. r_m*sin(θs)^2
        θs = @. asin(√(1/r_m)) + t_lines * (π/2 - asin(√(1/r_m)))
        r_lines = @. r_m*sin(θs)^2
        n_e_1_lines= @. 10^model_1.lgne_spl2d.(r_m, t_lines)
        n_h_1_lines= @. 10^model_1.lgnh_spl2d.(r_m, t_lines)

        T_e_1_lines= @. model_1.Te_spl2d.(r_m, t_lines)
        T_e_2_lines= @. model_2.Te_spl2d.(r_m, t_lines)
        
        n_e_2_lines= @. 10^model_2.lgne_spl2d.(r_m, t_lines)
        n_h_2_lines= @. 10^model_2.lgnh_spl2d.(r_m, t_lines)
        
        n_u_1_lines= @. 10^model_1.lgni_spl2d[u].(r_m, t_lines)
        n_l_1_lines= @. 10^model_1.lgni_spl2d[l].(r_m, t_lines)

        n_u_2_lines= @. 10^model_2.lgni_spl2d[u].(r_m, t_lines)
        n_l_2_lines= @. 10^model_2.lgni_spl2d[l].(r_m, t_lines)

        n_e_1_points= @. 10^model_1.lgne_spl2d.(r_m, t_points)
        n_h_1_points= @. 10^model_1.lgnh_spl2d.(r_m, t_points)

        T_e_1_points= @. model_1.Te_spl2d.(r_m, t_points)
        T_e_2_points= @. model_2.Te_spl2d.(r_m, t_points)
        
        n_e_2_points= @. 10^model_2.lgne_spl2d.(r_m, t_points)
        n_h_2_points= @. 10^model_2.lgnh_spl2d.(r_m, t_points)
        
        n_u_1_points= @. 10^model_1.lgni_spl2d[u].(r_m, t_points)
        n_l_1_points= @. 10^model_1.lgni_spl2d[l].(r_m, t_points)

        n_u_2_points= @. 10^model_2.lgni_spl2d[u].(r_m, t_points)
        n_l_2_points= @. 10^model_2.lgni_spl2d[l].(r_m, t_points)

        


        # plot!(plt, rs, (@. 1/(n_l_1s/n_u_1s*u^2/l^2 - 1)), lc = :red)
        # plot!(plt, rs, (@. 1/(n_l_2s/n_u_2s*u^2/l^2 - 1)), lc = :blue)
    
        plot!(T_e_plot, r_lines, T_e_1_lines, lc = :black, legend = false, ylabel = L"T_e,\ \mathrm{K}", xlabel = L"r,\ \mathrm{R}_\star")
        scatter!(T_e_plot, r_points, T_e_1_points, mc = :black, ms = 2)

        plot!(n_h_plot, r_lines, (@.  log10(n_h_1_lines) ), lc = :black, legend = false, ylabel = L"\lg\ n_H,\ \mathrm{cm}^{-3}", xlabel = L"r,\ \mathrm{R}_\star")
        scatter!(n_h_plot, r_points, (@. log10(n_h_1_points) ), mc = :black, ms = 2)

        if nonstat_undef
            plot!(f_plot, r_lines, (@. n_e_1_lines/n_h_1_lines), lc = :black, ylabel = L"n_e/n_H", xlabel = L"r,\ \mathrm{R}_\star")
            scatter!(f_plot, r_points, (@. n_e_1_points/n_h_1points), mc = :black, ms = 2, legend = false)

            plot!(S_plot, r_lines, (@. 1/(n_l_1_lines/n_u_1_lines*u^2/l^2 - 1)), lc = :black, ylabel = L"S_{%$u%$l},\ \mathrm{erg}\cdot\mathrm{cm}^{-2}\cdot\mathrm{s}^{-1}")
            scatter!(S_plot, r_points, (@. 1/(n_l_1_points/n_u_1_points*u^2/l^2 - 1)), mc = :black, ms = 2, legend = false, xlabel = L"r,\ \mathrm{R}_\star")

            plot!(n_l_plot, r_lines, @.( log10(n_l_1_lines) ), lc = :black, ylabel = L"\lg\ n_{%$l},\ \mathrm{cm}^{-3}")
            scatter!(n_l_plot, r_points, @.( log10(n_l_1_points) ), mc = :black, ms = 2, legend = false, xlabel = L"r,\ \mathrm{R}_\star")

            plot!(n_u_plot, r_lines, @.( log10(n_u_1_lines) ), lc = :black, ylabel = L"\lg\ n_{%$u},\ \mathrm{cm}^{-3}")
            scatter!(n_u_plot, r_points, @.( log10(n_u_1_points) ), mc = :black, ms = 2, legend = false, xlabel = L"r,\ \mathrm{R}_\star")
        else
            plot!(f_plot, r_lines, (@. n_e_1_lines/n_h_1_lines), lc = :red, ylabel = L"n_e/n_H")
            plot!(f_plot, r_lines, (@. n_e_2_lines/n_h_2_lines), lc = :blue, xlabel = L"r,\ \mathrm{R}_\star", legend = false)
            scatter!(f_plot, r_points, (@. n_e_1_points/n_h_1_points), mc = :red, ms = 2, legend = false)
            scatter!(f_plot, r_points, (@. n_e_2_points/n_h_2_points), mc = :blue, ms = 2)

            plot!(S_plot, r_lines, (@. 1/(n_l_1_lines/n_u_1_lines*u^2/l^2 - 1)), lc = :red, ylabel = L"S_{%$u%$l},\ \mathrm{erg}\cdot\mathrm{cm}^{-2}\cdot\mathrm{s}^{-1}")
            plot!(S_plot, r_lines, (@. 1/(n_l_2_lines/n_u_2_lines*u^2/l^2 - 1)), lc = :blue, xlabel = L"r,\ \mathrm{R}_\star", legend = false)
            scatter!(S_plot, r_points, (@. 1/(n_l_1_points/n_u_1_points*u^2/l^2 - 1)), mc = :red, ms = 2, legend = false)
            scatter!(S_plot, r_points, (@. 1/(n_l_2_points/n_u_2_points*u^2/l^2 - 1)), mc = :blue, ms = 2)

            plot!(n_l_plot, r_lines, @.( log10(n_l_1_lines) ), lc = :red, ylabel = L"\lg\ n_{%$l},\ \mathrm{cm}^{-3}")
            plot!(n_l_plot, r_lines, @.( log10(n_l_2_lines) ), lc = :blue, xlabel = L"r,\ \mathrm{R}_\star", legend = false)
            scatter!(n_l_plot, r_points, @.( log10(n_l_1_points) ), mc = :red, ms = 2, legend = false)
            scatter!(n_l_plot, r_points, @.( log10(n_l_2_points) ), mc = :blue, ms = 2)

            plot!(n_u_plot, r_lines, @.( log10(n_u_1_lines) ), lc = :red, ylabel = L"\lg\ n_{%$u},\ \mathrm{cm}^{-3}")
            plot!(n_u_plot, r_lines, @.( log10(n_u_2_lines) ), lc = :blue, xlabel = L"r,\ \mathrm{R}_\star", legend = false)
            scatter!(n_u_plot, r_points, @.( log10(n_u_1_points) ), mc = :red, ms = 2, legend = false)
            scatter!(n_u_plot, r_points, @.( log10(n_u_2_points) ), mc = :blue, ms = 2)
        end
    end
    plt = plot(f_plot, S_plot, n_u_plot, n_l_plot, layout = @layout([A B; C D]), size = (1000,1000), dpi = 300, margin = 5mm)

end


begin
    star_name = "coldhart94"
    star = Star(star_name)
    models_all = readdir("stars/$star_name")[2:end]
    int_T_max = 9000
    # int_T_maxs = [7000,8000,9000,10000,11000]
    T_e_str = string(int_T_max)
    plt = plot()
    mag_str = "4-6"
    n_T = length(T_max_strings)

    models = String[]
    regex = Regex("$(star_name)_[0-9]{2,2}_$(T_e_str)_$(mag_str)")
    for model_name in models_all
        occur = occursin(regex, model_name)
        if occur
            m = match(regex, model_name)
            if !isempty(m.match) && (isempty(models) || m.match != models[end])
                push!(models, m.match)
                println(m.match)
            end
        end
    end

    regex = Regex("$(star_name)_[0-9]{3,}_$(T_e_str)_$(mag_str)")
    for model_name in models_all
        occur = occursin(regex, model_name)
        if occur
            m = match(regex, model_name)
            if !isempty(m.match) && (isempty(models) || m.match != models[end])
                 push!(models, m.match)
                println(m.match)
            end
        end
    end

    u, l = 5, 2
    u2, l2 = 4, 2
    i = 40

    Ṁs, stat_eqwidths, nonstat_eqwidths = calceqwidths(models, u, l, i)
    rel = abs.(stat_eqwidths ./ nonstat_eqwidths)
    println(log10.(Ṁs))
    maximum_lgṀ = max(maximum_lgṀ, maximum(log10.(Ṁs)))
    # Ṁs, v0_stat_eqwidths, v0_nonstat_eqwidths = calceqwidths(models, u, l, i, suffix = "v0")

    # line_intencities_Ṁs = 10 .^ line_intencities_Ṁs
        
    # rel_rel = rel[nonstat_eqwidths .> 0.05 ]
    # v0_rel = abs.(v0_stat_eqwidths ./ v0_nonstat_eqwidths)
    rel_rel = []
    strs = []
    last_val = 0
    println(length(rel_rel))
    for k = 1:length(rel)
        
        eqw = nonstat_eqwidths[k]
        str = @sprintf("%4.2f", eqw)
        if (stat_eqwidths[k] > 0.1)
            last_val += 1
            push!(rel_rel, stat_eqwidths[k] / nonstat_eqwidths[k])
            push!(strs, "")
        else
            push!(rel_rel, NaN)
                # push!(strs, "")
        end
    end
    # push!(rel_rel, NaN)
    # push!(Ṁs, NaN)
    push!(strs, T_e_str)

    plot!(plt, Ṁs[1:last_val], log10.(rel_rel[1:last_val]), lc = :black, xaxis = :log, ylims = (-0.1,1), yticks = (log10.(1:10), string.(1:10)),
                label = false, xminorgrid = true)
        
    scatter!(plt, Ṁs[1:last_val], log10.(rel_rel[1:last_val]), mc = :black, ms = 3, xaxis = :log,
                label = false, xminorgrid = true)
        
    annotate!(plt, Ṁs[last_val], log10(rel_rel[last_val]), text(T_e_str, :bottom, 7))
    plot!(xlabel = L"\lg \dot{M}", ylabel = L"\mathrm{max}\{F^\mathrm{stat}_{\mathrm{H}_\alpha}\}/\mathrm{max}\{F^\mathrm{adv}_{\mathrm{H}_\alpha}\}")

    line_rel_plt = plt

    plots = []

    for i=1:last_val
        model = models[i]
        model_name_1 = "$(model)_stat_nonlocal"
        model_name_2 = "$(model)_nonstat_nonlocal"

        stat_model = SolidMagnetosphere(star, model_name_1)
        nonstat_model = try
            SolidMagnetosphere(star, model_name_2)
        catch
            continue
        end 
        state_plot = stateplot(stat_model, nonstat_model, u, l)
        plt = plot()

        plot!(plt, Ṁs[1:last_val], log10.(rel_rel[1:last_val]), lc = :black, xaxis = :log, ylims = (-0.1,1), yticks = (log10.(1:10), string.(1:10)),
                label = false, xminorgrid = true)
        
        scatter!(plt, Ṁs[1:last_val], log10.(rel_rel[1:last_val]), mc = :black, ms = 3, xaxis = :log,
                label = false, xminorgrid = true)

        scatter!(plt, [Ṁs[i]], [log10(rel_rel[i])], mc = :red, ms = 5)
        
        annotate!(plt, Ṁs[last_val], log10(rel_rel[last_val]), text(T_e_str, :bottom, 7))
        plot!(xlabel = L"\lg \dot{M}", ylabel = L"\mathrm{max}\{F^\mathrm{stat}_{\mathrm{H}_\alpha}\}/\mathrm{max}\{F^\mathrm{adv}_{\mathrm{H}_\alpha}\}")

        line_rel_plt = plt

        plt = plot(line_rel_plt, state_plot, layout = @layout([A B]), size = (1600, 800))
        push!(plots, plt)
    end
    
    plots

    # anim = @animate for plt in [plots; [plots[end]]]
    #     plot(plt)
    # end

    # gif(anim, "state.gif", fps = 1)
end
