using TTauUtils
using Plots
using Plots.Measures
using LaTeXStrings



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
    return sum(fs[fs .> 0])
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
    model_firstname = "$(star_name)_95_10500_2-3"
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
    i_angs = [10:20:80;]
    layout = @layout grid(2,2)
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
    star_name = "hart94"
    star = Star(star_name)
    model_firstname = "$(star_name)_90_8500_2-3"
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
    plt = plot(f_plot, S_plot, layout = @layout([A B]), size = (1000,600), dpi = 300, margin = 5mm)

    # line_name = linename(u, l)
    # i_ang = 45
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
    #     plot!(prof_plt, obs_v, obs_r, label = false, lc = :black, ls = :dash, xlims = (-300, 500))
    # end
    # plt = plot(prof_plt, size = (1000,600), dpi = 300, margin = 5mm)
end

# begin
#     star_name = "hart94"
#     star = Star(star_name)
#     models_all = readdir("stars/$star_name")[2:end]
#     T_e_str = "8500"
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

# begin
#     star_name = "hart94"
#     star = Star(star_name)
#     models_all = readdir("stars/$star_name")[2:end]
#     T_e_str = "8000"
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

#     u, l = 3, 2
#     i = 40

#     Ṁs = []
#     stat_eqwidths = []
#     nonstat_eqwidths = []

#     for model_name in models
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
#         line_name = linename(u, l)

#         Ṁ = model_1.Mdot
#         push!(Ṁs, Ṁ)
#         prof_1 = HydrogenProfile(star, model_1, "$(line_name)_$i")
#         push!(stat_eqwidths, eqwidth(prof_1))

#         if nonstat_undef
#             push!(nonstat_eqwidths, eqwidth(prof_1))
#         else
#             prof_2 = HydrogenProfile(star, model_2, "$(line_name)_$i")
#             push!(nonstat_eqwidths, eqwidth(prof_2))
#         end
#     end
#     plt = plot(Ṁs, abs.(stat_eqwidths), label = "stat", lc = :red)
#     scatter!(plt, Ṁs, abs.(stat_eqwidths), label = false, ms = 2, mc = :red)
#     plot!(plt, Ṁs, abs.(nonstat_eqwidths), label = "nonstat", lc = :blue)
#     scatter!(plt, Ṁs, abs.(nonstat_eqwidths), label = false, ms = 2, mc = :blue, xaxis = :log, yaxis = :log)
# end

begin
    star_name = "hart94"
    star = Star(star_name)
    models_all = readdir("stars/$star_name")[2:end]
    T_e_str = "10500"
    mag_str = "2-3"
    
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

    Ṁs = []
    stat_eqwidths = []
    nonstat_eqwidths = []

    for model_name in models
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
        line_name2 = linename(u2, l2)

        Ṁ = model_1.Mdot
        push!(Ṁs, Ṁ)
        prof_1 = HydrogenProfile(star, model_1, "$(line_name)_$i")
        prof_2 = HydrogenProfile(star, model_1, "$(line_name2)_$i")
        push!(stat_eqwidths, maxintensity(prof_1)/maxintensity(prof_2))

        if nonstat_undef
            push!(nonstat_eqwidths, maxintensity(prof_1)/maxintensity(prof_2))
        else
            prof_1 = HydrogenProfile(star, model_2, "$(line_name)_$i")
            prof_2 = HydrogenProfile(star, model_2, "$(line_name2)_$i")
            push!(nonstat_eqwidths, maxintensity(prof_1)/maxintensity(prof_2))
        end
    end
    plt = plot(Ṁs, abs.(stat_eqwidths), label = "stat", lc = :red)
    scatter!(plt, Ṁs, abs.(stat_eqwidths), label = false, ms = 2, mc = :red)
    plot!(plt, Ṁs, abs.(nonstat_eqwidths), label = "nonstat", lc = :blue)
    scatter!(plt, Ṁs, abs.(nonstat_eqwidths), label = false, ms = 2, mc = :blue, xaxis = :log)
end