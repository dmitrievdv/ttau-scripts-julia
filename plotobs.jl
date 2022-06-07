using Plots 


function plotmodel(pars, names, id)
    # plt = plot(xtics = (-300:100:600, string.(-300:100:600)))
    model_name, profile_name = names[id]
    model = if checkmodel(model_name, star)
        TTauUtils.Models.loadmodel(star, model_name)
    else
        calcmodel(model_name, star, pars[id, 1:4]...)
    end
    profile_nophot = if checkprofile(profile_name, model_name, star)
        HydrogenProfile(star, model, profile_name[1:5])
    else
        calcprofile(star, model_name, profile_name, pars[id, 5], 3=>2)
        HydrogenProfile(star, model, profile_name[1:5])
    end
    profile = HydrogenProfile(star, model, profile_name)
    v_mod, r_mod = getvandr(profile)
    r_obs_mod = velocitiesobstomod(v_mod, r_mod, v_obs, r_obs)
    v_mag, r_mag = getvandr(profile_nophot)
    # r_gauss_mod = r_mod .+ hotspot_gauss_model(v_mod, pars[id,6:end-1])
    # plot!(plt, v_mod, r_gauss_mod, lc = :red, label = "all")
    # plot!(plt, v_mod, hotspot_gauss_model(v_mod, pars[id,6:end-1]) .+ 1, lc = :red, la = 0.5, ls =:dash, label = "hotspot")
    plot!(plt, v_mod, r_obs_mod, xlims = (-300, 600), ylims = (0.4, 1.4), lc = :black, label = "obs")
    plot!(plt, v_mod, r_mod, lc = :red, label = "mag + phot")
    plot!(plt, v_mag, r_mag, lc = :orange, la = 1, label = "mag")
    r_obs_mod = velocitiesobstomod(v_mod, r_mod, v_obs, r_obs)
    plot!(plt, v_mod, r_obs_mod .- r_mod .+ 1, lc = :black, ls = :dash, la = 0.3, label = "mod - obs")
end

function plotmodelerr(pars, names, id)
    plt = plot()
    model_name, profile_name = names[id]
    println(names[id], " ", pars[id,5:end])
    model = TTauUtils.Models.loadmodel(star, model_name)
    profile_nophot = HydrogenProfile(star, model, profile_name[1:5])
    profile = HydrogenProfile(star, model, profile_name)
    v_mod, r_mod = getvandr(profile)
    v_mag, r_mag = getvandr(profile_nophot)
    r_gauss_mod = r_mod .+ hotspot_gauss_model(v_mod, pars[id,6:end-1])
    r_obs_mod = velocitiesobstomod(v_mod, r_gauss_mod, v_obs, r_obs)
    plot!(plt, v_mod, r_gauss_mod .- r_obs_mod, lc = :red, label = "mod - obs")
end

function plotδ(pars, names, δ_cut; la = 0.5)
    plt = plot()
    # min_δ = min(min_stat_δ, min_nonstat_δ)
    best_names = []
    best_pars = []
    for i=1:length(names)
        model_name, profile_name = names[i]
        type = split(model_name, '_')[4]
        # if type == "stat"
        #     min_δ = min_stat_δ
        # elseif type == "nonstat"
        #     min_δ = min_nonstat_δ
        # end
        if pars[i,2] < 10000; continue; end
        if pars[i,end] < δ_cut
            println(names[i], " ", pars[i,5:end])
            model = TTauUtils.Models.loadmodel(star, model_name)
            profile_nophot = HydrogenProfile(star, model, profile_name[1:5])
            profile = HydrogenProfile(star, model, profile_name)
            v_mod, r_mod = getvandr(profile)
            v_mag, r_mag = getvandr(profile_nophot)
            r_gauss_mod = r_mod .+ hotspot_gauss_model(v_mod, abs.(pars[i,6:end-1]))
            r_obs_mod = velocitiesobstomod(v_mod, r_mod, v_obs, r_obs)
            plot!(plt, v_mod, r_obs_mod, xlims = (-300, 600), ylims = (0.4, 1.4), lc = :skyblue, lw = 2, label = "obs")
            plot!(plt, v_mod, r_gauss_mod, lc = :black, ls = :dash, la = la, lw = 2, label = "mag + synt")
            # plot!(plt, v_mod, hotspot_gauss_model(v_mod, abs.(pars[i,6:end-1])) .+ 1, lc = :orange, la = la, lw = 2, label = "chromo")
            # plot!(plt, v_mod, r_mod, lc = :orange, ls = :dash, la = la, lw = 2)
            plot!(plt, v_mag, r_mag, lc = :red, la = la, lw = 2, label = "mag")
            
            # push!(plots, plt)
            push!(best_names, [model_name, profile_name])
            push!(best_pars, [i; pars[i,:]])
        end
    end
    # plot!(plt, v_obs, r_obs, xlims = (-500, 600), ylims = (0.4, 1.2), lc = :black)

    for (name, pars) in zip(best_names, best_pars)
        println(name)
        println(pars)
    end
    plt
end

function plotheat(pars, x_grid, y_grid)
    gridded_pars, gridded_names = putongrid(pars, fill(["", ""], length(pars[:,1])), x_grid, y_grid)
    xs = if x_grid isa Tuple
        x_grid[1]
    else
        x_grid
    end
    ys = if y_grid isa Tuple
        y_grid[1]
    else
        y_grid
    end
    heatmap(xs, ys,  gridded_pars[9, :, :]', clims = (0, 0.1))
end

function plotδheatTM(gridded_pars, lgṀs, T_maxs, i_rm, i_W, i_ang; clims = (0,0.1))
    ang = angs[i_ang]
    r_mi = r_mis[i_rm]
    W = Ws[i_W]
    r_mo = r_mi + W
    title = L"i = %$(ang)\degree,\ r_\mathrm{mi} = %$(r_mi),\ r_\mathrm{mo} = %$(r_mo)"
    heatmap(T_maxs, lgṀs, log10.(gridded_pars[9,:,:,i_rm,i_W,i_ang]), title = title, clims = clims)
end

function plotδheat(gridded_pars, dims, xs, ys; clims = (0, 0.1))
    # xs = selectdim(gridded_pars, 1, dims[1])
    δ_i = size(gridded_pars)[1]
    griddedδ = selectdim(gridded_pars, 1, δ_i)
    hm_size = size(griddedδ, dims[1]), size(griddedδ, dims[2])
    n_1, n_2 = hm_size
    mindims = deleteat!([1:ndims(griddedδ);], dims)
    println(mindims)
    println(hm_size)
    hm = reshape(minimum(griddedδ, dims = mindims), hm_size)
    heatmap(ys, xs, hm, clims = clims)
end

function plotδfromgrid(gridded_pars, gridded_names, i_rm, i_W, i_ang, δ_cut)
    local names = vec(gridded_names[:,:,i_rm, i_W, i_ang])
    local pars = transpose(reshape(vec(gridded_pars[:,:,:,i_rm, i_W, i_ang]), (9, length(vec(gridded_pars[1,:,:,i_rm, i_W, i_ang])))))
    plotδ(pars, names, 0.02)
end

function addHb(star, model_name, profile_name; sinicorr = true)
    model = loadmodel(star, model_name)
    prof_Ha = HydrogenProfile(star, model, profile_name)
    prof_Ha_nos = HydrogenProfile(star, model, profile_name[1:5])
    prof_Ha = TTauUtils.addphotosphespecdoppler(prof_Ha_nos, 0.05, "spec/M_p5250g4.0z-5.00t1.0_a+0.40c0.00n0.00o+0.40r0.00s0.00_VIS.spec", sinicorr = sinicorr, Δλ = 1.98)
    sinicorr = occursin("sini", profile_name)
    i = parse(Float64, split(profile_name, "_")[2])
    prof_Hb_nos = HydrogenProfileDoppler(model, 4, 2, 1.0*i, 0.05, 0.05, 0.1, 200)
    prof_Hb = TTauUtils.addphotosphespecdoppler(prof_Hb_nos, 0.05, "spec/M_p5250g4.0z-5.00t1.0_a+0.40c0.00n0.00o+0.40r0.00s0.00_VIS.spec", sinicorr = sinicorr, Δλ = 1.4)
    Hb_name = "Hb"*profile_name[3:end]
    saveprofile(prof_Hb, Hb_name)
    v_Ha_mod, r_Ha_mod = getvandr(prof_Ha)
    v_Hb_mod, r_Hb_mod = getvandr(prof_Hb)
    v_Ha_nos, r_Ha_nos = getvandr(prof_Ha_nos)
    v_Hb_nos, r_Hb_nos = getvandr(prof_Hb_nos)
    v_Ha_obs, r_Ha_obs = readobservation("spec/RZPsc_16-11-2013_proc.dat")
    v_Hb_obs, r_Hb_obs = readobservation("spec/RZPsc_Hb_16-11-2013_proc.dat")
    Ha_plt = plot(v_Ha_obs, r_Ha_obs, legend = false)
    plot!(Ha_plt, v_Ha_mod, r_Ha_mod .+ 0.01)
    plot!(Ha_plt, v_Ha_nos, r_Ha_nos .+ 0.01)
    Hb_plt = plot(v_Hb_obs, r_Hb_obs, legend = false)
    plot!(Hb_plt, v_Hb_mod, r_Hb_mod .- 0.03)
    plot!(Hb_plt, v_Hb_nos, r_Hb_nos .- 0.03)
    plot(Ha_plt, Hb_plt, layout = @layout([A B]), size = (1400, 700), dpi = 600)
end

function plotgauss(pars, names, δ_cut; la = 0.1)
    
    plt = plot(clims = (-11,-9))
    # min_δ = min(min_stat_δ, min_nonstat_δ)
    best_names = []
    best_pars = []
    min_δ = minimum(pars[:,end])
    δ_range = δ_cut - min_δ
    alphacolor(δ) = (1.0 - (δ - min_δ)/δ_range)
    for i=1:length(names)
        model_name, profile_name = names[i]
        type = split(model_name, '_')[4]
        # if type == "stat"
        #     min_δ = min_stat_δ
        # elseif type == "nonstat"
        #     min_δ = min_nonstat_δ
        # end
        if pars[i,end] < δ_cut
            # println(names[i], " ", pars[i,5:end])
            # push!(plots, plt)
            scatter!(plt, [pars[i,6]], [pars[i,7]], marker_z=log10.(pars[i,1]), ma = (alphacolor(pars[i,end]))^(2), label = false)
            push!(best_names, [model_name, profile_name])
            push!(best_pars, [i; pars[i,:]])
        end
    end
    # best_pars = hcat(best_pars...)
    # scatter!(plt, best_pars[7,:], best_pars[8,:], marker_z=log10.(best_pars[2,:]))
    # show(best_pars)
    # for (name, pars) in zip(best_names, best_pars)
    #     println(name)
    #     println(pars)
    # end
    plt
end