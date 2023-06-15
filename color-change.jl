using TTauUtils
using Plots
using Plots.PlotMeasures
using LaTeXStrings
using BSON
using DelimitedFiles

star = Star("RYLupi")

V_filter = TTauUtils.Eclipses.V
B_filter = TTauUtils.Eclipses.B

c = TTauUtils.ConstantsInCGS.c

B_ν = c/B_filter.λ_0
V_ν = c/V_filter.λ_0

function calcΔV(f_spot, f_a, T_spot, T_star, S_f, ΔV0, scatter_f, σ_f, τ_w)
    V_filter = TTauUtils.Eclipses.V
    V_ν = c/V_filter.λ_0
    B_V_star = TTauUtils.Stars.planckfunction(V_ν, T_star)
    B_V_spot = TTauUtils.Stars.planckfunction(V_ν, T_spot)
    E_V_br = B_V_star*(1 - f_a) + B_V_spot*f_a
    E_V0_br = B_V_star*(1 - f_spot) + B_V_spot*f_spot
    eclipse_dark = 10^(-0.4*ΔV0)
    # return 2.5*log10(1 + S_f*(E_V_br/E_V0_br/eclipse_dark/(1 + scatter_f) + (eclipse_dark*(1 + scatter_f) - scatter_f)^σ_f/eclipse_dark/(1 + scatter_f)))
    return -2.5*log10(1 + S_f*E_V_br/E_V0_br/eclipse_dark/(1+scatter_f)*(exp(-τ_w) - (eclipse_dark*(1 + scatter_f) - scatter_f)^σ_f))
end

function calcΔB(f_spot, f_a, T_spot, T_star, S_f, ΔV0, ΔBV0, scatter_f, σ_f, τ_w)
    B_filter = TTauUtils.Eclipses.B
    B_ν = c/B_filter.λ_0
    B_B_star = TTauUtils.Stars.planckfunction(B_ν, T_star)
    B_B_spot = TTauUtils.Stars.planckfunction(B_ν, T_spot)
    E_B_br = B_B_star*(1 - f_a) + B_B_spot*f_a
    E_B0_br = B_B_star*(1 - f_spot) + B_B_spot*f_spot
    ΔB0 = ΔBV0 + ΔV0
    eclipse_dark = 10^(-0.4*ΔB0)
    # return 2.5*log10(1 + S_f*(E_B_br/E_B0_br/eclipse_dark/(1 + scatter_f) + (eclipse_dark*(1 + scatter_f) - scatter_f)^σ_f/eclipse_dark/(1 + scatter_f)))
    return -2.5*log10(1 + S_f*E_B_br/E_B0_br/eclipse_dark/(1+scatter_f)*(exp(-τ_w) - (eclipse_dark*(1 + scatter_f) - scatter_f)^σ_f))
end

function savedata(filename, eclipse_pars, grid, Vs, BVs)
    savedir = "color-change-saves"
    mkpath("$savedir/$filename")
    writedlm("$savedir/$filename/start.dat", eclipse_pars)
    writedlm("$savedir/$filename/grid.dat", grid)
    bson("$savedir/$filename/data.bson", Dict(:V => Vs, :BV => BVs))
end

function readdata(filename)
    savedir = "color-change-saves"
    pars = readdlm("$savedir/$filename/start.dat")
    grid_data = readdlm("$savedir/$filename/grid.dat")
    n_grid = size(grid_data)[1]
    grid = Vector{Float64}[]
    for i_grid = 1:n_grid
        ith_grid = grid_data[i_grid, :] 
        push!(grid, Float64.(ith_grid[ith_grid .!= ""]))
    end
    data = BSON.load("$savedir/$filename/data.bson")
    Vs = data[:V]
    BVs = data[:BV]
    return pars, grid, Vs, BVs
end

function findΔVsΔBVs(ΔV0, ΔBV0, T_star, T_spots, f_as_in, S_fs_in, f_spots_in, σ_fs_in, scatter_f, τ_w, scatter_vb; log_scale = true)   
    n_T_spot = length(T_spots)
    n_f_a = length(f_as_in)
    n_S_f = length(S_fs_in)
    n_f_spot = length(f_spots_in)
    n_σ_f = length(σ_fs_in)

    f_as = if log_scale
        10 .^ f_as_in
    else
        f_as_in
    end 

    S_fs = if log_scale
        10 .^ S_fs_in
    else
        S_fs_in
    end 

    f_spots = if log_scale
        10 .^ f_spots_in
    else
        f_spots_in
    end 

    σ_fs = if log_scale
        10 .^ σ_fs_in
    else
        σ_fs_in
    end 

    ΔVs = zeros(n_T_spot, n_f_a, n_S_f, n_f_spot, n_σ_f)
    ΔBVs = zeros(n_T_spot, n_f_a, n_S_f, n_f_spot, n_σ_f)

    for i_T_spot = 1:n_T_spot
        T_spot = T_spots[i_T_spot]
        println(T_spot)
        for i_f_a = 1:n_f_a
            f_a = f_as[i_f_a]
            for i_S_f = 1:n_S_f
                S_f = S_fs[i_S_f]
                for i_f_spot = 1:n_f_spot
                    f_spot = f_spots[i_f_spot]
                    if (f_a*S_f ≤ f_spot) & ((1-f_a)*S_f ≤ (1 - f_spot))
                        for i_σ_f = 1:n_σ_f
                            σ_f = σ_fs[i_σ_f]
                            try 
                                ΔV = calcΔV(f_spot, f_a, T_spot, T_star, S_f, ΔV0, scatter_f, σ_f, τ_w)
                                ΔB = calcΔB(f_spot, f_a, T_spot, T_star, S_f, ΔV0, ΔBV0, scatter_vb*scatter_f, σ_f, τ_w)
                                ΔBV = ΔB - ΔV
                                ΔVs[i_T_spot, i_f_a, i_S_f, i_f_spot, i_σ_f] = ΔV
                                ΔBVs[i_T_spot, i_f_a, i_S_f, i_f_spot, i_σ_f] = ΔBV
                            catch e
                                ΔVs[i_T_spot, i_f_a, i_S_f, i_f_spot, i_σ_f] = 1e2
                                ΔBVs[i_T_spot, i_f_a, i_S_f, i_f_spot, i_σ_f] = 1e2
                            end
                            
                        end
                    else
                        ΔVs[i_T_spot, i_f_a, i_S_f, i_f_spot, :] .= 1e2
                        ΔBVs[i_T_spot, i_f_a, i_S_f, i_f_spot, :] .= 1e2
                    end
                end
            end
        end
    end
    return ΔVs, ΔBVs
end

T_star = star.T
# T_star = 4000

V_max = 10.2
BV_min = 0.9

τ_w = 0
ΔV0 = 12.8 - V_max
ΔBV0 = 1.4 - BV_min
scatter_f = 0.07
scatter_vb = 1.0

T_spot_grid = 5e3:5e2:15e3
lgf_a_grid = -1.5:0.05:0
lgS_f_grid = -2.5:0.05:0
lgf_spot_grid = -2.5:0.05:0
lgσ_f_grid = -1:0.1:1

T_spots = collect(T_spot_grid)
# f_spot = 3e-2
lgf_as = collect(lgf_a_grid)
lgS_fs = collect(lgS_f_grid)
lgf_spots = collect(lgf_spot_grid)
lgσ_fs = collect(lgσ_f_grid)

name = "RYLup-flipC"

not_calc = try readdata(name)
    false
catch e
    true
end

ignore_data = true

ΔVs, ΔBVs = if ignore_data | not_calc
    Vs, BVs = findΔVsΔBVs(ΔV0, ΔBV0, T_star, T_spots, lgf_as, lgS_fs, lgf_spots, lgσ_fs, scatter_f, τ_w, scatter_vb, log_scale = true)
    savedata(name, [τ_w, ΔV0, ΔBV0, scatter_f, scatter_vb], 
                   [T_spots, lgf_as, lgS_fs, lgf_spots, lgσ_fs],
                    Vs, BVs)
    Vs, BVs
else
    pars, grids, Vs, BVs = readdata(name)
    τ_w, ΔV0, ΔBV0, sactter_f, scatter_vb = pars
    T_spots, lgf_as, lgS_fs, lgf_spots, lgσ_fs = grids
    Vs, BVs
end
# f_as = [0:0.02:1;]
# f_spots = [0.001:0.001:0.02;]
# S_fs = [0.001:0.001:0.02;]
# σ_fs = [0.1:0.05:1;]



ΔVobs = -0.01; δΔV = 0.05
ΔBVobs = -0.55; δΔBV = √2*0.05


begin

choose(ΔV, ΔBV) = √(((ΔV - ΔVobs)/δΔV)^2  + ((ΔBV - ΔBVobs)/δΔBV)^2)
chooseV(ΔV) = abs((ΔV - ΔVobs)/δΔV) 
chooseBV(ΔBV) = abs((ΔBV - ΔBVobs)/δΔBV) 

resid_ΔV_ΔBV = choose.(ΔVs, ΔBVs)
resid_ΔV = chooseV.(ΔVs)
resid_ΔBV = chooseBV.(ΔBVs)
i = 5

function flatten2dinds(arr, dimensions)
    if length(dimensions) != 2
        throw(ErrorException("dimensions length must be 2"))
        return arr
    end
    n_dims = ndims(arr)
    array_shape = size(arr)
    dims_to_flat = zeros(Int, n_dims - 2)
    flat_shape_arr = zeros(Int, 2)
    flat_shape_index = 1
    dims_to_flat_index = 1
    for i=1:n_dims
        if (flat_shape_index > 2)
            dims_to_flat[dims_to_flat_index] = i
            dims_to_flat_index += 1
        else
            if i == dimensions[flat_shape_index]
                flat_shape_arr[flat_shape_index] = array_shape[i]
                flat_shape_index += 1
            else
                dims_to_flat[dims_to_flat_index] = i
                dims_to_flat_index += 1
            end  
        end
    end
    flat_shape = Tuple(flat_shape_arr)
    min_arr_inds = findmin(arr, dims = dims_to_flat)[2]
    return min_arr_inds, flat_shape
end

function flatten2d(arr, dimensions, min_arr)
    min_arr_inds, flat_shape = flatten2dinds(min_arr, dimensions)
    reshape(arr[min_arr_inds], flat_shape)
end

flatten2d(arr, dimensions) = flatten2d(arr, dimensions, arr)

par_names = [L"T_\mathrm{sp}", L"\lg f_a", L"\lg S_f", L"\lg f", L"\lg σ_a"] 
par_vals = [T_spots,  lgf_as, lgS_fs, lgf_spots, lgσ_fs]

# par_names = [L"T_{sp}", L"f_w", L"S_f", L"f", L"σ_f"] 
# par_vals = [T_spots,  f_as, S_fs, f_spots, σ_fs]

color_lims = (0, 4)
flat_dims = (1,4)
flat_par_names = par_names[collect(flat_dims)]
flat_par_vals = par_vals[collect(flat_dims)]
flat_inds, flat_shape = flatten2dinds(resid_ΔV_ΔBV, flat_dims)
resid_ΔV_ΔBV_flat = reshape(resid_ΔV_ΔBV[flat_inds], flat_shape)
resid_ΔV_flat = reshape(resid_ΔV[flat_inds], flat_shape)
resid_ΔBV_flat = reshape(resid_ΔBV[flat_inds], flat_shape)
ΔVs_flat = reshape(ΔVs[flat_inds], flat_shape)
ΔBVs_flat = reshape(ΔBVs[flat_inds], flat_shape)

plt_x = flat_par_vals[2]
plt_y = flat_par_vals[1]
x_label = flat_par_names[2]
y_label = flat_par_names[1]
residplt_spot = heatmap(plt_x, plt_y, resid_ΔV_ΔBV_flat, clims = color_lims, c = cgrad([:black, :white]))
residVplt_spot = heatmap(plt_x, plt_y, resid_ΔV_flat, clims = color_lims, colorbar = :false)
residBVplt_spot = heatmap(plt_x, plt_y, resid_ΔBV_flat, clims = color_lims, colorbar = :false)
Vplt_spot = heatmap(plt_x, plt_y, ΔVs_flat)
BVplt_spot = heatmap(plt_x, plt_y, ΔBVs_flat)

# plot(residVplt, residBVplt, residplt, layout = layout = grid(1, 3, widths=[0.3 ,0.3, 0.4]), size = (1300,400))
plot!(residplt_spot, xlabel = x_label, ylabel = y_label)

flat_dims = (2,3)
flat_par_names = par_names[collect(flat_dims)]
flat_par_vals = par_vals[collect(flat_dims)]
flat_inds, flat_shape = flatten2dinds(resid_ΔV_ΔBV, flat_dims)
resid_ΔV_ΔBV_flat = reshape(resid_ΔV_ΔBV[flat_inds], flat_shape)
resid_ΔV_flat = reshape(resid_ΔV[flat_inds], flat_shape)
resid_ΔBV_flat = reshape(resid_ΔBV[flat_inds], flat_shape)
ΔVs_flat = reshape(ΔVs[flat_inds], flat_shape)
ΔBVs_flat = reshape(ΔBVs[flat_inds], flat_shape)

plt_x = flat_par_vals[2]
plt_y = flat_par_vals[1]
x_label = flat_par_names[2]
y_label = flat_par_names[1]
residplt_w = heatmap(plt_x, plt_y, resid_ΔV_ΔBV_flat, clims = color_lims, c = cgrad([:black, :white]))
residVplt_w = heatmap(plt_x, plt_y, resid_ΔV_flat, clims = color_lims, colorbar = :false)
residBVplt_w = heatmap(plt_x, plt_y, resid_ΔBV_flat, clims = color_lims, colorbar = :false)
Vplt = heatmap(plt_x, plt_y, ΔVs_flat)
BVplt = heatmap(plt_x, plt_y, ΔBVs_flat)

# plot(residVplt, residBVplt, residplt, layout = layout = grid(1, 3, widths=[0.3 ,0.3, 0.4]), size = (1300,400))
plot!(residplt_w, xlabel = x_label, ylabel = y_label)

flat_dims = (3,5)
flat_par_names = par_names[collect(flat_dims)]
flat_par_vals = par_vals[collect(flat_dims)]
flat_inds, flat_shape = flatten2dinds(resid_ΔV_ΔBV, flat_dims)
resid_ΔV_ΔBV_flat = reshape(resid_ΔV_ΔBV[flat_inds], flat_shape)
resid_ΔV_flat = reshape(resid_ΔV[flat_inds], flat_shape)
resid_ΔBV_flat = reshape(resid_ΔBV[flat_inds], flat_shape)
ΔVs_flat = reshape(ΔVs[flat_inds], flat_shape)
ΔBVs_flat = reshape(ΔBVs[flat_inds], flat_shape)

plt_x = flat_par_vals[2]
plt_y = flat_par_vals[1]
y_label = flat_par_names[2]
x_label = flat_par_names[1]
residplt_sig = heatmap(plt_y, plt_x, resid_ΔV_ΔBV_flat', clims = color_lims, c = cgrad([:black, :white]))
residVplt_sig = heatmap(plt_x, plt_y, resid_ΔV_flat, clims = color_lims, colorbar = :false)
residBVplt_sig = heatmap(plt_x, plt_y, resid_ΔBV_flat, clims = color_lims, colorbar = :false)
Vplt = heatmap(plt_x, plt_y, ΔVs_flat)
BVplt = heatmap(plt_x, plt_y, ΔBVs_flat)

# plot(residVplt, residBVplt, residplt, layout = layout = grid(1, 3, widths=[0.3 ,0.3, 0.4]), size = (1300,400))
plot!(residplt_sig, xlabel = x_label, ylabel = y_label)

flat_dims = (3,4)
flat_par_names = par_names[collect(flat_dims)]
flat_par_vals = par_vals[collect(flat_dims)]
flat_inds, flat_shape = flatten2dinds(resid_ΔV_ΔBV, flat_dims)
resid_ΔV_ΔBV_flat = reshape(resid_ΔV_ΔBV[flat_inds], flat_shape)
resid_ΔV_flat = reshape(resid_ΔV[flat_inds], flat_shape)
resid_ΔBV_flat = reshape(resid_ΔBV[flat_inds], flat_shape)
ΔVs_flat = reshape(ΔVs[flat_inds], flat_shape)
ΔBVs_flat = reshape(ΔBVs[flat_inds], flat_shape)

plt_x = flat_par_vals[2]
plt_y = flat_par_vals[1]
x_label = flat_par_names[2]
y_label = flat_par_names[1]
residplt_sf = heatmap(plt_x, plt_y, resid_ΔV_ΔBV_flat, clims = color_lims, c = cgrad([:black, :white]))
residVplt_sf = heatmap(plt_x, plt_y, resid_ΔV_flat, clims = color_lims, colorbar = :false)
residBVplt_sf = heatmap(plt_x, plt_y, resid_ΔBV_flat, clims = color_lims, colorbar = :false)
Vplt = heatmap(plt_x, plt_y, ΔVs_flat)
BVplt = heatmap(plt_x, plt_y, ΔBVs_flat)

# plot(residVplt, residBVplt, residplt, layout = layout = grid(1, 3, widths=[0.3 ,0.3, 0.4]), size = (1300,400))
plot!(residplt_sf, xlabel = x_label, ylabel = y_label)

allresidplt = plot(residplt_spot, residplt_w, residplt_sf, residplt_sig, layout = grid(2,2), size = (1200, 600),
                     plot_title = L"\Delta V = %$ΔVobs,\ \Delta (B-V) = %$ΔBVobs,\ \tau_a = %$τ_w", leftmargin = 20px, bottommargin = 20px)

allresidplt
end

begin
    V_err = δΔV
    BV_err = δΔBV

    target_V = ΔVobs
    target_BV = ΔBVobs
    hist_dataV = vec(ΔVs[:,:,:,:,:])
    hist_dataBV = vec(ΔBVs[:,:,:,:,:])

    normal_V_data = hist_dataV[hist_dataV .< 50]
    normal_BV_data = hist_dataBV[hist_dataBV .< 50]

    box_width = 0.05
    
    start_val, end_val = if ΔVobs > 0
        0.0, 5.0
    else
        -5.0, 0.0
    end

    xlims, ylims = if ΔVobs > 0
        (0,1.2), (0,1)
    else
        (-1.2,0), (-1,0)
    end

    V_boxes = [start_val:box_width:end_val-box_width;] .+ box_width/2
    BV_boxes = [start_val:box_width:end_val-box_width;] .+ box_width/2
    n_V_boxes = length(V_boxes)
    n_BV_boxes = length(BV_boxes)

    function inboxcount2d(x_boxes, y_boxes, x_array, y_array)
        n_x_box = length(x_boxes)
        n_y_box = length(y_boxes)
        counts = zeros(Int, (n_x_box, n_y_box))
        x_box_width = x_boxes[2] - x_boxes[1]
        y_box_width = y_boxes[2] - y_boxes[1]
        start_x = x_boxes[1] - x_box_width/2
        start_y = y_boxes[1] - y_box_width/2
        println(start_x, " ", start_y)
        n_array = length(x_array)
        for i_arr = 1:n_array
            x = x_array[i_arr]
            y = y_array[i_arr]
            i_box_x = ceil(Int, (x-start_x)/x_box_width)
            i_box_y = ceil(Int, (y-start_y)/y_box_width)
            if (i_box_x > n_x_box) | (i_box_y > n_y_box) | (i_box_x < 1) | (i_box_y < 1)
                continue
            end
            counts[i_box_x, i_box_y] += 1
        end
        return counts
    end

    counts = inboxcount2d(V_boxes, BV_boxes, normal_V_data, normal_BV_data) 

    # scatter(normal_V_data, normal_BV_data, xlims = (0,1), ylims = (0,1), xlabel = L"\Delta V", ylabel = L"\Delta (B-V)", ms = 0.1, aspect_ratio = :equal)
    hist_plt = heatmap(BV_boxes, V_boxes, log10.(counts), aspect_ratio = :equal, c = cgrad([:white, :gray20]), xlims = xlims, ylims = ylims,
            clims = (0, 5), ylabel = L"\Delta V", xlabel = L"\Delta (B-V)", colorbar_title = L"\log N", xticks = -1.2:0.2:0)
    scatter!(hist_plt, [target_BV], [target_V], xerr = [BV_err], yerr = [V_err], lc = :black, mc = :black)
end

savefig(hist_plt, "color-change-saves/$name-hist.pdf")
savefig(allresidplt, "color-change-saves/$name-resid.pdf")