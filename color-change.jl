using TTauUtils
using Plots
using Plots.PlotMeasures
using LaTeXStrings

star = Star("RYLupi")

T_star = star.T

# T_star = 4000

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

ΔV0 = 1.5
ΔBV0 = 0.3
T_spots = [5e3:5e2:15e3;]
# f_spot = 3e-2
lgf_as = [-1.5:0.05:0;]
lgS_fs = [-2.5:0.05:0;]
lgf_spots = [-2.5:0.05:0;]
lgσ_fs = [-1:0.1:1;]
scatter_f = 0.06
scatter_vb = 1.3

# f_as = [0:0.02:1;]
# f_spots = [0.001:0.001:0.02;]
# S_fs = [0.001:0.001:0.02;]
# σ_fs = [0.1:0.05:1;]

τ_w = 10

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
                                ΔVs[i_T_spot, i_f_a, i_S_f, i_f_spot, i_σ_f] = 1e-2
                                ΔBVs[i_T_spot, i_f_a, i_S_f, i_f_spot, i_σ_f] = 1e-2
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

ΔVs, ΔBVs = findΔVsΔBVs(ΔV0, ΔBV0, T_star, T_spots, lgf_as, lgS_fs, lgf_spots, lgσ_fs, scatter_f, τ_w, scatter_vb, log_scale = true)
begin
ΔVobs = 0.3; δΔV = 1
ΔBVobs = 0.3; δΔBV = 1
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

color_lims = (0, 0.25)
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
residplt_spot = heatmap(plt_x, plt_y, resid_ΔV_ΔBV_flat, clims = color_lims)
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
residplt_w = heatmap(plt_x, plt_y, resid_ΔV_ΔBV_flat, clims = color_lims)
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
residplt_sig = heatmap(plt_y, plt_x, resid_ΔV_ΔBV_flat', clims = color_lims)
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
residplt_sf = heatmap(plt_x, plt_y, resid_ΔV_ΔBV_flat, clims = color_lims)
residVplt_sf = heatmap(plt_x, plt_y, resid_ΔV_flat, clims = color_lims, colorbar = :false)
residBVplt_sf = heatmap(plt_x, plt_y, resid_ΔBV_flat, clims = color_lims, colorbar = :false)
Vplt = heatmap(plt_x, plt_y, ΔVs_flat)
BVplt = heatmap(plt_x, plt_y, ΔBVs_flat)

# plot(residVplt, residBVplt, residplt, layout = layout = grid(1, 3, widths=[0.3 ,0.3, 0.4]), size = (1300,400))
plot!(residplt_sf, xlabel = x_label, ylabel = y_label)

allresidplt = plot(residplt_spot, residplt_w, residplt_sf, residplt_sig, layout = grid(2,2), size = (1200, 600),
                     plot_title = L"\Delta V = %$ΔVobs,\ \Delta (B-V) = %$ΔBVobs,\ \tau = %$τ_w", leftmargin = 20px, bottommargin = 20px)

normal_data = ΔVs .< 50

hist_dataV = vec(ΔVs[normal_data])
hist_dataBV = vec(ΔBVs[normal_data])

allresidplt
end