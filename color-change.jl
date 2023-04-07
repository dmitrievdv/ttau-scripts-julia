using TTauUtils
using Plots
using LaTeXStrings

star = Star("RYLupi")

T_star = star.T

T_star = 4000

V_filter = TTauUtils.Eclipses.V
B_filter = TTauUtils.Eclipses.B

c = TTauUtils.ConstantsInCGS.c

B_ν = c/B_filter.λ_0
V_ν = c/V_filter.λ_0

function calcΔV(f_spot, f_a, T_spot, T_star, S_f, ΔV0, scatter_f, σ_f)
    V_filter = TTauUtils.Eclipses.V
    V_ν = c/V_filter.λ_0
    B_V_star = TTauUtils.Stars.planckfunction(V_ν, T_star)
    B_V_spot = TTauUtils.Stars.planckfunction(V_ν, T_spot)
    E_V_br = B_V_star*(1 - f_a) + B_V_spot*f_a
    E_V0_br = B_V_star*(1 - f_spot) + B_V_spot*f_spot
    eclipse_dark = 10^(-0.4*ΔV0)*(1 + scatter_f)
    return 2.5*log10(1 + S_f*(E_V_br/E_V0_br/eclipse_dark + (eclipse_dark - scatter_f)^σ_f/eclipse_dark))
end

function calcΔB(f_spot, f_a, T_spot, T_star, S_f, ΔV0, ΔBV0, scatter_f, σ_f)
    B_filter = TTauUtils.Eclipses.B
    B_ν = c/B_filter.λ_0
    B_B_star = TTauUtils.Stars.planckfunction(B_ν, T_star)
    B_B_spot = TTauUtils.Stars.planckfunction(B_ν, T_spot)
    E_B_br = B_B_star*(1 - f_a) + B_B_spot*f_a
    E_B0_br = B_B_star*(1 - f_spot) + B_B_spot*f_spot
    ΔB0 = ΔBV0 + ΔV0
    eclipse_dark = 10^(-0.4*ΔB0)*(1 + scatter_f)
    return 2.5*log10(1 + S_f*(E_B_br/E_B0_br/eclipse_dark + (eclipse_dark - scatter_f)^σ_f/eclipse_dark))
end

ΔV0 = 2
ΔBV0 = 0.3
T_spots = [5e3:5e2:20e3;]
# f_spot = 3e-2
lgf_as = [-2:0.05:0;]
lgS_fs = [-3:0.05:0;]
lgf_spots = [-4:0.1:0;]
lgσ_fs = [-1:0.2:1;]
# scatter_f = 0.06

function findΔVsΔBVs(ΔV0, ΔBV0, T_star, T_spots, lgf_as, lgS_fs, lgf_spots, lgσ_fs, scatter_f)   
    n_T_spot = length(T_spots)
    n_f_a = length(lgf_as)
    n_S_f = length(lgS_fs)
    n_f_spot = length(lgf_spots)
    n_σ_f = length(lgσ_fs)

    ΔVs = zeros(n_T_spot, n_f_a, n_S_f, n_f_spot, n_σ_f)
    ΔBVs = zeros(n_T_spot, n_f_a, n_S_f, n_f_spot, n_σ_f)

    for i_T_spot = 1:n_T_spot
        T_spot = T_spots[i_T_spot]
        for i_f_a = 1:n_f_a
            f_a = 10^lgf_as[i_f_a]
            for i_S_f = 1:n_S_f
                S_f = 10^lgS_fs[i_S_f]
                for i_f_spot = 1:n_f_spot
                    f_spot = 10^lgf_spots[i_f_spot]
                    for i_σ_f = 1:n_σ_f
                        σ_f = 10^lgσ_fs[i_σ_f]
                        ΔV = calcΔV(f_spot, f_a, T_spot, T_star, S_f, ΔV0, 0.06, σ_f)
                        ΔB = calcΔB(f_spot, f_a, T_spot, T_star, S_f, ΔV0, ΔBV0, 0.06, σ_f)
                        ΔBV = ΔB - ΔV
                        ΔVs[i_T_spot, i_f_a, i_S_f, i_f_spot, i_σ_f] = ΔV
                        ΔBVs[i_T_spot, i_f_a, i_S_f, i_f_spot, i_σ_f] = ΔBV
                    end
                end
            end
        end
    end
    return ΔVs, ΔBVs
end

ΔVs, ΔBVs = findΔVsΔBVs(ΔV0, ΔBV0, T_star, T_spots, lgf_as, lgS_fs, lgf_spots, lgσ_fs, scatter_f)
begin
ΔVobs = 0.2; δΔV = 0.1
ΔBVobs = 0.3; δΔBV = 0.1
choose(ΔV, ΔBV) = √(((ΔV - ΔVobs)/δΔV)^2  + ((ΔBV - ΔBVobs)/δΔBV)^2)
chooseV(ΔV) = abs((ΔV - ΔVobs)/δΔV) 
chooseBV(ΔBV) = abs((ΔBV - ΔBVobs)/δΔBV) 

resid_ΔV_ΔBV = choose.(ΔVs, ΔBVs)
resid_ΔV = chooseV.(ΔVs)
resid_ΔBV = chooseBV.(ΔBVs)
i = 5

function flatten2d(arr, dimensions)
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
    min_arr = minimum(arr, dims = dims_to_flat)
    reshape(min_arr, flat_shape)
end

par_names = [L"T_{sp}", L"\lg f_w", L"\lg S_f", L"\lg f", L"\lg σ_f"] 
par_vals = [T_spots,  lgf_as, lgS_fs, lgf_spots, lgσ_fs]

flat_dims = (1,4)
flat_par_names = par_names[collect(flat_dims)]
flat_par_vals = par_vals[collect(flat_dims)]
resid_ΔV_ΔBV_flat = flatten2d(resid_ΔV_ΔBV, flat_dims)
resid_ΔV_flat = flatten2d(resid_ΔV, flat_dims)
resid_ΔBV_flat = flatten2d(resid_ΔBV, flat_dims)
ΔVs_flat = flatten2d(ΔVs, flat_dims)
ΔBVs_flat = flatten2d(ΔBVs, flat_dims)

plt_x = flat_par_vals[2]
plt_y = flat_par_vals[1]
x_label = flat_par_names[2]
y_label = flat_par_names[1]
residplt_spot = heatmap(plt_x, plt_y, resid_ΔV_ΔBV_flat, clims = (0,5))
residVplt_spot = heatmap(plt_x, plt_y, resid_ΔV_flat, clims = (0,5), colorbar = :false)
residBVplt_spot = heatmap(plt_x, plt_y, resid_ΔBV_flat, clims = (0,5), colorbar = :false)
Vplt_spot = heatmap(plt_x, plt_y, ΔVs_flat)
BVplt_spot = heatmap(plt_x, plt_y, ΔBVs_flat)

# plot(residVplt, residBVplt, residplt, layout = layout = grid(1, 3, widths=[0.3 ,0.3, 0.4]), size = (1300,400))
plot!(residplt_spot, xlabel = x_label, ylabel = y_label)

flat_dims = (2,3)
flat_par_names = par_names[collect(flat_dims)]
flat_par_vals = par_vals[collect(flat_dims)]
resid_ΔV_ΔBV_flat = flatten2d(resid_ΔV_ΔBV, flat_dims)
resid_ΔV_flat = flatten2d(resid_ΔV, flat_dims)
resid_ΔBV_flat = flatten2d(resid_ΔBV, flat_dims)
ΔVs_flat = flatten2d(ΔVs, flat_dims)
ΔBVs_flat = flatten2d(ΔBVs, flat_dims)

plt_x = flat_par_vals[2]
plt_y = flat_par_vals[1]
x_label = flat_par_names[2]
y_label = flat_par_names[1]
residplt_w = heatmap(plt_x, plt_y, resid_ΔV_ΔBV_flat, clims = (0,5))
residVplt_w = heatmap(plt_x, plt_y, resid_ΔV_flat, clims = (0,5), colorbar = :false)
residBVplt_w = heatmap(plt_x, plt_y, resid_ΔBV_flat, clims = (0,5), colorbar = :false)
Vplt = heatmap(plt_x, plt_y, ΔVs_flat)
BVplt = heatmap(plt_x, plt_y, ΔBVs_flat)

# plot(residVplt, residBVplt, residplt, layout = layout = grid(1, 3, widths=[0.3 ,0.3, 0.4]), size = (1300,400))
plot!(residplt_w, xlabel = x_label, ylabel = y_label)

flat_dims = (3,5)
flat_par_names = par_names[collect(flat_dims)]
flat_par_vals = par_vals[collect(flat_dims)]
resid_ΔV_ΔBV_flat = flatten2d(resid_ΔV_ΔBV, flat_dims)
resid_ΔV_flat = flatten2d(resid_ΔV, flat_dims)
resid_ΔBV_flat = flatten2d(resid_ΔBV, flat_dims)
ΔVs_flat = flatten2d(ΔVs, flat_dims)
ΔBVs_flat = flatten2d(ΔBVs, flat_dims)

plt_x = flat_par_vals[2]
plt_y = flat_par_vals[1]
x_label = flat_par_names[2]
y_label = flat_par_names[1]
residplt_sig = heatmap(plt_x, plt_y, resid_ΔV_ΔBV_flat, clims = (0,5))
residVplt_sig = heatmap(plt_x, plt_y, resid_ΔV_flat, clims = (0,5), colorbar = :false)
residBVplt_sig = heatmap(plt_x, plt_y, resid_ΔBV_flat, clims = (0,5), colorbar = :false)
Vplt = heatmap(plt_x, plt_y, ΔVs_flat)
BVplt = heatmap(plt_x, plt_y, ΔBVs_flat)

# plot(residVplt, residBVplt, residplt, layout = layout = grid(1, 3, widths=[0.3 ,0.3, 0.4]), size = (1300,400))
plot!(residplt_sig, xlabel = x_label, ylabel = y_label)

plot(residplt_spot, residplt_w, layout = grid(2,1), size = (600, 600))
end