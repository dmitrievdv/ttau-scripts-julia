using TTauUtils
using Plots

star = Star("RYLupi")

T_star = star.T

V_filter = TTauUtils.Eclipses.V
B_filter = TTauUtils.Eclipses.B

c = TTauUtils.ConstantsInCGS.c

B_ν = c/B_filter.λ_0
V_ν = c/V_filter.λ_0

function calcΔV(f_spot, f_a, T_spot, T_star, S_f, ΔV0)
    V_filter = TTauUtils.Eclipses.V
    V_ν = c/V_filter.λ_0
    B_V_star = TTauUtils.Stars.planckfunction(V_ν, T_star)
    B_V_spot = TTauUtils.Stars.planckfunction(V_ν, T_spot)
    E_V_br = B_V_star*(1 - f_a) + B_V_spot*f_a
    E_V0_br = B_V_star*(1 - f_spot) + B_V_spot*f_spot
    return 2.5*log10(1 + E_V_br/E_V0_br*S_f*10^(0.4*ΔV0))
end

function calcΔB(f_spot, f_a, T_spot, T_star, S_f, ΔV0, ΔBV0)
    B_filter = TTauUtils.Eclipses.B
    B_ν = c/B_filter.λ_0
    B_B_star = TTauUtils.Stars.planckfunction(B_ν, T_star)
    B_B_spot = TTauUtils.Stars.planckfunction(B_ν, T_spot)
    E_B_br = B_B_star*(1 - f_a) + B_B_spot*f_a
    E_B0_br = B_B_star*(1 - f_spot) + B_B_spot*f_spot
    ΔB0 = ΔBV0 + ΔV0
    return 2.5*log10(1 + E_B_br/E_B0_br*S_f*10^(0.4*ΔB0))
end

ΔV0 = 2
ΔBV0 = 0.4

T_spots = [5e3:5e2:16e3;]
f_spot = 3e-2
f_as = [0.0:0.1:1;]
S_fs = [0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1]

n_T_spot = length(T_spots)
n_f_a = length(f_as)
n_S_f = length(S_fs)

ΔVs = zeros(n_T_spot, n_f_a, n_S_f)
ΔBVs = zeros(n_T_spot, n_f_a, n_S_f)

for i_T_spot = 1:n_T_spot
    T_spot = T_spots[i_T_spot]
    for i_f_a = 1:n_f_a
        f_a = f_as[i_f_a]
        for i_S_f = 1:n_S_f
            S_f = S_fs[i_S_f]
            ΔV = calcΔV(f_spot, f_a, T_spot, T_star, S_f, ΔV0)
            ΔB = calcΔB(f_spot, f_a, T_spot, T_star, S_f, ΔV0, ΔBV0)
            ΔBV = ΔB - ΔV
            ΔVs[i_T_spot, i_f_a, i_S_f] = ΔV
            ΔBVs[i_T_spot, i_f_a, i_S_f] = ΔBV
        end
    end
end

begin
ΔVobs = 0.0; δΔV = 0.2
ΔBVobs = 0.7; δΔBV = 0.1
choose(ΔV, ΔBV) = √(((ΔV - ΔVobs)/δΔV)^2  + ((ΔBV - ΔBVobs)/δΔBV)^2)
chooseV(ΔV) = abs((ΔV - ΔVobs)/δΔV) 
chooseBV(ΔBV) = abs((ΔBV - ΔBVobs)/δΔBV) 

resid_ΔV_ΔBV = choose.(ΔVs, ΔBVs)
resid_ΔV = chooseV.(ΔVs)
resid_ΔBV = chooseBV.(ΔBVs)
i = 5

residplt = heatmap(resid_ΔV_ΔBV[:,:,i], clims = (0,10))
residVplt = heatmap(resid_ΔV[:,:,i], clims = (0,10))
residBVplt = heatmap(resid_ΔBV[:,:,i], clims = (0,10))
Vplt = heatmap(ΔVs[:,:,i])
BVplt = heatmap(ΔBVs[:,:,i])

plot(residVplt, residBVplt, Vplt, BVplt, layout = @layout([C D; A B]), size = (700,700))
end