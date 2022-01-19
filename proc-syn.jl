using Printf
using TTauUtils

star = Star("RZPsc")
ν_Hb = TTauUtils.HydrogenPopulations.linefrequency(3,2)
λ_Hb = 2.998e18/ν_Hb
I_c = TTauUtils.starcontinuum(star, ν_Hb)
I_c_syn = I_c*π*2.998e18/λ_Hb^2

obs_file = "spec/RZ_Psc_Hb_syn_unwid.txt"

λ_obs_data = Float64[]
r_obs_data = Float64[]

lines = readlines(obs_file)

for line in lines
    try
        λ, r = parse.(Float64, split(line))
        println("$λ $r")
        push!(λ_obs_data, λ)
        push!(r_obs_data, r)
    catch
        continue
    end 
end

# v_obs_data = @. (λ_obs_data - 6562.8)/6562.8*2.998e5

# v_obs = v_obs_data[(v_obs_data .> -600) .& (v_obs_data .< 600)]
# r_obs = r_obs_data[(v_obs_data .> -600) .& (v_obs_data .< 600)]

# plot(v_obs, r_obs)


out = open("spec/RZ_Psc_Hb_syn_unwid_corr.dat", "w")

for (λ,r) in zip(λ_obs_data, r_obs_data)
    @printf(out, "%.4f %.5e %.5f\n", λ, r*I_c_syn, r)
end
close(out)