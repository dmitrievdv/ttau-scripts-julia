using Printf

obs_file = "spec/2013_10_21_Ha_helioV2.txt"

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

v_obs_data = @. (λ_obs_data - 6562.8)/6562.8*2.998e5

v_obs = v_obs_data[(v_obs_data .> -600) .& (v_obs_data .< 600)]
r_obs = r_obs_data[(v_obs_data .> -600) .& (v_obs_data .< 600)]

plot(v_obs, r_obs)


out = open("spec/RZPsc_21-10-2013_proc.dat", "w")

for (v,r) in zip(v_obs, r_obs)
    @printf(out, "%.4f %.5f\n", v, r)
end
close(out)