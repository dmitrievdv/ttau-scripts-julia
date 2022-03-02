using TTauUtils
using Plots

function getvandr(profile :: TTauUtils.HydrogenProfile)
    u = profile.upper_level; l = profile.lower_level
    ν_0 = TTauUtils.HydrogenPopulations.linefrequency(u,l)
    v_prof = (ν_0 .- profile.frequencies) / ν_0*3e5 
    r_prof = profile.profile
    if v_prof[1] > 0.0
        reverse!(v_prof)
        reverse!(r_prof)
    end
    return v_prof, r_prof
end

model_name = "99_10500_5-9_nonstat_nonlocal"
star = Star("RZPsc")

model = loadmodel(star, model_name)

i = 30; α = 15

profiles = HydrogenProfile[]
spec_profiles = HydrogenProfile[]

for ϕ = 0:10:180
    profile = TTauUtils.HydrogenProfileDoppler(model, 3, 2, Orientation(i, α, ϕ), 0.1, 0.1, 0.1, 100, blue_v_max = 300, red_v_max = 600)
    saveprofile(profile, "Ha_rot_$i-$α-$ϕ")
    push!(profiles, profile)
    prof_spec = TTauUtils.addphotosphespecdoppler(profile, 0.3, "spec/RZ_Psc_Ha_syn_unwid_corr.dat")
    saveprofile(prof_spec, "Ha_rot_$(i)-$(α)-$(ϕ)_phot3crude")
    push!(spec_profiles, prof_spec)
end

plt = plot()
for prof in spec_profiles
    v, r = getvandr(prof)
    plot!(plt, v, r, label = "$(round(Int, prof.orientation.ϕ/pi*180))")
end
plt