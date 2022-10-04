using TTauUtils

models_file = "checmodels.dat"
model_names = readlines(models_file)
star = Star("RZPsc")

for model_name in model_names
    model_old = loadmodel(star, model_name)
    alg = model_old.alg
    r_mi, r_mo = model_old.geometry.r_mi, model_old.geometry.r_mo
    Ṁ = model_old.Mdot
    T_max = model_old.T_max

    if alg == "HartNHCoolStat"
        model_new = StationarySolidMagnetosphere(model_name[1:end-9], star, r_mi, r_mo, Ṁ, T_max, 10, n_t = 20, progress_output = false)
    elseif alg == "HartNHCoolNonStat"
        model_new = NonStationarySolidMagnetosphere(model_name[1:end-9], star, r_mi, r_mo, Ṁ, T_max, 10, n_t = 20, progress_output = false)
    end
    model_new = addnonlocal(model_new)
    savemodel(model_new)
end