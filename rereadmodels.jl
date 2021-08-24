using TTauUtils

function oldSolidMagnetosphere(star, name; dir = "stars")
    file_dir = dir*'/'*star.name*'/'*name
    oldSolidMagnetosphere(name, dir = file_dir, file_ext = "dat")
end

function oldSolidMagnetosphere(name; dir = ".", file_ext = "model")
    lines = readlines(dir*'/'*name*"."*file_ext)    
    if split(lines[1])[2] != "SolidMagnetosphere"; throw(WrongModelType); end
    alg = split(lines[1])[5]
    star = Star(split(lines[3])[3])
    r_mi, r_mo = parse.(Float64, split(lines[4])[4:5])
    geometry = TTauUtils.GeometryAndOrientations.DipoleGeometry(r_mi, r_mo)
    M_dot = parse(Float64, split(lines[5])[5])
    v_start = parse(Float64, split(lines[6])[5])
    T_max = parse(Float64, split(lines[7])[8])
    kinematics = TTauUtils.Models.SolidDipoleKin(TTauUtils.Stars.starradiiincm(star), TTauUtils.Stars.starvescincms(star), v_start*1e5)

    T_spot = TTauUtils.Stars.calcmagspottemperature(star, M_dot, r_mi, r_mo)
    spotted_star = TTauUtils.Stars.MagnetosphereSpotStar(star, T_spot, r_mi, r_mo)
    
    r_m_grid = parse.(Float64, split(lines[10]))
    t_grid = parse.(Float64, split(lines[13]))[1:end-1]
    # t_grid[end] = 1.0

    n_levels = parse(Int, split(lines[15])[end])

    n_r_m = length(r_m_grid); n_t = length(t_grid)
    Te_table = zeros(n_r_m, n_t)
    nh_table = zeros(n_r_m, n_t)
    ne_table = zeros(n_r_m, n_t)
    ni_tables = zeros(n_r_m, n_t, n_levels)


    for i ∈ 1:n_r_m
        for j ∈ 1:n_t
            line_n = 16 + i + (i-1)*(n_t+1) + (j-1)
            line_data = parse.(Float64, split(lines[line_n]))
            Te_table[i,j] = line_data[4]
            nh_table[i,j] = line_data[5]
            ne_table[i,j] = line_data[6]
            for l=1:n_levels
                ni_tables[i,j,l] = line_data[6+l]
            end
        end
    end
    
    Te_spl2d, lgnh_spl2d, lgne_spl2d, lgni_spl2d = TTauUtils.Models.dipolemagnetospheresplines(r_m_grid, t_grid, Te_table, nh_table, ne_table, ni_tables)

    SolidMagnetosphere(name, alg, M_dot, T_max, v_start, geometry, kinematics, spotted_star, 
                         r_m_grid, t_grid, lgnh_spl2d, lgne_spl2d, Te_spl2d, lgni_spl2d)
end

star_name = "RZPsc"
star = Star(star_name)
models = readdir("stars/$star_name")[2:end]

for model_name in models
    println(model_name)
    model = oldSolidMagnetosphere(star, model_name)
    savemodel(model)
end