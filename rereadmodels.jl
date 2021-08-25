using TTauUtils
using Plots.Measures

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

star_name = "hart94"
star = Star(star_name)
models = readdir("stars/$star_name")[2:end]

for model_name in models
    println(model_name)
    model = oldSolidMagnetosphere(star, model_name)
    savemodel(model)
end

# begin 
#     star_name = "hart94"
#     star = Star(star_name)
#     model_name = "hart94_95_8500_2-3_stat_nonlocal"
#     model = SolidMagnetosphere(star, model_name)
#     r_ms = model.r_m_grid
#     t_points = model.t_grid
#     t_lines = collect(range(0, 1, length = 100))

#     T_e_plot = plot()
#     n_h_plot = plot()
#     f_plot = plot()

#     for r_m in r_ms
#         θ_lines = @. asin(√(1/r_m)) + t_lines * (π/2 - asin(√(1/r_m)))
#         r_lines = @. r_m*sin(θ_lines)^2
#         θ_points = @. asin(√(1/r_m)) + t_points * (π/2 - asin(√(1/r_m)))
#         r_points = @. r_m*sin(θ_points)^2
#         n_e_lines = @. 10^model.lgne_spl2d.(r_m, ts_lines)
#         n_h_lines = @. 10^model.lgnh_spl2d.(r_m, ts_lines)
#         T_e_lines = @. model.Te_spl2d.(r_m, ts_lines)

#         n_e_points = @. 10^model.lgne_spl2d.(r_m, ts_points)
#         n_h_points = @. 10^model.lgnh_spl2d.(r_m, ts_points)
#         T_e_points = @. model.Te_spl2d.(r_m, ts_points)

#         plot!(T_e_plot, r_lines, T_e_lines, lc = :black, legend = false, ylabel = L"T_e,\ \mathrm{K}", xlabel = L"r,\ \mathrm{R}_\star", margin = 5mm)
#         scatter!(T_e_plot, r_points, T_e_points, mc = :black, ms = 2, legend = false)
#         plot!(n_h_plot, r_lines, (@.  log10(n_h_lines) ), lc = :black, legend = false, ylabel = L"\lg\ n_H,\ \mathrm{cm}^{-3}", xlabel = L"r,\ \mathrm{R}_\star", margin = 5mm)
#         scatter!(n_h_plot, r_points, (@.  log10(n_h_points)), mc = :black, ms = 2, legend = false)
#         plot!(f_plot, r_lines, (@.  n_e_lines/n_h_lines), lc = :black, legend = false, ylabel = L"n_e/n_H", xlabel = L"r,\ \mathrm{R}_\star")
#         scatter!(f_plot, r_points, (@.  n_e_points/n_h_points), mc = :black, ms = 2, legend = false)
#     end
#     plt = plot(T_e_plot, n_h_plot, layout = (@layout grid(1,2)), size = (1000, 500), dpi = 300)
#     savefig(plt, "Te_nH.png")
# end

# begin 
#     star_name = "hart94"
#     star = Star(star_name)
#     model_name = "hart94_90_8000_2-3_stat_nonlocal"
#     model = SolidMagnetosphere(star, model_name)
#     r_ms = model.r_m_grid
#     t_points = model.t_grid
#     t_lines = collect(range(0, 1, length = 100))

#     T_e_plot = plot()
#     n_h_plot = plot()
#     f_plot = plot()

#     for r_m in r_ms
#         θ_lines = @. asin(√(1/r_m)) + t_lines * (π/2 - asin(√(1/r_m)))
#         r_lines = @. r_m*sin(θ_lines)^2
#         θ_points = @. asin(√(1/r_m)) + t_points * (π/2 - asin(√(1/r_m)))
#         r_points = @. r_m*sin(θ_points)^2
#         n_e_lines = @. 10^model.lgne_spl2d.(r_m, ts_lines)
#         n_h_lines = @. 10^model.lgnh_spl2d.(r_m, ts_lines)
#         T_e_lines = @. model.Te_spl2d.(r_m, ts_lines)

#         n_e_points = @. 10^model.lgne_spl2d.(r_m, ts_points)
#         n_h_points = @. 10^model.lgnh_spl2d.(r_m, ts_points)
#         T_e_points = @. model.Te_spl2d.(r_m, ts_points)

#         plot!(T_e_plot, r_lines, T_e_lines, lc = :black, legend = false, ylabel = L"T_e,\ \mathrm{K}", xlabel = L"r,\ \mathrm{R}_\star", margin = 5mm)
#         scatter!(T_e_plot, r_points, T_e_points, mc = :black, ms = 2, legend = false)
#         plot!(n_h_plot, r_lines, (@.  log10(n_h_lines) ), lc = :black, legend = false, ylabel = L"\lg\ n_H,\ \mathrm{cm}^{-3}", xlabel = L"r,\ \mathrm{R}_\star", margin = 5mm)
#         scatter!(n_h_plot, r_points, (@.  log10(n_h_points)), mc = :black, ms = 2, legend = false)
#         plot!(f_plot, r_lines, (@.  n_e_lines/n_h_lines), lc = :black, legend = false, ylabel = L"n_e/n_H", xlabel = L"r,\ \mathrm{R}_\star")
#         scatter!(f_plot, r_points, (@.  n_e_points/n_h_points), mc = :black, ms = 2, legend = false)
#     end
#     plt = plot(T_e_plot, n_h_plot, layout = (@layout grid(1,2)), size = (1000, 500), dpi = 300)
#     savefig(plt, "Te_nH.png")
# end