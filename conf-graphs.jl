# include("RZPsc.jl")

# include("compobs.jl")
# include("plotobs.jl")
# using StaticArrays
# using LinearAlgebra
# # using TTauUtils
# # using Dierckx
# # using Interpolations
# # using Statistics
# # using LsqFit
# # using Plots
# # using LaTeXStrings
# # using SpecialFunctions
# # using Printf

mkpath("/home/coloboquito/work/papers-in-work/rzpsc/")
conf_path = "/home/coloboquito/work/papers-in-work/rzpsc/"

function makeschemegrid(best_pars, best_names, CaNa_pars, CaNa_names, grid...)
    n_pars = length(best_pars[1,:])
    n_grid = length(grid)
    println(n_grid)

    grid_sigfig = fill(8, n_grid)
    par_indeces = [1:n_grid;]
    log_axis = fill(false, n_grid)
    pars_length = zeros(Int, n_grid)
    par_axes = fill([], n_grid)
    
    for i = 1:n_grid
        par_grid = grid[i]
        if par_grid isa Array
            par_axes[i] = par_grid
            pars_length[i] = length(par_grid)
        elseif par_grid isa Tuple
            pars_length[i] = length(par_grid[1])
            par_axes[i] = par_grid[1]
            if length(par_grid) > 1
                options = par_grid[2]
                if haskey(options, :index)
                    par_indeces[i] = options.index
                end
                if haskey(options, :sigfig)
                    grid_sigfig[i] = options.sigfig
                end
                if haskey(options, :axis)
                    log_axis[i] = (options.axis == :log)
                end
            end        
        end
    end

    gridded_pars = zeros(n_pars, pars_length...)
    grid_array = zeros(pars_length...)
    min_δ = minimum(best_pars[:, end])
    thres_δ = maximum(best_pars[:, end])
    for index in keys(grid_array)
        gridded_pars[end, index] = 1e10
        for i in 1:n_grid
            i_ind = index[i]
            par = if log_axis[i]
                10.0 .^ par_axes[i][i_ind]
            else
                par_axes[i][i_ind]
            end
            gridded_pars[i, index] = par
        end
    end

    gridded_names = fill(["", ""], (pars_length...))
    # println(grid_sigfig)
    # println(par_indeces)
    # println(pars_length)
    n_models = length(best_pars[:,1])
    for i=1:n_models
        cur_pars = best_pars[i, par_indeces]
        grid_indeces = zeros(Int, n_grid)
        on_grid = true
        for j = 1:n_grid
            par_axis = if log_axis[j]
                10.0 .^ par_axes[j]
            else
                par_axes[j]
            end
            par_value = cur_pars[j]
            exp10 = floor(log10(par_value))
            grid_index_found = findall(x -> abs((x-par_value)/10.0^exp10) < 10.0^(-grid_sigfig[j]+1), par_axis)
            index_found = !isempty(grid_index_found)
            if index_found
                grid_indeces[j] = grid_index_found[1]
            else
                on_grid = false
                break
            end
        end
        if on_grid
            # if best_pars[i, end] < gridded_pars[end, grid_indeces...]
                gridded_pars[:, grid_indeces...] = best_pars[i,:]
                gridded_names[grid_indeces...] = [best_names[i][1], best_names[i][2]]
                gridded_pars[end, grid_indeces...] = 1.0
            # end
        end
    end

    n_models = length(CaNa_pars[:,1])
    for i=1:n_models
        cur_pars = CaNa_pars[i, par_indeces]
        grid_indeces = zeros(Int, n_grid)
        on_grid = true
        for j = 1:n_grid
            par_axis = if log_axis[j]
                10.0 .^ par_axes[j]
            else
                par_axes[j]
            end
            par_value = cur_pars[j]
            exp10 = floor(log10(par_value))
            grid_index_found = findall(x -> abs((x-par_value)/10.0^exp10) < 10.0^(-grid_sigfig[j]+1), par_axis)
            index_found = !isempty(grid_index_found)
            if index_found
                grid_indeces[j] = grid_index_found[1]
            else
                on_grid = false
                break
            end
        end
        if on_grid
            # if CaNa_pars[i, end] < gridded_pars[end, grid_indeces...]
                gridded_pars[:, grid_indeces...] = CaNa_pars[i,:]
                gridded_names[grid_indeces...] = [CaNa_names[i][1], CaNa_names[i][2]]
                gridded_pars[end, grid_indeces...] = 0.0
            # end
        end
    end

    return gridded_pars, gridded_names
end

# v_obs, r_obs = readobservation("spec/RZPsc_16-11-2013_proc.dat")
# stat_pars, stat_names = readmodels(star, "spec/RZPsc_16-11-2013_proc.dat", "stat_nonlocal", prof_suffix = "phot3crude")
# # stat_pars = addfluxconstant(0.01, stat_pars, stat_names)
# vsini_stat_pars, vsini_stat_names = readmodels(star, "spec/RZPsc_16-11-2013_proc.dat", "stat_nonlocal", prof_suffix = "phot3crude-sini")

# x = 1
# gridded_stat_pars, gridded_stat_names = putongrid(stat_pars, stat_names, (lgṀs, (sigfig = 3, axis = :log)), T_maxs, r_mis, Ws, angs); ""
# gridded_vsini_stat_pars, gridded_vsini_stat_names = putongrid(vsini_stat_pars, vsini_stat_names, (lgṀs, (sigfig = 3, axis = :log)), T_maxs, r_mis, Ws, angs); ""
# x = correctgridforcorotation!(gridded_stat_pars, corotationradius(star))
# x = correctgridforcorotation!(gridded_vsini_stat_pars, corotationradius(star))
plotheat(Matrix(stat_pars), (T_maxs, (index = 2,)), (lgṀs, (index = 1, axis = :log, sigfig = 3)), 
         clims = (0, 0.1), xlabel = L"T_\mathrm{max},\ \mathrm{K}", ylabel = L"\lg\dot{M},\ \mathrm{M_\odot/yr}")
savefig(conf_path*"stat_all_MT.pdf")

plotheat(Matrix(nonstat_pars), (T_maxs, (index = 2,)), (lgṀs, (index = 1, axis = :log, sigfig = 3)), 
         clims = (0, 0.1), xlabel = L"T_\mathrm{max},\ \mathrm{K}", ylabel = L"\lg\dot{M},\ \mathrm{M_\odot/yr}")
savefig(conf_path*"nonstat_all_MT.pdf")

plotheat(Matrix(best_nonstat_pars), (T_maxs, (index = 2,)), (lgṀs, (index = 1, axis = :log, sigfig = 3)), 
clims = (0, 0.1), xlabel = L"T_\mathrm{max},\ \mathrm{K}", ylabel = L"\lg\dot{M},\ \mathrm{M_\odot/yr}")
savefig(conf_path*"nonstat_best_MT.pdf")

plotheat(Matrix(best_stat_pars), (T_maxs, (index = 2,)), (lgṀs, (index = 1, axis = :log, sigfig = 3)), 
clims = (0, 0.1), xlabel = L"T_\mathrm{max},\ \mathrm{K}", ylabel = L"\lg\dot{M},\ \mathrm{M_\odot/yr}")
savefig(conf_path*"stat_best_MT.pdf")

plotheat(Matrix(CaNa_best_nonstat_pars), (T_maxs, (index = 2,)), (lgṀs, (index = 1, axis = :log, sigfig = 3)), 
clims = (0, 0.1), xlabel = L"T_\mathrm{max},\ \mathrm{K}", ylabel = L"\lg\dot{M},\ \mathrm{M_\odot/yr}")
savefig(conf_path*"nonstat_cana_MT.pdf")

plotheat(Matrix(CaNa_best_stat_pars), (T_maxs, (index = 2,)), (lgṀs, (index = 1, axis = :log, sigfig = 3)), 
clims = (0, 0.1), xlabel = L"T_\mathrm{max},\ \mathrm{K}", ylabel = L"\lg\dot{M},\ \mathrm{M_\odot/yr}")
savefig(conf_path*"stat_cana_MT.pdf")

plotheat(Matrix(nonstat_pars), (r_mis, (index = 3,)), (angs, (index = 5,)), 
         clims = (0, 0.1), xlabel = L"R_\mathrm{in},\ \mathrm{R_{star}}", ylabel = L"i,\ \degree")
savefig(conf_path*"nonstat_all_RI.pdf")

plotheat(Matrix(best_nonstat_pars), (r_mis, (index = 3,)), (angs, (index = 5,)), 
         clims = (0, 0.1), xlabel = L"R_\mathrm{in},\ \mathrm{R_{star}}", ylabel = L"i,\ \degree")
savefig(conf_path*"nonstat_best_RI.pdf")

plotheat(Matrix(CaNa_best_nonstat_pars), (r_mis, (index = 3,)), (angs, (index = 5,)), 
         clims = (0, 0.1), xlabel = L"R_\mathrm{in},\ \mathrm{R_{star}}", ylabel = L"i,\ \degree")
savefig(conf_path*"nonstat_cana_RI.pdf")

plotheat(Matrix(stat_pars), (r_mis, (index = 3,)), (angs, (index = 5,)), 
         clims = (0, 0.1), xlabel = L"R_\mathrm{in},\ \mathrm{R_{star}}", ylabel = L"i,\ \degree")
savefig(conf_path*"stat_all_RI.pdf")

plotheat(Matrix(best_stat_pars), (r_mis, (index = 3,)), (angs, (index = 5,)), 
         clims = (0, 0.1), xlabel = L"R_\mathrm{in},\ \mathrm{R_{star}}", ylabel = L"i,\ \degree")
savefig(conf_path*"stat_best_RI.pdf")

plotheat(Matrix(CaNa_best_stat_pars), (r_mis, (index = 3,)), (angs, (index = 5,)), 
         clims = (0, 0.1), xlabel = L"R_\mathrm{in},\ \mathrm{R_{star}}", ylabel = L"i,\ \degree")
savefig(conf_path*"stat_cana_RI.pdf")


plotheat(Matrix(stat_pars), (r_mis, (index = 3,)), (angs, (index = 5,)), 
         clims = (0, 0.1), xlabel = L"R_\mathrm{in},\ \mathrm{R_{star}}", ylabel = L"i,\ \degree")
savefig(conf_path*"stat_all_RW.pdf")

plotheat(Matrix(best_stat_pars), (r_mis, (index = 3,)), (angs, (index = 5,)), 
         clims = (0, 0.1), xlabel = L"R_\mathrm{in},\ \mathrm{R_{star}}", ylabel = L"i,\ \degree")
savefig(conf_path*"stat_best_RW.pdf")

plotheat(Matrix(CaNa_best_stat_pars), (r_mis, (index = 3,)), (angs, (index = 5,)), 
         clims = (0, 0.1), xlabel = L"R_\mathrm{in},\ \mathrm{R_{star}}", ylabel = L"i,\ \degree")
savefig(conf_path*"stat_cana_RW.pdf")

plotheat(Matrix(nonstat_pars), (r_mis, (index = 3,)), (Ws, (index = 4,)), 
         clims = (0, 0.1), xlabel = L"R_\mathrm{in},\ \mathrm{R_{star}}", ylabel = L"W,\ \mathrm{R_{star}}")
savefig(conf_path*"nonstat_all_RW.pdf")

plotheat(Matrix(best_nonstat_pars), (r_mis, (index = 3,)), (Ws, (index = 4,)), 
         clims = (0, 0.1), xlabel = L"R_\mathrm{in},\ \mathrm{R_{star}}", ylabel = L"W,\ \mathrm{R_{star}}")
savefig(conf_path*"nonstat_best_RW.pdf")

plotheat(Matrix(CaNa_best_nonstat_pars), (r_mis, (index = 3,)), (Ws, (index = 4,)), 
         clims = (0, 0.1), xlabel = L"R_\mathrm{in},\ \mathrm{R_{star}}", ylabel = L"W,\ \mathrm{R_{star}}")
savefig(conf_path*"nonstat_cana_RW.pdf")


δ, id = findmin(stat_pars[:,9])
plotmodel(stat_pars, stat_names, id, xlabel = L"v,\ \mathrm{km/s}", ylabel = L"I/I_c")
savefig(conf_path*"best_stat_prof.pdf")

δ, id = findmin(CaNa_best_nonstat_pars[:,9])
plotmodel(CaNa_best_nonstat_pars, CaNa_best_nonstat_names, id, xlabel = L"v,\ \mathrm{km/s}", ylabel = L"I/I_c")
savefig(conf_path*"best_prof.pdf")

scheme_gridded_stat_pars, scheme_gridded_stat_names = makeschemegrid(best_stat_pars, best_stat_names, 
                                                        CaNa_best_stat_pars, CaNa_best_stat_names, 
                                                        (lgṀs, (sigfig = 3, axis = :log)), T_maxs, r_mis, Ws, angs); ""
scheme_gridded_nonstat_pars, scheme_gridded_nonstat_names = makeschemegrid(best_nonstat_pars, best_nonstat_names, 
                                                        CaNa_best_nonstat_pars, CaNa_best_nonstat_names, 
                                                        (lgṀs, (sigfig = 3, axis = :log)), T_maxs, r_mis, Ws, angs); ""
x = correctnonstatgrid!(scheme_gridded_nonstat_pars, scheme_gridded_nonstat_names, scheme_gridded_stat_pars, scheme_gridded_stat_names)

# x = correctgridforcorotation!(gridded_nonstat_pars, corotationradius(star))

# scheme_stat_pars, scheme_stat_names = flattengrid(scheme_stat_pars, scheme_stat_names)
# scheme_nonstat_pars, scheme_nonstat_names = flattengrid(scheme_nonstat_pars, scheme_nonstat_names)