using StaticArrays
using LinearAlgebra
# using TTauUtils
# using Dierckx
# using Interpolations
# using Statistics
# using LsqFit
# using Plots
# using LaTeXStrings
# using SpecialFunctions
# using Printf
using Plots
# using Distributed

n_procs = 4
# addprocs(n_procs)
# @everywhere using TTauUtils
using TTauUtils

# include("compobs.jl")


# stat_pars, stat_names = loadparameters("paper-grid_RZPsc-corr_stat.dat", 4, 4)

r_mis = [2.0:1:10.0;]
Ws = [1:0.2:4;]
T_maxs = [7000:1000.0:15000;]
lgṀs = [-11:0.2:-8.4;]
angs = [35:5.0:60;]

# gridded_stat_pars, gridded_stat_names = putongrid(stat_pars, stat_names, (lgṀs, (sigfig = 3, axis = :log)), T_maxs, r_mis, Ws, angs, empty = 1.0); ""
# int = interpolate((lgṀs, T_maxs, r_mis, Ws, angs), gridded_stat_pars[end,:,:,:,:,:], Gridded(Linear()))

function splitgrid(grid; empty = 1.0)
    n_pars = length(size(grid))
    grid_n = 2^(floor(Int, log2(size(grid)[1] - 1)) + 1) + 1
    splitted_grid = fill(empty, Tuple(repeat([grid_n], n_pars)))
    for index in keys(grid)
        index_arr = collect(Tuple(index))
        new_index_arr = 2*index_arr .- 1
        new_index = CartesianIndex(new_index_arr...)
        splitted_grid[new_index] = grid[index]
    end
    return splitted_grid
end

splittedindex(index) = 2*index - CartesianIndex(fill(1, length(index))...)

function findindecestocalculate(neighbours, good_indeces, calculated_indeces, n_split, n_pars)
    indeces_to_calculate = CartesianIndex[]
    for good_index in good_indeces
        for neighbour in neighbours
            index = good_index + neighbour
            if index in calculated_indeces 
                continue
            end
            if index in indeces_to_calculate
                continue
            end
            outside_grid = false
            for i = 1:n_pars
                if (index[i] < 1) | (index[i] > n_split)
                    outside_grid = true
                    break
                end
            end
            if outside_grid
                continue
            end
            push!(indeces_to_calculate, index)
        end
    end
    return indeces_to_calculate
end

function findtocalculate!(good_index, neighbours, calculated)
    n_splits = size(calculated)
    n_pars = length(n_splits)
    indeces_to_calculate = CartesianIndex[]
    for neighbour in neighbours
        index = good_index + neighbour
        outside_grid = false
        for i_par = 1:n_pars
            if (index[i_par] < 1) | (index[i_par] > n_splits[i_par])
                outside_grid = true
                break
            end
        end
        if outside_grid
            continue
        end
        if !calculated[index]
            calculated[index] = true
            push!(indeces_to_calculate, index)
        end
    end
    return indeces_to_calculate
end

function newmodelalgint(model, axes, splits; threshold = 1.5, dims = (1,2), empty = 1.0)
    pars_min = axes[:,1]
    pars_max = axes[:,2]
    pars_Δ = pars_max .- pars_min
    n_pars = length(pars_max)
    initial_grid = fill(empty, Tuple(fill(2, n_pars)))
    n_calc = 0
    for index in keys(initial_grid)
        index_arr = collect(Tuple(index))
        pars = pars_min .+ pars_Δ .* (index_arr .- 1)
        initial_grid[index] = model(pars...)
        n_calc += 1
    end
    # return initial_grid
    min_δ = minimum(abs.(initial_grid))
    thres_δ = threshold*min_δ

    neighbours = CartesianIndex[]
    index_arr = zeros(Int, n_pars)
    for i = 1:n_pars
        index_arr[i] = 1
        index = CartesianIndex(index_arr...)
        push!(neighbours, index)
        push!(neighbours, -index)
        index_arr[i] = 0
    end
    for dir = 0:2^n_pars-1
        index_arr = digits(dir, base = 2, pad = n_pars)*2 .- 1
        push!(neighbours, CartesianIndex(index_arr...))
    end
    plots = []
    split_n = 1
    grid = initial_grid
    push!(plots, plotgrid(grid, dims))
    while split_n ≤ splits
        splitted_grid = splitgrid(grid)
        calculated = fill(false, size(splitted_grid))
        # to_calculate = fill(false, size(splitted_grid))
        indeces_to_calculate = CartesianIndex[]
        for index in keys(grid)
            δ_grid = grid[index]
            splitted_index = splittedindex(index)
            splitted_grid[splitted_index] = grid[index]
            calculated[splitted_index] = true
            if δ_grid < thres_δ
                # println(splitted_index)
                # print(sum(calculated), " ")
                found = findtocalculate!(splitted_index, neighbours, calculated)
                # print(length(found), "(", length(neighbours), ") ")
                # println(sum(calculated))
                append!(indeces_to_calculate, found)
            end
        end
        push!(plots, plotgrid(splitted_grid, dims))
        while (length(indeces_to_calculate) > 0)
            for index in indeces_to_calculate
                index_arr = collect(Tuple(index))
                pars = pars_min .+ pars_Δ .* (index_arr .- 1) ./ 2^split_n
                splitted_grid[index] = model(pars...)
                n_calc += 1
                # calculated[index] = true
            end
            min_δ = minimum(abs.(splitted_grid))
            thres_δ = threshold*min_δ
            calculated_indeces = deepcopy(indeces_to_calculate)
            indeces_to_calculate = CartesianIndex[]
            for index in calculated_indeces
                if splitted_grid[index] < thres_δ
                    # println(index)
                    # print(sum(calculated), " ")
                    found = findtocalculate!(index, neighbours, calculated)
                    # print(length(found), "(", length(neighbours), ") ")
                    # println(sum(calculated))
                    append!(indeces_to_calculate, found)
                end
            end
            # println(thres_δ)
            # println(length(indeces_to_calculate))
            push!(plots, plotgrid(splitted_grid, dims))
        end
        split_n += 1
        grid = splitted_grid
        print(n_calc, " ")
        println(sum(calculated))
        # return calculated_indeces
    end
    for grid_index in keys(grid)
        if grid[grid_index] < 0
           grid[grid_index] = -grid[grid_index] 
        end
    end
    # println(n_calc)
    return grid, plots
end

function plotgrid(grid, dims; clims = (0, 0.1))
    # xs = selectdim(gridded_pars, 1, dims[1])
    hm_size = size(grid, dims[1]), size(grid, dims[2])
    n_1, n_2 = hm_size
    mindims = deleteat!([1:ndims(grid);], dims)
    # println(mindims)
    # println(hm_size)
    hm = reshape(minimum(grid, dims = mindims), hm_size)
    heatmap(hm, clims = clims)
end

function linename(u, l)
    line_name = ""
    if l == 2
        line_name *= 'H'
    end
    if l == 3
        line_name *= "Pa"
    end
    if l == 4
        line_name *= "Br"
    end
    line_name *= 'a' + u-1-l
end

function gennamesandpars(indeces_to_calculate, n_splits, gridname, par_axes, line)
    pars_min = par_axes[:,1]
    pars_max = par_axes[:,2]
    pars_Δ = pars_max .- pars_min
    n_pars = length(pars_min)
    n_ind = length(indeces_to_calculate)
    pars = zeros(n_pars, n_ind)
    names = String[]
    rational_index = zeros(Rational{Int}, n_pars)
    for i_ind in 1:n_ind
        index = indeces_to_calculate[i_ind]
        model_name = gridname
        for i_par = 1:n_pars-1
            rational_index[i_par] = (index[i_par] - 1)//(n_splits[i_par] - 1)
            pars[i_par, i_ind] = pars_min[i_par] + rational_index[i_par] * pars_Δ[i_par]
            model_name = model_name*"_$(rational_index[i_par].num)-$(rational_index[i_par].den)"
        end
        i_par = n_pars
        rational_index[i_par] = (index[i_par] - 1)//(n_splits[i_par] - 1)
        pars[i_par, i_ind] = pars_min[i_par] + rational_index[i_par] * pars_Δ[i_par]
        profile_name = linename(line...)*"_$(rational_index[i_par].num)-$(rational_index[i_par].den)"
        push!(names, "$model_name/$profile_name")
    end
    return names, pars
end

function computemodels(names, pars, star)
    model_names = String[]
    n_pars, n_profs = size(pars)
    model_name, profile_name = split(names[1], "/")
    push!(model_names, model_name)
    model_pars = pars[:, 1]
    for i_prof = 2:n_profs
        model_name, profile_name = split(names[i_prof], "/")
        if !(model_name in model_names)
            push!(model_names, model_name)
            model_pars = [model_pars pars[:, i_prof]]
        end
    end

    n_models = length(model_names)
    n_iter = n_models ÷ n_procs

    for i_iter = 1:n_iter
        iter_start = (i_iter-1)*n_procs + 1
        iter_end_outside = (i_iter*n_procs > n_models)
        iter_end = iter_end_outside*n_models + !iter_end_outside*i_iter*n_procs
        println(i_iter, " ", model_names[iter_start:iter_end])
        # distributed model loop...
        # saving mdels...
    end

    

    n_iter = n_profs ÷ n_procs

    for i_iter = 1:n_iter
        iter_start = (i_iter-1)*n_procs + 1
        iter_end_outside = (i_iter*n_procs > n_profs)
        iter_end = iter_end_outside*n_profs + !iter_end_outside*i_iter*n_procs
        iter_range = iter_start:iter_end
        iter_length = length(iter_range)
        iter_models = String[] # HydrogenModel in the future...
        profiles_iter_model_indeces = zeros(Int, iter_length)

        for i_prof_iter = 1:iter_length
            i_prof = iter_range[i_prof_iter]
            model_name, profile_name = split(names[i_prof], "/")
            if !(model_name in iter_models)
                push!(iter_models, model_name) # loadmodel in the future...
                profiles_iter_model_indeces[i_prof_iter] = length(iter_models)
            end
        end

        print("$i_iter: ")
        for i_prof_iter = 1:iter_length  # distributed in the future...
            i_prof = iter_range[i_prof_iter]
            i_model = profiles_iter_model_indeces[i_prof_iter]
            print("($i_prof, $i_model, $(names[i_prof])) ")
        end
        println("")
        # saving profiles...
    end
end

function newmodelalg(star, gridname, par_axes, split, line; threshold = 2)
    pars_min = par_axes[:,1]
    pars_max = par_axes[:,2]
    pars_Δ = pars_max .- pars_min
    n_pars = length(pars_max)
    initial_grid = fill(empty, Tuple(fill(2, n_pars)))
    n_calc = 0

    neighbours = CartesianIndex[]
    index_arr = zeros(Int, n_pars)
    for i = 1:n_pars
        index_arr[i] = 1
        index = CartesianIndex(index_arr...)
        push!(neighbours, index)
        push!(neighbours, -index)
        index_arr[i] = 0
    end
    for dir = 0:2^n_pars-1
        index_arr = digits(dir, base = 2, pad = n_pars)*2 .- 1
        push!(neighbours, CartesianIndex(index_arr...))
    end
    n_splits = size(initial_grid)
    indeces_to_calculate = keys(initial_grid)
    names, pars = gennamesandpars(indeces_to_calculate, n_splits, gridname, par_axes, line)
    computemodels(names, pars, star)

    # return names, pars

    # for index in keys(initial_grid)
    #     index_arr = collect(Tuple(index))
    #     pars = pars_min .+ pars_Δ .* (index_arr .- 1)
    #     initial_grid[index] = model(pars...)
    #     n_calc += 1
    # end
    # return initial_grid
    # min_δ = minimum(abs.(initial_grid))
    # thres_δ = threshold*min_δ
end

star = Star("RZPsc")
parameter_grid = [lgṀs, T_maxs, r_mis, Ws, angs]
par_axes = Float64[-11 -8;    # lgṀ
                  7000 15000; # T_max
                     2 10;    # R_in
                     1 4;     # W
                    30 60;    # i
                  ]
newmodelalg(star, "grid", par_axes, 2, (3,2))