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

include("compobs.jl")


stat_pars, stat_names = loadparameters("paper-grid_RZPsc_stat.dat", 4, 4)

r_mis = [2.0:1:10.0;]
Ws = [1:0.2:4;]
T_maxs = [7000:1000:15000;]
lgṀs = [-11:0.2:-8.4;]
angs = [35:5:60;]

gridded_stat_pars, gridded_stat_names = putongrid(stat_pars, stat_names, (lgṀs, (sigfig = 3, axis = :log)), T_maxs, r_mis, Ws, angs); ""
int = interpolate((lgṀs, T_maxs, r_mis, Ws, angs), gridded_stat_pars[end,:,:,:,:,:], Gridded(Linear()))

function splitgrid(grid)
    n_pars = length(size(grid))
    grid_n = 2^(floor(Int, log2(size(grid)[1] - 1)) + 1) + 1
    splitted_grid = fill(-1.0, Tuple(repeat([grid_n], n_pars)))
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

function newmodelalg(model, axes, splits; threshold = 1.5)
    pars_min = axes[:,1]
    pars_max = axes[:,2]
    pars_Δ = pars_max .- pars_min
    n_pars = length(pars_max)
    initial_grid = fill(-1.0, Tuple(fill(2, n_pars)))
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

    split_n = 1
    grid = initial_grid
    while split_n ≤ splits
        splitted_grid = splitgrid(grid)
        n_split = size(splitted_grid)[1]
        calculated_indeces = CartesianIndex[]
        for index in keys(grid)
            push!(calculated_indeces, splittedindex(index))
        end
        
        good_indeces = splittedindex.(findall(x -> (x < thres_δ) & (x > 0), grid))
        indeces_to_calculate = findindecestocalculate(neighbours, good_indeces, calculated_indeces, n_split, n_pars)
        # println(length(indeces_to_calculate))
        # return indeces_to_calculate
        while length(indeces_to_calculate) > 0
            for index in indeces_to_calculate
                index_arr = collect(Tuple(index))
                pars = pars_min .+ pars_Δ .* (index_arr .- 1) ./ 2^split_n
                splitted_grid[index] = model(pars...)
                n_calc += 1
            end
            append!(calculated_indeces, indeces_to_calculate)
            min_δ = minimum(abs.(splitted_grid))
            thres_δ = threshold*min_δ
            println(thres_δ)
            good_indeces = indeces_to_calculate[findall(i -> splitted_grid[i] < thres_δ, indeces_to_calculate)]
            indeces_to_calculate = findindecestocalculate(neighbours, good_indeces, calculated_indeces, n_split, n_pars)
            # println(length(indeces_to_calculate))
        end
        split_n += 1
        grid = splitted_grid
        # return calculated_indeces
    end
    for grid_index in keys(grid)
        if grid[grid_index] < 0
           grid[grid_index] = -grid[grid_index] 
        end
    end
    println(n_calc)
    return grid
end

function plotgrid(grid, dims, xs, ys; clims = (0, 0.1))
    # xs = selectdim(gridded_pars, 1, dims[1])
    hm_size = size(grid, dims[1]), size(grid, dims[2])
    n_1, n_2 = hm_size
    mindims = deleteat!([1:ndims(grid);], dims)
    println(mindims)
    println(hm_size)
    hm = reshape(minimum(grid, dims = mindims), hm_size)
    heatmap(hm, clims = clims)
end