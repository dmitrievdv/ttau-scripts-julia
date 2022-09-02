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

