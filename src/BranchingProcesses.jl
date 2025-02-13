module BranchingProcesses

using Distributions
using AbstractTrees
using ColorSchemes
using Gadfly
using LaTeXStrings
using Printf
using SimpleNonlinearSolve
using MatrixEquations
using LinearAlgebra
using MLBase
using StatsBase
#using KrylovKit
using Arpack
using Zygote
using InvertedIndices

export SimTree, binarysplit, gathertipdata
export oupsim, branchingoupsim, varsteady, varclone, varratiofun, stdclone, oupinference
export mvoupsim, mvbranchingoupsim, mvcovarsteady, mvcovarclone, mvdeltaclone, isreversible
export driftlsqfun, driftlsqfun_grad, learndriftmodesreversible, invertdriftmodesreversible, collectdriftmodesreversible, diffusiondriftmodesreversible, learndriftmatrixreversible, learndriftreversible2
export mvoupsimclamp
export plotbranchingprocess
export normalizecols, normalizecols!, scalecols!, geomsum, invgeomsum, stdconfint, varconfint, cleanroc


include("single_variable_oup.jl")

include("multi_variable_oup.jl")

include("multi_variable_oup_perturb.jl")

include("multi_variable_oup_learning.jl")

include("tree_tools.jl")

include("plotting_tools.jl")

include("utils.jl")

end # module BranchingProcesses
