module FASolverAvalanche

using Parameters # For creating structs with default values

# For convinience. These types will be used throughout.
global const FLOAT_TYPE = Ref{DataType}(Float64)
global const INT_TYPE = Ref{DataType}(Int64)
global const g = [0.0,0.0,-9.81]

# Dependencies 
import LinearAlgebra.mul! 
import LinearAlgebra.I as identity
import LinearAlgebra: norm2   
import Interpolations: linear_interpolation
import ExtendableSparse: ExtendableSparseMatrix, ILUZeroPreconditioner, factorize!, flush!
import SparseArrays: SparseMatrixCSC
import LinearSolve as LinSolv
import ForwardDiff.Dual, ForwardDiff.Partials, ForwardDiff.value 

## UTILS SUBMODULE ## 
include("./module/utils/utils.jl")
#####################

## PARSER SUBMODULE ##
include("./module/parser/parser.jl") 
######################


include("./module/cache.jl") # Cache structure for interpolations 
include("./module/init.jl") # Initializing Library 
include("./module/geometry.jl") # For precomputing mesh geometry (doesn't require interpolations)
include("./module/mesh.jl") # Mesh Computations (area, quality, etc.)
include("./module/reordering.jl") # Reordering of the mesh 
include("./module/precomputations.jl") # Pre-processing and structuring mesh 
include("./module/interpolators.jl") # Interpolators for solving at edge problems
include("./module/quality.jl") # Non-orthogonal quality estimation
include("./module/initialConditions.jl") # Initialize the discretized system
include("./module/solver.jl") # Solution of the Savage Hutter Model 
include("./module/vizualization.jl") # Saving the file to disk
end
