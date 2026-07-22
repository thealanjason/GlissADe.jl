module GlissADe

using Parameters # For creating structs with default values

# Some Useful constants
const FLOAT_TYPE = Ref{DataType}(Float64)
const INT_TYPE = Ref{DataType}(Int64)
const THREADS = Ref{Bool}(true)
const STATS = Ref{Bool}(true)
const PLOTS = Ref{Bool}(false)

# const g = [0.0, 0.0, -9.81]

# Dependencies
# import LinearAlgebra.mul!
# import LinearAlgebra.I as identity
# import LinearAlgebra: norm2
# import Interpolations: linear_interpolation
import ExtendableSparse: ExtendableSparseMatrix, ILUZeroPreconditioner, factorize!, flush!
import SparseArrays: SparseMatrixCSC
# import LinearSolve as LinSolv
# import ForwardDiff.Dual, ForwardDiff.Partials, ForwardDiff.value
import JLD2: load_object, save_object

## UTILS SUBMODULE ##
include("./module/utils/utils.jl")
#####################

## PARSER SUBMODULE ##
include("./module/parser/parser.jl")
######################

## GEOMETRY SUBMODULE ##
include("./module/geometry/geometry.jl")
######################

## CACHE SUBMODULE ##
include("./module/cache/cache.jl")
######################

## INIT SUBMODULE ##
include("./module/init/init.jl")
######################

## MESH SUBMODULE ##
include("./module/mesh/mesh.jl")
######################

## REORDERING SUBMODULE ##
include("./module/reordering/reordering.jl")
######################

## PRECOMPUTATIONS SUBMODULE ##
include("./module/precomputations/precomputations.jl")
######################

## INTERPOLATORS SUBMODULE ##
include("./module/interpolators/interpolators.jl")
######################

## QUALITY SUBMODULE ##
include("./module/quality/quality.jl")
######################

## INITIALCONDITIONS SUBMODULE ##
include("./module/initialConditions/initialConditions.jl")
######################

## SOLVER SUBMODULE ##
include("./module/solver/solver.jl") # Solution of the Savage Hutter Model
######################

## VISUALIZATION SUBMODULE ##
include("./module/visualization/vizualization.jl") # Saving the file to disk
######################

end
