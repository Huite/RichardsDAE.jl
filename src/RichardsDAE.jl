module RichardsDAE

using Revise
using LinearAlgebra
using Statistics
using SparseArrays
using DataFrames

include("types.jl")
include("forcing.jl")
include("parameters.jl")
include("bdf.jl")
include("state.jl")
include("equations/waterbalance.jl")
include("equations/headbasedbdf1.jl")
include("equations/reducedbdf1.jl")
include("equations/daebdf1.jl")

include("solver/linear.jl")
include("solver/timestep.jl")
include("solver/relax.jl")
include("solver/newton.jl")


include("model.jl")

include("constitutive/haverkamp.jl")

include("case_helpers.jl")

end
