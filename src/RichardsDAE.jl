module RichardsDAE

using Revise
using LinearAlgebra
using Statistics
using SparseArrays
using DataFrames
using DataInterpolations

include("types.jl")
include("forcing.jl")
include("parameters.jl")
include("bdf.jl")
include("state.jl")
include("equations/waterbalance.jl")
include("equations/headbased.jl")
include("equations/reduced.jl")
include("equations/mixeddae.jl")

include("solver/linear.jl")
include("solver/timestep.jl")
include("solver/relax.jl")
include("solver/newton.jl")


include("model.jl")

include("constitutive/haverkamp.jl")
include("constitutive/mualemvangenuchten.jl")

include("case_helpers.jl")

end
