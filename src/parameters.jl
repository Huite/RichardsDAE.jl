
struct RichardsParameters{F,C,T,B}
    formulation::F
    constitutive::Vector{C}
    Δz::Float64
    forcing::MeteorologicalForcing
    bottomboundary::B
    topboundary::T
    n::Int
    currentforcing::Vector{Float64}  # P, ET

    function RichardsParameters(
        formulation::F;
        constitutive::Vector{C},
        Δz,
        forcing,
        bottomboundary::B,
        topboundary::T,
    ) where {F<:RichardsFormulation,C,T,B}
        new{F,C,T,B}(
            formulation,
            constitutive,
            Δz,
            forcing,
            bottomboundary,
            topboundary,
            length(constitutive),
            zeros(Float64, 2),
        )
    end
end


function Base.show(io::IO, rp::RichardsParameters)
    C = eltype(rp.constitutive)
    T = typeof(rp.topboundary)
    B = typeof(rp.bottomboundary)

    # Get clean type names
    rp_name = string(Base.typename(typeof(rp)).name)
    c_name = string(Base.typename(C).name)
    t_name = string(Base.typename(T).name)
    b_name = string(Base.typename(B).name)

    println(io, "$rp_name{$c_name,$t_name,$b_name}:")
    println(io, "  Grid: $(rp.n) layers, Δz = $(rp.Δz)")
    println(io, "  Constitutive: $(c_name)")
    println(io, "  Bottom boundary: ", rp.bottomboundary)
    println(io, "  Top boundary: ", rp.topboundary)
    println(io, "  Meteorological forcing: $(length(rp.forcing.t)) time steps")
    print(
        io,
        "  Current forcing: P = $(rp.currentforcing[1]), ET = $(rp.currentforcing[2])",
    )
end


const MixedDAEParameters{C,T,B} = RichardsParameters{DAEMixedBDF1,C,T,B}
const HeadBasedParameters{C,T,B} = RichardsParameters{HeadBasedBDF1,C,T,B}
const ReducedDAEParameters{C,T,B} = RichardsParameters{ReducedBDF1,C,T,B}