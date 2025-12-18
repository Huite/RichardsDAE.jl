struct RichardsParameters{F,O,C,T,B}
    formulation::F
    bdforder::O
    constitutive::Vector{C}
    Δz::Float64
    forcing::MeteorologicalForcing
    bottomboundary::B
    topboundary::T
    n::Int
    currentforcing::Vector{Float64}  # P, ET

    function RichardsParameters(
        formulation::F,
        bdforder::O;
        constitutive::Vector{C},
        Δz,
        forcing,
        bottomboundary::B,
        topboundary::T,
    ) where {F<:RichardsFormulation,O<:BDF,C,T,B}
        new{F,O,C,T,B}(
            formulation,
            bdforder,
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

const MixedDAEParameters{O,C,T,B} = RichardsParameters{MixedDAE,O,C,T,B}
const HeadBasedParameters{O,C,T,B} = RichardsParameters{HeadBased,O,C,T,B}
const ReducedDAEParameters{O,C,T,B} = RichardsParameters{ReducedDAE,O,C,T,B}

nunknown(p::RichardsParameters) = p.n
nunknown(p::MixedDAEParameters) = 2 * p.n
