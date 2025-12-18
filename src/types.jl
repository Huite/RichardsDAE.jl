abstract type Parameters end
abstract type State end
abstract type ConstitutiveRelationships end
abstract type TimeStepper end

abstract type RichardsFormulation end
struct HeadBased <: RichardsFormulation end
struct MixedDAE <: RichardsFormulation end
struct ReducedDAE <: RichardsFormulation end

abstract type BDF end
struct BDF1 <: BDF end
struct BDF2 <: BDF end
struct BDF3 <: BDF end
struct BDF4 <: BDF end
struct BDF5 <: BDF end

max_order(_::BDF1) = 1
max_order(_::BDF2) = 2
max_order(_::BDF3) = 3
max_order(_::BDF4) = 4
max_order(_::BDF5) = 5

function aqueous_saturation(ψ, C::ConstitutiveRelationships)
    return moisture_content(ψ, C) / C.θs
end
