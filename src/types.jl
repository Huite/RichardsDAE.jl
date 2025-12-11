abstract type Parameters end
abstract type State end
abstract type ConstitutiveRelationships end
abstract type TimeStepper end

abstract type RichardsFormulation end
struct HeadBasedBDF1 <: RichardsFormulation end
struct DAEMixedBDF1 <: RichardsFormulation end
struct ReducedBDF1 <: RichardsFormulation end
struct ReducedBDF2 <: RichardsFormulation end
