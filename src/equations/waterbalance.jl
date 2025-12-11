function bottomflux(ψ, parameters::RichardsParameters, boundary::Nothing)
    return 0.0
end

function bottomboundary_jacobian(ψ, parameters, boundary::Nothing)
    return 0.0
end

function topflux(ψ, parameters::RichardsParameters, boundary::Nothing)
    return 0.0
end

function topboundary_jacobian(ψ, parameters, boundary::Nothing)
    return 0.0
end

function bottomboundary_coefficient(ψ, parameters::RichardsParameters, boundary::Nothing)
    return 0.0
end

function topboundary_coefficient(ψ, parameters::RichardsParameters, boundary::Nothing)
    return 0.0
end

# Precipitation

function forcingflux(ψ, parameters::RichardsParameters)
    return parameters.currentforcing[1]
end

function forcing_jacobian(ψ, parameters::RichardsParameters)
    return 0.0
end

# Store k value since it never changes
struct HeadBoundary
    ψ::Float64
    k::Float64
end

function HeadBoundary(ψ, constitutive::ConstitutiveRelationships)
    return HeadBoundary(ψ, conductivity(ψ, constitutive))
end

function bottomflux(ψ, parameters::RichardsParameters, boundary::HeadBoundary)
    kmean = 0.5 * (conductivity(ψ[1], parameters.constitutive[1]) + boundary.k)
    Δψ = boundary.ψ - ψ[1]
    Δz = 0.5 * parameters.Δz
    return kmean * (Δψ / Δz - 1)
end

function bottomboundary_jacobian(ψ, parameters::RichardsParameters, boundary::HeadBoundary)
    kmean = 0.5 * (conductivity(ψ[1], parameters.constitutive[1]) + boundary.k)
    Δψ = boundary.ψ - ψ[1]
    dk = 0.5 * dconductivity(ψ[1], parameters.constitutive[1])
    Δz = 0.5 * parameters.Δz
    return -(kmean / Δz) + dk * (Δψ / Δz - 1)
end

function bottomboundary_coefficient(
    ψ,
    parameters::RichardsParameters,
    boundary::HeadBoundary,
)
    kmean = 0.5 * (conductivity(ψ[1], parameters.constitutive[1]) + boundary.k)
    Δz = 0.5 * parameters.Δz
    return -kmean / Δz
end

function topflux(ψ, parameters::RichardsParameters, boundary::HeadBoundary)
    kmean = 0.5 * (conductivity(ψ[end], parameters.constitutive[end]) + boundary.k)
    Δψ = boundary.ψ - ψ[end]
    Δz = 0.5 * parameters.Δz
    return kmean * (Δψ / Δz + 1)
end

function topboundary_jacobian(ψ, parameters::RichardsParameters, boundary::HeadBoundary)
    kmean = 0.5 * (conductivity(ψ[end], parameters.constitutive[end]) + boundary.k)
    Δψ = boundary.ψ - ψ[end]
    dk = 0.5 * dconductivity(ψ[end], parameters.constitutive[end])
    Δz = 0.5 * parameters.Δz
    return -(kmean / Δz) + dk * (Δψ / Δz + 1)
end

function topboundary_coefficient(ψ, parameters::RichardsParameters, boundary::HeadBoundary)
    kmean = 0.5 * (conductivity(ψ[end], parameters.constitutive[end]) + boundary.k)
    Δz = 0.5 * parameters.Δz
    return -kmean / Δz
end

# Free drainage

struct FreeDrainage end

function bottomflux(ψ, parameters::RichardsParameters, boundary::FreeDrainage)
    return -conductivity(ψ[1], parameters.constitutive[1])
end

function bottomboundary_jacobian(ψ, parameters::RichardsParameters, boundary::FreeDrainage)
    return -dconductivity(ψ[1], parameters.constitutive[1])
end

function bottomboundary_coefficient(
    ψ,
    parameters::RichardsParameters,
    boundary::FreeDrainage,
)
    return 0.0
end

# Full column

function waterbalance!(∇q, ψ, parameters::RichardsParameters)
    (; constitutive, Δz, bottomboundary, topboundary, n) = parameters
    @. ∇q = 0.0
    Δz⁻¹ = 1.0 / Δz

    # Internodal flows
    k_lower = conductivity(ψ[1], constitutive[1])
    for i = 1:(n-1)
        upper = i + 1
        k_upper = conductivity(ψ[upper], constitutive[upper])
        k_inter = 0.5 * (k_lower + k_upper)
        Δψ = ψ[upper] - ψ[i]
        q = k_inter * (Δψ * Δz⁻¹ + 1)
        ∇q[i] += q
        ∇q[upper] -= q
        k_lower = k_upper
    end

    # Boundary conditions
    qbot = bottomflux(ψ, parameters, bottomboundary)
    qtop = topflux(ψ, parameters, topboundary) + forcingflux(ψ, parameters)
    ∇q[1] += qbot
    ∇q[end] += qtop
    return qbot, qtop
end

function dwaterbalance!(J, ψ, parameters::RichardsParameters)
    (; constitutive, Δz, bottomboundary, topboundary, n) = parameters

    dFᵢdψᵢ = J.d  # derivatives of F₁, ... Fₙ with respect to ψ₁, ... ψₙ
    dFᵢ₊₁dψᵢ = J.dl  # derivatives of F₂, ... Fₙ with respect to ψ₁, ... ψₙ₋₁
    dFᵢ₋₁dψᵢ = J.du  # derivatives of F₁, ... Fₙ₋₁ with respect to ψ₂, ... ψₙ
    Δz⁻¹ = 1.0 / Δz

    # First compute the off-diagonal terms -- relating to the internodal flows.
    k_lower = conductivity(ψ[1], constitutive[1])
    dk_lower = dconductivity(ψ[1], constitutive[1])
    for i = 1:(n-1)
        upper = i + 1
        k_upper = conductivity(ψ[upper], constitutive[upper])
        dk_upper = dconductivity(ψ[upper], constitutive[upper])
        conductance = 0.5 * (k_lower + k_upper) * Δz⁻¹
        Δψ = ψ[upper] - ψ[i]
        dFᵢ₊₁dψᵢ[i] = conductance - dk_lower * (Δψ * Δz⁻¹ + 1) * 0.5
        dFᵢ₋₁dψᵢ[i] = conductance + dk_upper * (Δψ * Δz⁻¹ + 1) * 0.5
        k_lower = k_upper
        dk_lower = dk_upper
    end

    dFᵢdψᵢ .= 0.0
    @views dFᵢdψᵢ[1:(end-1)] .-= dFᵢ₊₁dψᵢ
    @views dFᵢdψᵢ[2:end] .-= dFᵢ₋₁dψᵢ

    J.d[1] += bottomboundary_jacobian(ψ, parameters, bottomboundary)
    J.d[end] += topboundary_jacobian(ψ, parameters, topboundary)
    J.d[end] += forcing_jacobian(ψ, parameters)
    return
end
