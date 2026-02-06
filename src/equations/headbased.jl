function residual!(rhs, state::RichardsState, parameters::HeadBasedParameters, Δt)
    waterbalance!(state.∇q, state.u, parameters)
    Δz = parameters.Δz
    bdf = state.bdf
    for i = 1:parameters.n
        ψ = state.u[i]
        C = specific_moisture_capacity(ψ, parameters.constitutive[i])

        ψ_bdf = bdf.a[1] * ψ
        for j = 1:bdf.order
            ψ_bdf += bdf.a[j+1] * bdf.u_prev[j, i]
        end
        rhs[i] = -(state.∇q[i] - Δz * C * ψ_bdf)
    end
    return
end

"""
This is actually a Picard formulation, since Newton is too unstable for the head-based formulation.
"""
function jacobian!(J, state, parameters::HeadBasedParameters, Δt)
    (; constitutive, Δz, bottomboundary, topboundary, n) = parameters
    ψ = state.u

    Cᵢ = J.d
    Cᵢ₊₁ = J.dl
    Cᵢ₋₁ = J.du
    Δz⁻¹ = 1.0 / Δz

    a1 = state.bdf.a[1]

    for i = 1:n
        C = specific_moisture_capacity(ψ[i], constitutive[i])
        Cᵢ[i] = -(Δz * C) * a1
    end

    k_lower = conductivity(ψ[1], constitutive[1])
    for i = 1:(n-1)
        upper = i + 1
        k_upper = conductivity(ψ[upper], constitutive[upper])
        conductance = 0.5 * (k_lower + k_upper) * Δz⁻¹
        Cᵢ₊₁[i] = conductance
        Cᵢ₋₁[i] = conductance
        # Next iteration
        k_lower = k_upper
    end

    # Then add to the diagonal term
    @views Cᵢ[1:(end-1)] .-= Cᵢ₊₁
    @views Cᵢ[2:end] .-= Cᵢ₋₁

    # Boundary conditions
    cbot = bottomboundary_coefficient(ψ, parameters, bottomboundary)
    ctop = topboundary_coefficient(ψ, parameters, topboundary)
    Cᵢ[1] += cbot
    Cᵢ[end] += ctop
    return
end
