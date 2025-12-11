function residual!(rhs, state::RichardsState, parameters::HeadBasedParameters, Δt)
    waterbalance!(state.∇q, state.ψ, parameters)
    Δz = parameters.Δz
    for i = 1:parameters.n
        ψ = state.ψ[i]
        C = specific_moisture_capacity(ψ, parameters.constitutive[i])
        #Sa = aqueous_saturation(ψ, parameters.constitutive[i])
        #Ss = parameters.constitutive[i].Ss
        #rhs[i] = -(state.∇q[i] - Δz * (C + Sa * Ss) * (ψ - state.ψ_old[i]) / Δt)
        rhs[i] = -(state.∇q[i] - Δz * C * (ψ - state.ψ_old[i]) / Δt)
    end
    return
end

"""
This is actually a Picard formulation, since Newton is too unstable for the head-based formulation.
"""
function jacobian!(J, state, parameters::HeadBasedParameters, Δt)
    (; constitutive, Δz, bottomboundary, topboundary, n) = parameters
    ψ = state.ψ

    Cᵢ = J.d
    Cᵢ₊₁ = J.dl
    Cᵢ₋₁ = J.du
    Δz⁻¹ = 1.0 / Δz

    for i = 1:n
        C = specific_moisture_capacity(ψ[i], constitutive[i])
        #Sa = aqueous_saturation(ψ[i], constitutive[i])
        #Ss = constitutive[i].Ss
        #Cᵢ[i] = -(Δz * (C + Sa * Ss)) / Δt
        Cᵢ[i] = -(Δz * C) / Δt
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
