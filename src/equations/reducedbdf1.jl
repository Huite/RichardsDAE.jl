function residual!(rhs, state::RichardsState, parameters::ReducedDAEParameters, Δt)
    waterbalance!(state.∇q, state.ψ, parameters)
    Δz = parameters.Δz
    for i = 1:parameters.n
        θ = moisture_content(state.ψ[i], parameters.constitutive[i])
        rhs[i] = -(state.∇q[i] - Δz * (θ - state.θ_old[i]) / Δt)
    end
    return
end

"""
    Construct the Jacobian matrix for the Richards equation finite difference system.
    Sets coefficients for the tridiagonal matrix representing ∂F/∂ψ from the perspective 
    of cell i, with connections to cells i-1 and i+1.

    Use Δt = ∞ for steady-state simulations.
"""
function jacobian!(J, state, parameters::ReducedDAEParameters, Δt)
    dwaterbalance!(J, state.ψ, parameters)
    Δz = parameters.Δz
    for i = 1:parameters.n
        C = specific_moisture_capacity(state.ψ[i], parameters.constitutive[i])
        #Sa = aqueous_saturation(state.ψ[i], parameters.constitutive[i])
        #Ss = parameters.constitutive[i].Ss
        #J.d[i] -= (Δz * (C + Sa * Ss)) / Δt
        J.d[i] -= (Δz * C) / Δt
    end
    return
end
