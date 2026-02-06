function residual!(rhs, state::RichardsState, parameters::ReducedDAEParameters, Δt)
    waterbalance!(state.∇q, state.u, parameters)
    Δz = parameters.Δz
    bdf = state.bdf
    for i = 1:parameters.n
        θ = moisture_content(state.u[i], parameters.constitutive[i])
        θ_bdf = bdf.a[1] * θ
        for j = 1:bdf.order
            ψ_prev = bdf.u_prev[j, i]
            θ_prev = moisture_content(ψ_prev, parameters.constitutive[i])
            θ_bdf += bdf.a[j+1] * θ_prev
        end
        rhs[i] = -(state.∇q[i] - Δz * θ_bdf)
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
    dwaterbalance!(J, state.u, parameters)
    Δz = parameters.Δz
    a1 = state.bdf.a[1]
    for i = 1:parameters.n
        C = specific_moisture_capacity(state.u[i], parameters.constitutive[i])
        J.d[i] -= (Δz * C) * a1
    end
    return
end
