function residual!(rhs, state::RichardsDAEState, parameters::MixedDAEParameters, Δt)
    n = parameters.n
    Δz = parameters.Δz
    ψ = @view state.u[1:n]
    θ = @view state.u[(n+1):2*n]
    θ_old = @view state.u_old[(n+1):2*n]

    ∇q = state.∇q
    waterbalance!(∇q, ψ, parameters)
    for i = 1:n
        # -R1: differential residual
        rhs[i] = -(∇q[i] - Δz * (θ[i] - θ_old[i]) / Δt)
        # -R2: algebraic residual
        rhs[n+i] = -(θ[i] - moisture_content(ψ[i], parameters.constitutive[i]))
    end
    return
end

function jacobian!(J, state::RichardsDAEState, parameters::MixedDAEParameters, Δt)
    n = parameters.n
    Δz = parameters.Δz
    ψ = @view state.u[1:n]
    fluxJ = state.fluxJ

    dwaterbalance!(fluxJ, ψ, parameters)
    fill!(J.nzval, 0.0)

    # top-left block: A
    # Copy over from tridiagonal matrix
    for i = 1:n
        J[i, i] = fluxJ.d[i]
        if i < n
            J[i+1, i] = fluxJ.dl[i]
            J[i, i+1] = fluxJ.du[i]
        end
    end

    # top-right block: ∂F₁/∂θ = -(Δz/Δt)I
    ΔzΔt⁻¹ = Δz / Δt
    for i = 1:n
        J[i, n+i] = -ΔzΔt⁻¹
    end

    # bottom blocks: ∂F₂/∂ψ = -C*, ∂F₂/∂θ = I
    for i = 1:n
        h = parameters.constitutive[i]
        ψi = ψ[i]
        C = specific_moisture_capacity(ψi, h)
        #Sa = aqueous_saturation(ψi, h)
        #Ss = h.Ss

        # bottom-left: -C*
        #J[n+i, i] = -(C + Sa * Ss)
        J[n+i, i] = -C
        # bottom-right: I
        J[n+i, n+i] = 1.0
    end

    return
end
