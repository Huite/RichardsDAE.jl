function compute_savedflows!(state, parameters::RichardsParameters, Δt)
    ψ = primary(state)
    n = parameters.n

    qbot_new = bottomflux(ψ, parameters, parameters.bottomboundary)
    qtop_new = topflux(ψ, parameters, parameters.topboundary) + forcingflux(ψ, parameters)

    if state.bdf.order == 1
        state.flows[1] += Δt * qbot_new
        state.flows[2] += Δt * qtop_new
    else
        # Get ψ from previous timestep (first column of history)
        # Trapezoidal approximation
        ψ_prev = @view state.bdf.u_prev[1, 1:n]
        qbot_prev = bottomflux(ψ_prev, parameters, parameters.bottomboundary)
        qtop_prev =
            topflux(ψ_prev, parameters, parameters.topboundary) +
            forcingflux(ψ_prev, parameters)

        state.flows[1] += 0.5 * Δt * (qbot_new + qbot_prev)
        state.flows[2] += 0.5 * Δt * (qtop_new + qtop_prev)
    end
    return
end


struct RichardsState <: State
    u::Vector{Float64}
    bdf::BDFWorkSpace
    ∇q::Vector{Float64}
    flows::Vector{Float64}
end

"""Return the primary state."""
function primary(state::RichardsState)
    return state.u
end


function prepare_state(p::RichardsParameters, initial)
    bdf = BDFWorkSpace(p)
    bdf.u_prev[1, :] .= initial
    return RichardsState(
        copy(initial),  # ψ
        bdf,
        zero(initial),
        zeros(2),  # qbottom, qtop
    )
end

function apply_update!(state::RichardsState, linearsolver, a)
    @. state.u += a * linearsolver.ϕ
    return
end

function copy_state!(state::RichardsState, parameters::RichardsParameters)
    copyto!(state.bdf.u_prev[:, 1], state.u)
    return
end

function rewind!(state::RichardsState)
    copyto!(state.u, state.bdf.u_prev[:, 1])
    return
end

# DAE form

struct RichardsDAEState <: State
    u::Vector{Float64}
    bdf::BDFWorkSpace
    ∇q::Vector{Float64}
    flows::Vector{Float64}
    fluxJ::Tridiagonal{Float64,Vector{Float64}}
end

function primary(state::RichardsDAEState)
    u = state.u
    n = Int(length(u) // 2)
    ψ = @view u[1:n]
    return ψ
end

function prepare_state(p::MixedDAEParameters, initial)
    n = p.n
    θ = [moisture_content(ψ, p.constitutive[i]) for (i, ψ) in enumerate(initial)]
    bdf = BDFWorkSpace(p)
    bdf.u_prev[1, 1:n] .= initial
    bdf.u_prev[1, (n+1):end] .= θ
    return RichardsDAEState(
        [copy(initial); θ],  # [ψ, θ]
        bdf,
        zero(initial),  # ∇q
        zeros(2),  # qbottom, qtop
        Tridiagonal(zeros(n-1), zeros(n), zeros(n-1)),
    )
end

function apply_update!(state::RichardsDAEState, linearsolver, a)
    @. state.u += a * linearsolver.ϕ
    return
end

function copy_state!(state::RichardsDAEState, parameters)
    copyto!(state.u_old, state.u)
    return
end

function rewind!(state::RichardsDAEState)
    copyto!(state.u, state.bdf.u_prev[1, :])
end

# Initial states

@kwdef struct InitialHydrostatic
    watertable::Float64
end

@kwdef struct InitialConstant
    ψ::Float64
end

function initialψ(I::InitialHydrostatic, Δz, Δztotal, n)
    z = collect(Δz:Δz:Δztotal)  # heights from bottom
    return I.watertable .- z
end

function initialψ(I::InitialConstant, Δz, Δztotal, n)
    return fill(I.ψ, n)
end
