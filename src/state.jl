function compute_savedflows!(state, parameters::RichardsParameters, Δt)
    ψ = primary(state)
    state.flows[1] += Δt * bottomflux(ψ, parameters, parameters.bottomboundary)
    state.flows[2] +=
        Δt * (
            topflux(ψ, parameters, parameters.topboundary) +
            forcingflux(ψ, parameters)
        )
    return
end

struct RichardsState <: State
    ψ::Vector{Float64}
    ψ_old::Vector{Float64}
    θ_old::Vector{Float64}
    ∇q::Vector{Float64}
    flows::Vector{Float64}
end

"""Return the primary state."""
function primary(state::RichardsState)
    return state.ψ
end


function prepare_state(p::RichardsParameters, initial)
    return RichardsState(
        copy(initial),  # ψ
        copy(initial),  # ψ_old,
        zero(initial),
        zero(initial),
        zeros(2),  # qbottom, qtop
    )
end

function apply_update!(state::RichardsState, linearsolver, a)
    @. state.ψ += a * linearsolver.ϕ
    return
end

function copy_state!(state::RichardsState, parameters::RichardsParameters)
    copyto!(state.ψ_old, state.ψ)
    state.θ_old .= moisture_content.(state.ψ_old, parameters.constitutive)
    return
end

function rewind!(state::RichardsState)
    copyto!(state.ψ, state.ψ_old)
end

# DAE form

struct RichardsDAEState <: State
    u::Vector{Float64}
    u_old::Vector{Float64}
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
    return RichardsDAEState(
        [copy(initial); θ],  # [ψ, θ]
        [copy(initial); copy(θ)],  # [ψ_old, θ_old]
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
    copyto!(state.u, state.u_old)
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
