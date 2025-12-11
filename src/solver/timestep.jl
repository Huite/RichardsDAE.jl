struct FixedTimeStepper <: TimeStepper
    Δt0::Float64
end

# Called in implicit model
function compute_timestep_size(timestepper::FixedTimeStepper, _, converged::Bool, _)
    if !converged
        error("Failed to converge. Consider using adaptive time stepping.")
    end
    return timestepper.Δt0
end

# Called in explicit model
function compute_timestep_size(
    timestepper::FixedTimeStepper,
    state::State,
    parameters::RichardsParameters,
    Δt,
)
    return timestepper.Δt0
end

# [implicit]
struct AdaptiveTimeStepper <: TimeStepper
    Δt0::Float64
    n_increase::Int
    increase::Float64
    n_decrease::Int
    decrease::Float64
    failure::Float64
    Δtmin::Float64
    Δtmax::Float64
    function AdaptiveTimeStepper(;
        Δt0 = nothing,
        n_increase = 5,
        increase = 1.25,
        n_decrease = 15,
        decrease = 0.9,
        failure = 0.5,
        Δtmin = 1e-6,
        Δtmax = 1.0,
    )
        if isnothing(Δt0)
            Δt0 = sqrt(Δtmin * Δtmax)
        end
        return new(Δt0, n_increase, increase, n_decrease, decrease, failure, Δtmin, Δtmax)
    end
end

# [implicit]
"""
Modify time step based on convergence behavior.
"""
function compute_timestep_size(
    timestepper::AdaptiveTimeStepper,
    Δt,
    converged,
    n_newton_iter,
)
    if !converged
        newΔt = Δt * timestepper.failure
    elseif n_newton_iter > timestepper.n_decrease
        newΔt = Δt * timestepper.decrease
    elseif n_newton_iter < timestepper.n_increase
        newΔt = Δt * timestepper.increase
    else
        newΔt = Δt
    end
    if newΔt < timestepper.Δtmin
        error("Time step below Δtmin")
    end
    return min(newΔt, timestepper.Δtmax)
end
