struct FixedTimeStepper <: TimeStepper
    Δt0::Float64
end

# Called in implicit model
function compute_timestep_size(timestepper::FixedTimeStepper, _, converged::Bool, _)
    if !converged
        error("Failed to converge. Consider a smaller step size.")
    end
    return timestepper.Δt0
end
