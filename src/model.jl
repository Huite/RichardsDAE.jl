struct Model{P<:RichardsParameters,S<:State,T<:TimeStepper,NLS}
    parameters::P  # Physical parameters
    state::S  # State and dependent variables
    solver::NLS  # non-linear solver
    tspan::Tuple{Float64,Float64}
    saveat::Vector{Float64}  # frequency
    saved::Matrix{Float64}  # output
    savedflows::Matrix{Float64}  # output
    timestepper::T
end

function Base.show(io::IO, model::Model)
    p_name = Base.typename(typeof(model.parameters)).name
    println(io, "Model:")
    println(io, "  Parameters: ", p_name)
    println(io, "  Time span: ", model.tspan)
    println(io, "  Solver: ", typeof(model.solver))
    println(io, "    Linear solver: ", typeof(model.solver).parameters[1])
    if typeof(model.solver) isa NewtonSolver
        println(io, "    Relaxation: ", typeof(model.solver).parameters[2])
    end
    println(io, "  Time stepper: ", typeof(model.timestepper))
    println(io, "  Save points: ", length(model.saveat), " points")
    if !isempty(model.saveat)
        println(io, "    Range: [", first(model.saveat), ", ", last(model.saveat), "]")
    end
    print(io, "  Output size: ", size(model.saved))
end

"""
Create saveat vector of times.

    If nothing is provided for saveat, use the forcing times.
    Constrain saveat times such that `tstart <= tsave <= tend`.
"""
function create_saveat(saveat, forcing, tspan)::Vector{Float64}
    if isnothing(saveat)
        saveat = forcing.t
    end

    # Remove times beyond the time span as defined by tspan.
    tstart, tend = tspan
    first_index = searchsortedfirst(saveat, tstart)  # Find first element >= tstart
    last_index = searchsortedlast(saveat, tend)      # Find last element <= tend
    saveat = saveat[first_index:last_index]
    if isempty(saveat) || saveat[end] != tend
        push!(saveat, tend)
    end

    return saveat
end

function run!(model::Model)
    tstart, tend = model.tspan
    Δt = model.timestepper.Δt0
    t = tstart

    # Store the initial state.
    model.saved[:, 1] .= primary(model.state)
    # Accumulated flows are zero.
    model.savedflows[:, 1] .= 0.0

    tforce = model.parameters.forcing.t[1]
    save_index = 2  # t = 0.0 is already included.
    force_index = 1

    while (t < tend) && (!isapprox(t, tend))
        # New forcing
        if isapprox(t, tforce)
            force!(model.parameters, t)
            force_index += 1
            tforce =
                (force_index <= length(model.parameters.forcing.t)) ?
                model.parameters.forcing.t[force_index] : Inf
        end

        # Limit time step to not overshoot the next critical point
        tsave = model.saveat[save_index]
        t_critical = min(tsave, tforce, tend)
        if (t + Δt) > t_critical
            Δt = t_critical - t
        end
        Δt_actual = timestep!(model, Δt)
        t += Δt_actual

        # Store output
        if isapprox(t, tsave)
            model.saved[:, save_index] .= primary(model.state)
            model.savedflows[:, save_index] .= model.state.flows
            save_index += 1
        end
    end
end

function reset_and_run!(model::Model, initial)
    # Wipe results
    model.saved .= 0.0
    model.state.flows .= 0.0
    # Set initial state
    primary_state = primary(model.state)
    copyto!(primary_state, initial)
    run!(model)
    return
end

function Model(
    parameters::RichardsParameters,
    initial::Vector{Float64},
    solver,
    tspan,
    saveat,
    timestepper::TimeStepper,
)
    state = prepare_state(parameters, initial)
    saveat = create_saveat(saveat, parameters.forcing, tspan)
    nstate = length(primary(state))
    nsave = length(saveat)
    saved = zeros(nstate, nsave)
    savedflows = zeros(2, nsave)
    return Model(parameters, state, solver, tspan, saveat, saved, savedflows, timestepper)
end

"""
First order implicit (Euler Backward) time integration, with optional:

* Adaptive time stepping
* Line searches or backtracking
"""
function timestep!(model::Model, Δt)
    copy_state!(model.state, model.parameters)
    converged, n_iter = nonlinearsolve!(model.solver, model.state, model.parameters, Δt)

    while !converged
        Δt = compute_timestep_size(model.timestepper, Δt, converged, n_iter)
        rewind!(model.state)
        converged, n_iter = nonlinearsolve!(model.solver, model.state, model.parameters, Δt)
    end

    # Compute the flows based on the current solution
    compute_savedflows!(model.state, model.parameters, Δt)

    # After convergence, compute the recommended next step size based on solver performance?
    #    Δt_next = compute_next_time_step(model.timestepper, Δt, converged, n_newton_iter)
    #    return Δt, Δt_next
    return Δt
end

function nonlinearsolve!(nonlinearsolver, state, parameters, Δt)
    for i = 1:nonlinearsolver.maxiter
        residual!(nonlinearsolver.linearsolver.rhs, state, parameters, Δt)
        # Check the residual immediately for convergence.
        if converged(nonlinearsolver, primary(state))
            return true, i
        end
        jacobian!(nonlinearsolver.linearsolver.J, state, parameters, Δt)
        linearsolve!(nonlinearsolver.linearsolver)
        @show i, nonlinearsolver.linearsolver.ϕ
        relaxed_update!(
            nonlinearsolver.relax,
            nonlinearsolver.linearsolver,
            state,
            parameters,
            Δt,
        )
    end
    return false, nonlinearsolver.maxiter
end
