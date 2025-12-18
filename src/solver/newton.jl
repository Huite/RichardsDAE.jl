""" 
    Non-linear Newton-Raphson solver, with options for backtracking/line search.

Requires a state type with the following associated methods:

* residual!
* jacobian!
* copy_state!

"""

mutable struct NewtonSolver{LS<:LinearSolver,R<:Relaxation}
    linearsolver::LS
    relax::R
    maxiter::Int
    abstol::Float64
    reltol::Float64
    njacobian::Int
    nresidual::Int
    nlinsolve::Int
    maxresidual::Vector{Float64}
end

function NewtonSolver(
    linearsolver::LinearSolver;
    relax::Relaxation = ScalarRelaxation(0.0),
    maxiter::Int = 100,
    abstol::Float64 = 1e-6,
    reltol::Float64 = 1e-6,
)
    return NewtonSolver(linearsolver, relax, maxiter, abstol, reltol, 0, 0, 0, Float64[])
end

function Base.show(io::IO, solver::NewtonSolver)
    LS = typeof(solver.linearsolver)
    R = typeof(solver.relax)
    # Get short names for the types
    ls_name = string(LS)
    r_name = string(R)
    print(
        io,
        "NewtonSolver{$ls_name,$r_name}(maxiter=$(solver.maxiter), abstol=$(solver.abstol)), reltol=$(solver.reltol),",
    )
end

function converged(newton::NewtonSolver, state)
    residual = @view newton.linearsolver.rhs[1:length(state)]
    minr, maxr = extrema(residual)
    push!(newton.maxresidual, max(maxr, abs(minr)))
    return all(
        i -> abs(residual[i]) < newton.abstol + newton.reltol * abs(state[i]),
        eachindex(residual),
    )
end
