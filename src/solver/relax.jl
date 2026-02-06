abstract type Relaxation end
abstract type LineSearch <: Relaxation end

struct ScalarRelaxation <: Relaxation
    relax::Float64
    function ScalarRelaxation(relax::Float64)
        0 <= relax < 1 ||
            throw(ArgumentError("Relaxation parameter must be in [0,1): got $relax"))
        new(relax)
    end
end

function relaxed_update!(
    relaxation::ScalarRelaxation,
    linearsolver,
    state,
    parameters,
    Δt,
)::Bool
    apply_update!(state, linearsolver, 1.0 - relaxation.relax)
    return true
end

@kwdef struct SimpleLineSearch <: LineSearch
    a0::Float64 = 1.0  # initial step size
    b::Float64 = 0.5  # backtracking reduction factor
    c::Float64 = 1e-4  # Sufficient decrease for Armijo condition
    minstep::Float64 = 1e-2
    maxiter::Int = 5
end

function compute_step(ls::SimpleLineSearch, _, α₂, _, _, _)
    return α₂, max(ls.b * α₂, ls.minstep * α₂)
end

# [nonlinear_solve]
function relaxed_update!(ls::LineSearch, linearsolver, state, parameters, Δt)
    # α₀ = 0.0 (implicit)
    α₁ = 0.0
    α₂ = ls.a0
    # Compute the L2 norm of the residual to check for convergence
    L2₀ = norm(linearsolver.rhs)
    L2₁ = L2₀

    L2best = L2₁
    αbest = α₂

    for _ = 1:ls.maxiter
        # Take a step
        apply_update!(state, linearsolver, α₂)
        residual!(linearsolver.rhs, state, parameters, Δt)
        L2₂ = norm(linearsolver.rhs)

        # Armijo condition for sufficient decrease
        if L2₂ <= ((1 - ls.c * α₂) * L2₀)
            return true
        end
        if L2₂ < L2best
            L2best = L2₂
            αbest = α₂
        end

        # Undo the step by applying a NEGATIVE update
        apply_update!(state, linearsolver, -α₂)

        # Compute new step size
        α₁, α₂ = compute_step(ls, α₁, α₂, L2₀, L2₁, L2₂)
        L2₁ = L2₂
    end
    # Achieved maximum iterations, use αbest.
    apply_update!(state, linearsolver, αbest)
    return false
end
