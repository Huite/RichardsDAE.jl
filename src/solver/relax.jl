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

@kwdef struct QuadraticLineSearch <: LineSearch
    a0::Float64 = 1.0  # initial step size
    c::Float64 = 1e-4  # Sufficient decrease for Armijo condition
    maxiter::Int = 5
    low::Float64 = 0.01  # interpolation bounds
    high::Float64 = 0.9  # interpolation bounds
end

function quadratic_minimum(α₂, L2₀, L2₂)
    # Fit y = f(x) = ax² + bx + c
    # Then find its minimum: -b/2a
    #
    # f'(x) = 2ax + b, hence:
    # b = f'(0)
    # c = f(0)
    # a = (f(x₁) - bx₁ - c) / x₁²
    #
    # Hence: 
    # -b / 2a = f'(0) / (2 (f(x₁) - bx₁ - c) / x₁²) = 
    # -b / 2a = (f'(0) * x₁²) / 2 (f(x₁) - bx₁ - c)
    #
    # Note that for the L2 norm specifically: b = -L2₀

    grad₀ = -L2₀
    y1 = L2₂
    x1 = α₂
    c = L2₀
    b = grad₀
    newstep = -b * (x1^2) / (2 * y1 - b * x1 - c)
    return newstep
end

function cubic_minimum(α₁, α₂, L2₀, L2₁, L2₂)
    # Fit y = f(x) = ax³ + bx² + cx + d
    # Then find its minimum by solving f'(x) = 0
    # f'(x) = 3ax² + 2bx + c = 0
    #
    # We use function values at x = 0, x = α₁, x = α₂
    # And the derivative at x = 0
    #
    # Known values:
    # f(0) = L2₀ = d
    # f'(0) = -L2₀ = c
    # f(α₁) = L2₁
    # f(α₂) = L2₂


    # Initial slope at x = 0 (specifically for L2 norm)
    c = -L2₀

    # Compute the differences between actual and linear approximation
    diff₁ = L2₁ - L2₀ - c * α₁
    diff₂ = L2₂ - L2₀ - c * α₂

    # Calculate cubic coefficients using the differences
    # This matrix solution is derived from the system of equations:
    # a(α₁)³ + b(α₁)² = diff₁
    # a(α₂)³ + b(α₂)² = diff₂
    div = 1.0 / (α₁^2 * α₂^2 * (α₂ - α₁))
    a = (α₁^2 * diff₂ - α₂^2 * diff₁) * div
    b = (-α₁^3 * diff₂ + α₂^3 * diff₁) * div

    # Find the minimizer of the cubic
    # If coefficient of x³ is effectively zero, fall back to quadratic model
    if abs(a) < eps(Float64)
        # For a quadratic model: minimum at x = -c/(2b)
        newstep = -c / (2 * b)
    else
        # Solving the quadratic: 3ax² + 2bx + c = 0
        # Using the quadratic formula: x = (-2b ± √(4b² - 4*3a*c))/(2*3a)
        # Simplified: x = (-b ± √(b² - 3a*c))/(3a)
        # We take the formula that gives us a positive step
        d = max(b^2 - 3 * a * c, 0.0)
        # Avoids catastrophic cancellation when b and √d are close
        newstep = -c / (b * sqrt(d))
    end
    return newstep
end

function compute_step(ls::QuadraticLineSearch, α₁, α₂, L2₀, L2₁, L2₂)
    newstep = quadratic_minimum(α₂, L2₀, L2₂)
    if isnan(newstep)
        newstep = α₂ * ls.high
    end
    newstep = clamp(newstep, α₂ * ls.low, α₂ * ls.high)
    return α₂, newstep
end


@kwdef struct CubicLineSearch <: LineSearch
    a0::Float64 = 1.0  # initial step size
    c::Float64 = 1e-4  # Sufficient decrease for Armijo condition
    maxiter::Int = 5
    low::Float64 = 0.01  # interpolation bounds
    high::Float64 = 0.9  # interpolation bounds
end

function compute_step(ls::CubicLineSearch, α₁, α₂, L2₀, L2₁, L2₂)
    if isapprox(L2₀, L2₁)
        # Will trigger first iteration
        newstep = quadratic_minimum(α₂, L2₀, L2₂)
    else
        newstep = cubic_minimum(α₁, α₂, L2₀, L2₁, L2₂)
    end

    # Check for NaN, and bound the step size within low and high multipliers
    # of the current step.
    if isnan(newstep)
        newstep = α₂ * ls.high
    end
    newstep = clamp(newstep, α₂ * ls.low, α₂ * ls.high)

    return α₂, newstep
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
