"""
Use Thomas for efficiency, LU for stability.
LDLT is omitted, since the Newton Jacobian is not symmetric.
"""

abstract type LinearSolver end

# [nonlinear_solve]
"""Tridiagonal linear solver."""
struct LinearSolverThomas <: LinearSolver
    n::Int
    J::Tridiagonal{Float64,Vector{Float64}}
    rhs::Vector{Float64}
    ϕ::Vector{Float64}
    y::Vector{Float64}
    B::Vector{Float64}
    function LinearSolverThomas(n)
        J = Tridiagonal(zeros(n - 1), zeros(n), zeros(n - 1))
        new(n, J, zeros(n), zeros(n), zeros(n), zeros(n))
    end
end

function Base.show(io::IO, solver::LinearSolverThomas)
    print(io, "LinearSolverThomas(n=$(solver.n))")
end

# [nonlinear_solve]
"""Thomas algorithm."""
function linearsolve!(solver::LinearSolverThomas)
    (; n, J, rhs, ϕ, y, B) = solver

    B = J.d[1]
    ϕ[1] = rhs[1] / B

    for j = 2:n
        y[j] = J.du[j-1] / B
        B = J.d[j] - J.dl[j-1] * y[j]
        if abs(B) < 1.e-12
            # This should only happen on last element of forward pass for problems
            # with zero eigenvalue. In that case the algorithmn is still stable.
            error("Beta too small!")
        end
        ϕ[j] = (rhs[j] - J.dl[j-1] * ϕ[j-1]) / B
    end

    for j = 1:(n-1)
        k = n - j
        ϕ[k] = ϕ[k] - y[k+1] * ϕ[k+1]
    end
    return
end

struct LinearSolverLU{M,F} <: LinearSolver
    n::Int
    J::M
    F::F
    rhs::Vector{Float64}
    ϕ::Vector{Float64}
end

function Base.show(io::IO, solver::LinearSolverLU)
    print(io, "LinearSolverLU(n=$(solver.n))")
end

function linearsolve!(solver::LinearSolverLU)
    # Note: This allocates since it still allocates the workspace.
    # See: https://github.com/DynareJulia/FastLapackInterface.jl
    # Or just use: https://docs.sciml.ai/LinearSolve/stable/
    lu!(solver.F, solver.J)
    # Inplace for Tridiagonal since Julia 1.11.
    # Note: stores result in B, overwrites diagonals.
    ldiv!(solver.F, solver.rhs)
    copyto!(solver.ϕ, solver.rhs)
end

function LinearSolverLU(parameters::RichardsParameters)
    n = parameters.n
    J = Tridiagonal(zeros(n - 1), zeros(n), zeros(n - 1))
    return LinearSolverLU(n, J, lu(J; check = false), zeros(n), zeros(n))
end

function LinearSolverLU(parameters::MixedDAEParameters)
    n = parameters.n
    # Block structure:
    #      [ A    I/Δt ]
    #      [ -C     I  ]

    # 1D FD stencil for A (top-left block)
    Ai = [collect(1:(n-1)); collect(1:n); collect(2:n)]
    Aj = [collect(2:n); collect(1:n); collect(1:(n-1))]

    # Upper-right block (I/Δt)
    Iui = collect(1:n)
    Iuj = collect((n+1):2n)

    # Lower-left block (-C)
    Ci = collect((n+1):2n)
    Cj = collect(1:n)

    # Lower-right block (I)
    Ili = collect((n+1):2n)
    Ilj = collect((n+1):2n)

    row_idx = [Ai; Iui; Ci; Ili]
    col_idx = [Aj; Iuj; Cj; Ilj]

    V = ones(length(row_idx))
    J = sparse(row_idx, col_idx, V, 2n, 2n)
    return LinearSolverLU(2*n, J, lu(J; check = false), zeros(2*n), zeros(2*n))
end