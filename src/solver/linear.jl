abstract type LinearSolver end

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
