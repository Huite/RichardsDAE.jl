import RichardsDAE as RD
using Plots
using DataFrames
using CSV

function create_celia(formulation, bdf)
    # This is the original benchmark, with Dirichlet boundaries.
    # Note: units are centimeters and seconds!
    soil = RD.Haverkamp(
        a = 1.611e6,
        B = 3.96,
        y = 4.74,
        A = 1.175e6,
        ks = 0.00944,
        θs = 0.287,
        θr = 0.075,
        Ss = 0.0,
    )
    celia = RD.Case(
        formulation = formulation,
        bdf = bdf,
        soil = soil,
        Δz = 1.0,
        Δztotal = 40.0,
        tend = 360.0,
        save_dt = 120.0,
        ψ0 = RD.InitialConstant(-61.5),
        bottomboundary = RD.HeadBoundary(-61.5, soil),
        topboundary = RD.HeadBoundary(-20.5, soil),
        forcing = nothing,
    )
    return celia
end

function create_celia_fixedq(formulation, bdf)
    # This is the adapted benchmark with fixed rates.
    # Note: units are centimeters and seconds!
    soil = RD.Haverkamp(
        a = 1.611e6,
        B = 3.96,
        y = 4.74,
        A = 1.175e6,
        ks = 0.00944,
        θs = 0.287,
        θr = 0.075,
        Ss = 0.0,
    )
    celia = RD.Case(
        formulation = formulation,
        bdf = bdf,
        soil = soil,
        Δz = 1.0,
        Δztotal = 40.0,
        tend = 360.0,
        save_dt = 120.0,
        ψ0 = RD.InitialConstant(-61.5),
        bottomboundary = nothing,
        topboundary = nothing,
        forcing = RD.MeteorologicalForcing([0.0], [0.005], [0.0]),
    )
    return celia

end

function run(formulations)
    models = []
    for (formulation, bdf, timestepper) in formulations
        case = create_celia_fixedq(formulation, bdf)
        # Solver is formulation dependent: DAE Jacobian is larger.
        solver = RD.NewtonSolver(
            RD.LinearSolverLU(case.parameters),
            relax = RD.ScalarRelaxation(0.0),
            maxiter = 100,
            abstol = 1e-8,
            reltol = 1e-8,
        )
        model =
            RD.Model(case.parameters, case.ψ0, solver, case.tspan, case.saveat, timestepper)
        RD.run!(model)
        push!(models, model)
    end
    return models
end

# Comparison results for pressure head plot.

models = run((
    (RD.HeadBased(), RD.BDF1(), RD.FixedTimeStepper(120.0)),
    (RD.ReducedDAE(), RD.BDF1(), RD.FixedTimeStepper(120.0)),
    (RD.ReducedDAE(), RD.BDF2(), RD.FixedTimeStepper(120.0)),
    (RD.ReducedDAE(), RD.BDF1(), RD.FixedTimeStepper(1.0)),
))
labels = [
    "ψ-based (Δt=120.0)",
    "Reduced-BDF1 (Δt=120.0)",
    "Reduced-BDF2 (Δt=120.0)",
    "Reduced-BDF1 (Δt=1.0)",
]
data = Dict{String,Vector{Float64}}()
for (label, model) in zip(labels, models)
    n = 40
    finalψ = model.saved[1:n, end]
    data[label] = finalψ
end

headdf = DataFrame(data)
CSV.write("cases/output/haverkamp-final-head.csv", headdf)

# Mass conservation diagnostics.

models = run((
    (RD.HeadBased(), RD.BDF1(), RD.FixedTimeStepper(120.0)),
    (RD.MixedDAE(), RD.BDF1(), RD.FixedTimeStepper(120.0)),
    (RD.ReducedDAE(), RD.BDF1(), RD.FixedTimeStepper(120.0)),
    (RD.MixedDAE(), RD.BDF2(), RD.FixedTimeStepper(120.0)),
    (RD.ReducedDAE(), RD.BDF2(), RD.FixedTimeStepper(120.0)),
))
labels = (
    "HeadBased-BDF1",
    "MixedDAE-BDF1",
    "ReducedDAE-BDF1",
    "MixedDAE-BDF2",
    "ReducedDAE-BDF2",
)
waterbalances = [RD.waterbalance_dataframe(m) for m in models[1:5]];

results = DataFrame(
    model = labels,
    mass_balance_ratio = RD.massbalance_balance_ratio.(waterbalances),
    mass_rmse = RD.massbalance_rmse.(waterbalances),
)
CSV.write("cases/output/haverkamp-mass.csv", results)

