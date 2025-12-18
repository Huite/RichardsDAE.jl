import RichardsDAE as RD
using Plots
using DataFrames
using CSV

function create_celia(formulation, bdf)
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
        dt = 120.0,
        ψ0 = RD.InitialConstant(-61.5),
        bottomboundary = RD.HeadBoundary(-61.5, soil),
        topboundary = RD.HeadBoundary(-20.5, soil),
        forcing = nothing,
    )
    return celia
end

function create_celia_fixedq(formulation, bdf)
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
        dt = 120.0,
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
        case = create_celia(formulation, bdf)
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

models = run((
    (RD.HeadBased(), RD.BDF1(), RD.FixedTimeStepper(120.0)),
    (RD.MixedDAE(), RD.BDF1(), RD.FixedTimeStepper(120.0)),
    (RD.ReducedDAE(), RD.BDF1(), RD.FixedTimeStepper(120.0)),
    (RD.MixedDAE(), RD.BDF2(), RD.FixedTimeStepper(120.0)),
    (RD.ReducedDAE(), RD.BDF2(), RD.FixedTimeStepper(120.0)),
    (RD.ReducedDAE(), RD.BDF1(), RD.FixedTimeStepper(1.0)),
))


const COLORS = RD.okabe_ito_colors()

p = plot(xlabel = "Elevation (cm)", ylabel = "Pressure head (cm)", dpi = 300)
plot!(p, models[6].saved[:, end], label = "Dense", color = COLORS[:black], lw = 2)
plot!(p, models[1].saved[:, end], label = "ψ-based", color = COLORS[:dark_orange], lw = 2)
#plot!(models[2].saved[:, end], label="DAE-BDF1")
plot!(p, models[3].saved[:, end], label = "Celia-BDF1", color = COLORS[:light_blue], lw = 2)
#plot!(models[4].saved[:, end], label="DAE-BDF2")
plot!(p, models[5].saved[:, end], label = "Celia-BDF2", color = COLORS[:green], lw = 2)
savefig(p, "cases/output/Celia.pdf")
savefig(p, "cases/output/Celia.png")

labels = [
    "HeadBased-BDF1",
    "MixedDAE-BDF1",
    "ReducedDAE-BDF1",
    "MixedDAE-BDF2",
    "ReducedDAE-BDF2",
]
waterbalances = [RD.waterbalance_dataframe(m) for m in models[1:5]];

results = DataFrame(
    model = labels,
    mass_bias = RD.massbalance_bias.(waterbalances),
    mass_rmse = RD.massbalance_rmse.(waterbalances),
)
CSV.write("cases/output/celia.csv", results)


models[2].solver.njacobian
models[3].solver.njacobian
models[2].solver.nresidual
models[3].solver.nresidual
models[2].solver.nlinsolve
models[3].solver.nlinsolve


models[4].solver.njacobian
models[5].solver.njacobian
models[4].solver.nresidual
models[5].solver.nresidual
models[4].solver.nlinsolve
models[5].solver.nlinsolve


p = plot(
    xlabel = "Residual evaluation #",
    ylabel = "Maximum residual",
    yscale = :log10,
    dpi = 300,
)
plot!(p, models[2].solver.maxresidual, label = "DAE-BDF1", lw = 3, color = COLORS[:black])
plot!(
    p,
    models[3].solver.maxresidual,
    label = "Celia-BDF1",
    lw = 2,
    ls = :dash,
    color = COLORS[:light_blue],
)
savefig("BDF1-residual.png")
p = plot(
    xlabel = "Residual evaluation #",
    ylabel = "Maximum residual",
    yscale = :log10,
    dpi = 300,
)
plot!(p, models[4].solver.maxresidual, label = "DAE-BDF2", lw = 3, color = COLORS[:black])
plot!(
    p,
    models[5].solver.maxresidual,
    label = "Celia-BDF2",
    lw = 2,
    ls = :dash,
    color = COLORS[:green],
)
savefig("BDF2-residual.png")

p = plot(
    xlabel = "Residual evaluation #",
    ylabel = "Maximum residual",
    yscale = :log10,
    dpi = 300,
)
plot!(
    p,
    models[2].solver.maxresidual,
    label = "DAE-BDF1",
    lw = 2,
    color = COLORS[:light_blue],
)
plot!(
    p,
    models[4].solver.maxresidual,
    label = "DAE-BDF2",
    lw = 2,
    color = COLORS[:green],
    ls = :dash,
)
savefig(p, "BDF1-vs-BDF2.png")


labels = [
    "HeadBased-BDF1",
    "MixedDAE-BDF1",
    "ReducedDAE-BDF1",
    "MixedDAE-BDF2",
    "ReducedDAE-BDF2",
    "Dense-BDF1",
]
solver_stats = DataFrame(
    model = labels,
    njacobians = [m.solver.njacobian for m in models[1:6]],
    nresiduals = [m.solver.nresidual for m in models[1:6]],
    nlinsolves = [m.solver.nlinsolve for m in models[1:6]],
)
CSV.write("cases/output/solverstats.csv", solver_stats)
