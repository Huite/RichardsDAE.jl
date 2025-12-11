import RichardsDAE as RD
using Plots

function create_celia()
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
        #formulation=RD.HeadBasedBDF1(),
        #formulation = RD.ReducedBDF1(),
        formulation = RD.DAEMixedBDF1(),
        soil = soil,
        Δz = 1.0,
        Δztotal = 40.0,
        tend = 360.0,
        dt = 1.0,
        ψ0 = RD.InitialConstant(-61.5),
        bottomboundary = RD.HeadBoundary(-61.5, soil),
        topboundary = RD.HeadBoundary(-20.5, soil),
        forcing = nothing,
    )
    return celia
end

celia = create_celia()
solver = RD.NewtonSolver(
    RD.LinearSolverLU(celia.parameters),
    relax = RD.ScalarRelaxation(0.0),
    maxiter=100,
    abstol = 1e-8,
    reltol = 1e-8,
)
timestepper = RD.FixedTimeStepper(1.0)
#timestepper = RD.AdaptiveTimeStepper(Δt0=1.0)
model = RD.Model(celia.parameters, celia.ψ0, solver, celia.tspan, celia.saveat, timestepper)
RD.run!(model)

plot!(model.saved[:, end])

J = Matrix(model.solver.linearsolver.J)
cond(J)

# DAE cond: 1415
# Reduced cond: 3.75