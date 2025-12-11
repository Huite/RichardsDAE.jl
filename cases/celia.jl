import RichardsDAE as RD
using Plots

function create_celia(formulation)
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

function run(formulations)
    models = []
    for formulation in formulations
        case = create_celia(formulation)
        # Solver is formulation dependent: DAE Jacobian is larger.
        solver = RD.NewtonSolver(
            RD.LinearSolverLU(case.parameters),
            relax = RD.ScalarRelaxation(0.0),
            maxiter = 100,
            abstol = 1e-8,
            reltol = 1e-8,
        )
        timestepper = RD.FixedTimeStepper(120.0)
        model =
            RD.Model(case.parameters, case.ψ0, solver, case.tspan, case.saveat, timestepper)
        RD.run!(model)
        push!(models, model)
    end
    return models
end

models = run((
    RD.HeadBasedBDF1(),
    RD.DAEMixedBDF1(),
    RD.ReducedBDF1(),
    #    RD.ReducedBDF2(),
))

plot(models[1].saved[:, end])
plot!(models[2].saved[:, end])
plot!(models[3].saved[:, end])

df1 = RD.waterbalance_dataframe(models[1])
df2 = RD.waterbalance_dataframe(models[2])
df3 = RD.waterbalance_dataframe(models[3])

mb1 = RD.massbalance_bias(df1)
mb2 = RD.massbalance_bias(df2)
mb3 = RD.massbalance_bias(df3)

mbrmse1 = RD.massbalance_rmse(df1)
mbrmse2 = RD.massbalance_rmse(df2)
mbrmse3 = RD.massbalance_rmse(df3)
