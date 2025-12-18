import RichardsDAE as RD
using Plots
using DataFrames
using CSV


function create_millersand(formulation, bdf)
    sand = RD.ModifiedMualemVanGenuchten(
        a = 5.470,
        n = 4.264,
        l = 0.5,
        ks = 5.040,
        θr = 0.093,
        θs = 0.301,
        Ss = 1e-6,
        ψe = -1e-2,
    )
    sandspline = RD.SplineConstitutive(sand)
    millersand = RD.Case(
        formulation = formulation,
        bdf = bdf,
        soil = sandspline,
        Δz = 0.0125,
        Δztotal = 10.0,
        tend = 0.18,
        dt = 0.01,
        ψ0 = RD.InitialHydrostatic(watertable = 0.0),
        topboundary = RD.HeadBoundary(0.1, sandspline),
        bottomboundary = RD.HeadBoundary(0.0, sandspline),
        forcing = nothing,
    )
    return millersand
end

case = create_millersand(RD.ReducedDAE(), RD.BDF1())
solver = RD.NewtonSolver(
    RD.LinearSolverThomas(case.parameters),
    #RD.LinearSolverLU(case.parameters),
    relax = RD.SimpleLineSearch(),
    maxiter = 150,
    abstol = 1e-7,
    reltol = 1e-7,
)
timestepper = RD.FixedTimeStepper(0.001)#(Δt0=1e-2, Δtmin=1e-2, Δtmax=1.0e-2)
#timestepper = RD.AdaptiveTimeStepper(Δt0=1e-3, Δtmin=1e-6, Δtmax=1.0e-2)
model = RD.Model(case.parameters, case.ψ0, solver, case.tspan, case.saveat, timestepper)
RD.run!(model)

plot(model.saved[:, end])


J = model.solver.linearsolver.J
ψ = model.state.u
extrema(ψ)
k = [RD.conductivity(ψ[i], case.parameters.constitutive[i]) for i = 1:800]
plot(k)

extrema(ψ)

##

function create_millerloam(formulation, bdf)
    loam = RD.ModifiedMualemVanGenuchten(
        a = 3.600,
        n = 1.560,
        l = 0.5,
        ks = 0.250,
        θr = 0.078,
        θs = 0.430,
        Ss = 0.0,
        ψe = -1e-3,
    )
    loamspline = RD.SplineConstitutive(loam)
    millerloam = RD.Case(
        formulation = formulation,
        bdf = bdf,
        soil = loamspline,
        Δz = 0.0125,
        Δztotal = 5.0,
        tend = 2.25,
        dt = 0.75,
        ψ0 = RD.InitialHydrostatic(watertable = 0.0),
        topboundary = RD.HeadBoundary(0.1, loamspline),
        bottomboundary = RD.HeadBoundary(0.0, loamspline),
        forcing = nothing,
    )
    return millerloam
end


case = create_millerloam(RD.ReducedDAE(), RD.BDF2())
solver = RD.NewtonSolver(
    RD.LinearSolverThomas(case.parameters),
    #RD.LinearSolverLU(case.parameters),
    relax = RD.SimpleLineSearch(),
    maxiter = 150,
    abstol = 1e-7,
    reltol = 1e-7,
)
timestepper = RD.FixedTimeStepper(3.0)#(Δt0=1e-2, Δtmin=1e-2, Δtmax=1.0e-2)
#timestepper = RD.AdaptiveTimeStepper(Δt0=1e-3, Δtmin=1e-6, Δtmax=1.0e-2)
model = RD.Model(case.parameters, case.ψ0, solver, case.tspan, case.saveat, timestepper)
RD.run!(model)

plot!(model.saved[:, end])

##

function create_millerclayloam(formulation, bdf)
    clayloam = RD.ModifiedMualemVanGenuchten(
        a = 1.900,
        n = 1.310,
        l = 0.5,
        ks = 0.062,
        θr = 0.095,
        θs = 0.410,
        ψe = -1e-3,
        Ss = 0.0,
    )
    clayloamspline = RD.SplineConstitutive(clayloam)

    millerclayloam = RD.Case(
        formulation = formulation,
        bdf = bdf,
        soil = clayloam,#spline,
        Δz = 0.00625,
        Δztotal = 2.0,
        tend = 1.0,
        dt = 1.0,
        ψ0 = RD.InitialHydrostatic(watertable = 0.0),
        topboundary = RD.HeadBoundary(0.1, clayloamspline),
        bottomboundary = RD.HeadBoundary(0.0, clayloamspline),
        forcing = nothing,
    )
    return millerclayloam
end

case = create_millerclayloam(RD.ReducedDAE(), RD.BDF1())
solver = RD.NewtonSolver(
    RD.LinearSolverThomas(case.parameters),
    #RD.LinearSolverLU(case.parameters),
    relax = RD.SimpleLineSearch(),
    maxiter = 150,
    abstol = 1e-7,
    reltol = 1e-7,
)
timestepper = RD.FixedTimeStepper(1e-2)#(Δt0=1e-2, Δtmin=1e-2, Δtmax=1.0e-2)
#timestepper = RD.AdaptiveTimeStepper(Δt0=1e-3, Δtmin=1e-6, Δtmax=1.0e-2)
model = RD.Model(case.parameters, case.ψ0, solver, case.tspan, case.saveat, timestepper)
RD.run!(model)

plot!(model.saved[:, end])

##



##

function run(formulations)
    models = []
    for (formulation, bdf, timestepper) in formulations
        case = create_millersand(formulation, bdf)
        # Solver is formulation dependent: DAE Jacobian is larger.
        solver = RD.NewtonSolver(
            RD.LinearSolverLU(case.parameters),
            relax = RD.ScalarRelaxation(0.0),
            maxiter = 100,
            abstol = 1e-6,
            reltol = 1e-6,
        )
        model =
            RD.Model(case.parameters, case.ψ0, solver, case.tspan, case.saveat, timestepper)
        RD.run!(model)
        push!(models, model)
    end
    return models
end

models = run((
#    (RD.HeadBased(), RD.BDF1(), RD.AdaptiveTimeStepper(Δt0=1.0)),
    (RD.MixedDAE(), RD.BDF1(), RD.AdaptiveTimeStepper(Δt0 = 1.0)),
#    (RD.ReducedDAE(), RD.BDF1(), RD.AdaptiveTimeStepper(Δt0=1.0)),
#    (RD.MixedDAE(), RD.BDF2(), RD.AdaptiveTimeStepper(Δt0=1.0)),
#    (RD.ReducedDAE(), RD.BDF2(), RD.AdaptiveTimeStepper(Δt0=1.0)),
#    (RD.ReducedDAE(), RD.BDF3(), RD.AdaptiveTimeStepper(Δt0=1.0)),
))

p = plot()
plot!(p, models[1].saved[:, end])
#plot!(p, models[2].saved[:, end])
#plot!(p, models[3].saved[:, end])
#plot!(p, models[4].saved[:, end])
#plot!(p, models[5].saved[:, end])

models[1].solver.njacobian
models[2].solver.njacobian
models[3].solver.njacobian
models[4].solver.njacobian
models[5].solver.njacobian


models[1].solver.nresidual
models[2].solver.nresidual
models[3].solver.nresidual
models[4].solver.nresidual
models[5].solver.nresidual

