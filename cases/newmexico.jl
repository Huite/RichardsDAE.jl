import RichardsDAE as RD
using Plots
using DataFrames
using CSV
using Statistics


function create_newmexico(formulation, bdf)
    sand = RD.MualemVanGenuchten(
        a = 3.35,
        n = 2.0,
        l = 0.5,
        ks = 7.97,
        θr = 0.102,
        θs = 0.368,
        Ss = 0.0,
    )
    sandspline = RD.SplineConstitutive(sand)
    newmexico = RD.Case(
        formulation = formulation,
        bdf = bdf,
        soil = sand,
        Δz = 2.5e-3,
        Δztotal = 0.3,
        tend = 0.25,
        dt = 0.25,
        ψ0 = RD.InitialConstant(-1.0e1),
        topboundary = RD.HeadBoundary(-7.5e-1, sand),
        bottomboundary = RD.HeadBoundary(-1.0e1, sand),
        forcing = nothing,
    )
    return newmexico
end


function run(formulations)
    models = []
    for (formulation, bdf, timestepper) in formulations
        @show formulation
        case = create_newmexico(formulation, bdf)
        solver = RD.NewtonSolver(
            RD.LinearSolverLU(case.parameters),
            relax = RD.SimpleLineSearch(),
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

dts = (0.1, 0.01, 0.001, 0.001, 0.0001)
e = collect(-1:-0.5:-4)
dts = 10 .^ e
formulations = []
for bdf in (RD.BDF1(), RD.BDF2(), RD.BDF3(), RD.BDF4(), RD.BDF5())
    for dt in dts
        formulation = (RD.ReducedDAE(), bdf, RD.FixedTimeStepper(dt))
        push!(formulations, formulation)
    end
end
push!(formulations, (RD.ReducedDAE(), RD.BDF1(), RD.FixedTimeStepper(1.0e-6)))

models = run(formulations)


function rmse(model, refmodel)
    ψ = model.saved[:, end]
    ψref = refmodel.saved[:, end]
    error = ψ .- ψref
    return sqrt(mean(error .^ 2))
end


refmodel = models[end]
model_rmses = [rmse(model, refmodel) for model in models]
model_work = [model.solver.njacobian for model in models]
ntimes = length(dts)
nbdf = 5

work = reshape(model_work[1:(end-1)], (ntimes, nbdf))
error = reshape(model_rmses[1:(end-1)], (ntimes, nbdf))

p = scatter(
    work[:, 1],
    error[:, 1],
    label = "BDF1",
    xaxis = :log10,
    yaxis = :log10,
    ylabel = "RMSE",
    xlabel = "Work",
)
plot!(p, work[:, 1], error[:, 1], label = "BDF1")
scatter!(p, work[:, 2], error[:, 2], label = "BDF2")
plot!(p, work[:, 2], error[:, 2], label = "BDF2")
scatter!(p, work[:, 3], error[:, 3], label = "BDF3")
plot!(p, work[:, 3], error[:, 3], label = "BDF3")


p = scatter(
    work[:, 1],
    error[:, 1],
    label = "BDF1",
    xaxis = :log10,
    ylabel = "RMSE",
    xlabel = "Work",
)
plot!(p, work[:, 1], error[:, 1], label = "BDF1")
scatter!(p, work[:, 2], error[:, 2], label = "BDF2")
plot!(p, work[:, 2], error[:, 2], label = "BDF2")
scatter!(p, work[:, 3], error[:, 3], label = "BDF3")
plot!(p, work[:, 3], error[:, 3], label = "BDF3")

scatter!(p, model_work[6:10], model_rmses[6:10], label = "BDF2")
scatter!(p, model_work[11:15], model_rmses[11:15], label = "BDF3")
#scatter!(p, model_work[16:20], model_rmses[16:20], label="BDF4")
#scatter!(p, model_work[21:25], model_rmses[21:25], label="BDF5")

plot!(p, model_work[1:5], model_rmses[1:5], label = "BDF1")
plot!(p, model_work[6:10], model_rmses[6:10], label = "BDF2")
plot!(p, model_work[11:15], model_rmses[11:15], label = "BDF3")
#plot!(p, model_work[16:20], model_rmses[16:20], label="BDF4")
#plot!(p, model_work[21:25], model_rmses[21:25], label="BDF5")





case = create_newmexico(RD.ReducedDAE(), RD.BDF1())
solver = RD.NewtonSolver(
    RD.LinearSolverThomas(case.parameters),
    #RD.LinearSolverLU(case.parameters),
    relax = RD.SimpleLineSearch(),
    maxiter = 150,
    abstol = 1e-8,
    reltol = 1e-8,
)
timestepper = RD.FixedTimeStepper(0.1)#(Δt0=1e-2, Δtmin=1e-2, Δtmax=1.0e-2)
#timestepper = RD.AdaptiveTimeStepper(Δt0=1e-3, Δtmin=1e-6, Δtmax=1.0e-2)

model = RD.Model(case.parameters, case.ψ0, solver, case.tspan, case.saveat, timestepper)
RD.run!(model)


p = plot(models[1].saved[:, end])
plot!(p, models[end].saved[:, end])
plot!(p, models[14].saved[:, end])
