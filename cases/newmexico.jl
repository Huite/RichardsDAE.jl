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
    newmexico = RD.Case(
        formulation = formulation,
        bdf = bdf,
        soil = sand,
        Δz = 2.5e-3,
        Δztotal = 0.3,
        tend = 0.25,
        save_dt = 0.25,
        ψ0 = RD.InitialConstant(-1.0e1),
        topboundary = RD.HeadBoundary(-7.5e-1, sand),
        bottomboundary = RD.HeadBoundary(-1.0e1, sand),
        forcing = nothing,
    )
    return newmexico
end

function rmse(model, refmodel)
    ψ = model.saved[:, end]
    ψref = refmodel.saved[:, end]
    error = ψ .- ψref
    return sqrt(mean(error .^ 2))
end

function run(formulations)
    models = []
    for (formulation, bdf, timestepper) in formulations
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

# Run fixed time step benchmarks.

dts = 10 .^ collect(-1:-0.5:-4)
formulations = []
for bdf in (RD.BDF1(), RD.BDF2(), RD.BDF3())
    for dt in dts
        formulation = (RD.ReducedDAE(), bdf, RD.FixedTimeStepper(dt))
        push!(formulations, formulation)
    end
end
# Add a reference run with fine time discretization.
push!(formulations, (RD.ReducedDAE(), RD.BDF1(), RD.FixedTimeStepper(1.0e-6)))
models = run(formulations)

# Collect the output, store the final heads and create a work-error plot.

data = Dict{String,Vector{Float64}}()
labels = ("BDF1", "BDF2", "BDF3")
n = 120
for (model, (label, dt)) in zip(models, Iterators.product(labels, dts))
    data["Reduced-$(label) (dt=$(dt))"] = model.saved[1:n, end]
end
data["Reduced-BDF1 (dt=0.000001)"] = models[end].saved[1:n, end]
headdf = DataFrame(data)
CSV.write("cases/output/newmexico-final-head.csv", headdf)

refmodel = models[end]
model_rmses = [rmse(model, refmodel) for model in models]
njac = [model.solver.njacobian for model in models]
nres = [model.solver.nresidual for model in models]
nsolve = [model.solver.nlinsolve for model in models]
model_work = nsolve
ntimes = length(dts)
nbdf = 3
work = reshape(model_work[1:(end-1)], (ntimes, nbdf))
error = reshape(model_rmses[1:(end-1)], (ntimes, nbdf))

const COLORS = RD.okabe_ito_colors()
p = plot(
    xaxis = :log10,
    #    yaxis = :log10,
    ylabel = "RMSE (m)",
    xlabel = "Work",
)
plot!(p, work[:, 1], error[:, 1], label = "", color = COLORS[:dark_orange], lw = 2)
scatter!(p, work[:, 1], error[:, 1], label = "BDF1", color = COLORS[:dark_orange])
plot!(p, work[:, 2], error[:, 2], label = "", color = COLORS[:light_blue], lw = 2)
scatter!(p, work[:, 2], error[:, 2], label = "BDF2", color = COLORS[:light_blue])
plot!(p, work[:, 3], error[:, 3], label = "", color = COLORS[:green], lw = 2)
scatter!(p, work[:, 3], error[:, 3], label = "BDF3", color = COLORS[:green])
savefig(p, "cases/output/newmexico-errorwork.pdf")
