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

function generate_formulations(scheme, dts)
    formulations = []
    for bdf in (RD.BDF1(), RD.BDF2(), RD.BDF3())
        for dt in dts
            formulation = (scheme, bdf, RD.FixedTimeStepper(dt))
            push!(formulations, formulation)
        end
    end
    return formulations
end

# Run fixed time step benchmarks.
# Add a reference run with fine time discretization.

dts = 10 .^ collect(-1:-0.5:-4)
reduced_formulations = generate_formulations(RD.ReducedDAE(), dts)
mixed_formulations = generate_formulations(RD.MixedDAE(), dts)
reduced_models = run(reduced_formulations)
mixed_models = run(mixed_formulations)
reference_model = run([(RD.ReducedDAE(), RD.BDF1(), RD.FixedTimeStepper(1.0e-6))])[1]

# Collect the output, store the final heads and create a work-error plot.

data = Dict{String,Vector{Float64}}()
labels = ("BDF1", "BDF2", "BDF3")
n = 120
for (model, (label, dt)) in zip(reduced_models, Iterators.product(labels, dts))
    data["Reduced-$(label) (dt=$(dt))"] = model.saved[1:n, end]
end
for (model, (label, dt)) in zip(mixed_models, Iterators.product(labels, dts))
    data["Mixed-$(label) (dt=$(dt))"] = model.saved[1:n, end]
end
data["Reduced-BDF1 (dt=0.000001)"] = reference_model.saved[1:n, end]

headdf = DataFrame(data)
CSV.write("cases/output/newmexico-final-head.csv", headdf)

reduced_model_rmses = [rmse(model, reference_model) for model in reduced_models]
mixed_model_rmses = [rmse(model, reference_model) for model in mixed_models]
ntimes = length(dts)
nbdf = 3

const COLORS = RD.okabe_ito_colors()
p = plot(
    xaxis = :log10,
    ylabel = "RMSE (m)",
    xlabel = "Δt (d)",
    xticks = ([0.0001, 0.001, 0.01, 0.1], ["0.0001", "0.001", "0.01", "0.1"]),
)
reduced_error = reshape(reduced_model_rmses, (ntimes, nbdf))
mixed_error = reshape(mixed_model_rmses, (ntimes, nbdf))
scatter!(
    p,
    dts,
    mixed_error[:, 1],
    label = "Mixed-BDF1",
    color = COLORS[:dark_orange],
    markersize = 5,
    markerstrokewidth = 0,
)
scatter!(
    p,
    dts,
    mixed_error[:, 2],
    label = "Mixed-BDF2",
    color = COLORS[:light_blue],
    markersize = 5,
    markerstrokewidth = 0,
)
scatter!(
    p,
    dts,
    mixed_error[:, 3],
    label = "Mixed-BDF3",
    color = COLORS[:green],
    markersize = 5,
    markerstrokewidth = 0,
)
plot!(
    p,
    dts,
    reduced_error[:, 1],
    label = "Reduced-BDF1",
    color = COLORS[:dark_orange],
    lw = 2,
)
plot!(
    p,
    dts,
    reduced_error[:, 2],
    label = "Reduced-BDF2",
    color = COLORS[:light_blue],
    lw = 2,
)
plot!(p, dts, reduced_error[:, 3], label = "Reduced-BDF3", color = COLORS[:green], lw = 2)
savefig(p, "cases/output/newmexico-error-timestep.pdf")
