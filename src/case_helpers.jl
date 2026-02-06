struct Case
    parameters::RichardsParameters
    ψ0::Vector{Float64}
    tspan::Tuple{Float64,Float64}
    saveat::Vector{Float64}
end

function Case(;
    formulation,
    bdf,
    soil,
    Δz,
    Δztotal,
    tend,
    save_dt,
    ψ0,
    bottomboundary,
    topboundary,
    forcing,
)
    n = Int(Δztotal / Δz)
    if isnothing(forcing)
        forcing = MeteorologicalForcing([0.0], [0.0], [0.0])
    end
    saveat = collect(0.0:save_dt:tend)
    return Case(
        RichardsParameters(
            formulation,
            bdf,
            constitutive = fill(soil, n),
            Δz = Δz,
            forcing = forcing,
            bottomboundary = bottomboundary,
            topboundary = topboundary,
        ),
        initialψ(ψ0, Δz, Δztotal, n),
        (0.0, tend),
        saveat,
    )
end

function storage(ψ, parameters::RichardsParameters)
    θ = moisture_content.(ψ, parameters.constitutive)
    #S_elastic = [c.Ss * θi / c.θs for (c, θi) in zip(parameters.constitutive, θ)]
    #return parameters.Δz * (θ + S_elastic)
    return parameters.Δz * θ
end

function waterbalance_dataframe(model)
    # Compute cumulative flows back to reporting steps.
    # One smaller (vertex vs interval): no cumulative yet flow at t=0.
    parameters = (model.parameters)
    n = parameters.n
    S = [storage(ψ, parameters) for ψ in eachcol(model.saved[1:n, :])]
    return DataFrame(
        :t => model.saveat,
        :storage => sum.(S),
        :qbot => model.savedflows[1, :],
        :qtop => model.savedflows[2, :],
    )
end

function massbalance_balance_ratio(waterbalance)
    ΔS = waterbalance.storage[end] - waterbalance.storage[1]
    net_boundary_flux = waterbalance.qbot[end] + waterbalance.qtop[end]
    return net_boundary_flux / ΔS
end

function massbalance_rmse(waterbalance)
    ΔS = diff(waterbalance.storage)
    qb = diff(waterbalance.qbot)
    qt = diff(waterbalance.qtop)
    error = qt + qb - ΔS
    return sqrt(mean(error .^ 2))
end

function okabe_ito_colors()
    # Should be colorblind friendly.
    return Dict(
        :orange => "#E69F00",
        :light_blue => "#56B4E9",
        :green => "#009E73",
        :yellow => "#F0E442",
        :blue => "#0072B2",
        :dark_orange => "#D55E00",
        :pink => "#CC79A7",
        :black => "#000000",
    )
end
