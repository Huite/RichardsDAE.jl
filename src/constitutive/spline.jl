# [core]
struct SplineConstitutive{T} <: ConstitutiveRelationships
    θ::T
    k::T
    nknots::Int
    maxrelerror::Float64
    θs::Float64
    Ss::Float64
end

function Base.show(io::IO, sc::SplineConstitutive)
    splinename = Base.typename(typeof(sc.k)).name
    println(
        io,
        "SplineConstitutive(θ=$splinename, k=$splinename, nknots=$(sc.nknots), maxrelerror=$(sc.maxrelerror), θs=$(sc.θs), Ss=$(sc.Ss))",
    )
end

function conductivity(ψ, hsc::SplineConstitutive)
    return hsc.k(ψ)
end

function dconductivity(ψ, hsc::SplineConstitutive)
    return DataInterpolations.derivative(hsc.k, ψ, 1)
end

function moisture_content(ψ, hsc::SplineConstitutive)
    return hsc.θ(ψ)
end

function specific_moisture_capacity(ψ, hsc::SplineConstitutive)
    return DataInterpolations.derivative(hsc.θ, ψ, 1)
end

function logknots(ψmin, ψe, offset, nknots)
    # Use log-spaced interpolation points.
    # Approach saturation at ψe, but not quite.
    return vcat(
        -exp10.(range(log10(abs(ψmin)), log10(abs(ψe - offset)), length = nknots)),
        ψe .+ [offset, 2 * offset],
    )
end

"""
   SplineConstitutive(parameters; kwargs...) -> spline

Create a spline-based approximation of soil hydraulic functions for numerically stable 
evaluation in Richards equation solvers.

This function fits PCHIP (Piecewise Cubic Hermite Interpolating Polynomial) splines to 
the moisture content θ(ψ) and hydraulic conductivity K(ψ) functions, providing smooth 
derivatives that avoid numerical issues with steep gradients near saturation.

# Arguments
- `parameters`: Soil hydraulic parameter struct (e.g., `MualemVanGenuchten`, `Haverkamp`)

# Keyword Arguments
- `relative_error=1e-4`: Target maximum relative error for spline approximation
- `offset=1e-3`: Offset from air entry pressure ψe to avoid singularities
- `ψmin=-1e4`: Initial guess for minimum pressure head [L]. Automatically adjusted if 
  conductivity becomes negligible (< 1e-12)
- `nknots=100`: Initial number of knot points for spline interpolation
- `iter=10`: Maximum refinement iterations to achieve target accuracy

# Returns
- `spline`: `SplineConstitutive` object with interpolated θ(ψ) and K(ψ) functions

# Details
The algorithm uses logarithmically-spaced knot points concentrated near the air entry 
pressure where gradients are steepest. If the target accuracy is not met, the number 
of knots is increased by 50% and the process repeats.

The minimum pressure `ψmin` is automatically adjusted to avoid regions where hydraulic 
conductivity approaches zero, preventing division-by-zero errors in relative error 
calculations.

The maximum relative error is returned for diagnostic purposes under `.maxrelerror`.
"""
function SplineConstitutive(
    parameters;
    relative_error = 1e-4,
    offset = 1e-3,
    ψmin = -1e4,
    nknots = 100,
    iter = 10,
)
    local θspline, kspline, t, maxerror

    if hasproperty(parameters, :ψe)
        ψe = parameters.ψe
    else
        ψe = 0.0
    end

    # Test ψmin; it's unsuitable if it returns zero since we'll end
    # end up dividing by it later. Should converge rapidly.
    for _ = 1:100
        if conductivity(ψmin, parameters) > 1e-12
            break
        else
            ψmin *= 0.5
        end
    end

    for _ = 1:iter
        t = logknots(ψmin, ψe, offset, nknots)
        θ = moisture_content.(t, Ref(parameters))
        θspline = DataInterpolations.PCHIPInterpolation(
            θ,
            t,
            extrapolation = DataInterpolations.ExtrapolationType.Constant,
        )
        k = conductivity.(t, Ref(parameters))
        kspline = DataInterpolations.PCHIPInterpolation(
            k,
            t,
            extrapolation = DataInterpolations.ExtrapolationType.Constant,
        )

        # Test at ten times the density.
        ψtest =
            -exp10.(
                range(
                    log10(abs(ψmin)),
                    log10(abs(ψe - 2.0 * offset)),
                    length = nknots * 10,
                ),
            )
        θref = moisture_content.(ψtest, Ref(parameters))
        kref = conductivity.(ψtest, Ref(parameters))

        θerror = @. (θspline(ψtest) - θref) / θref
        kerror = @. (kspline(ψtest) - kref) / kref

        maxerror = max(maximum(abs.(θerror)), maximum(abs.(kerror)))

        if maxerror < relative_error
            break
        else
            nknots = ceil(Int, nknots * 1.5)
        end
    end
    return SplineConstitutive(
        θspline,
        kspline,
        nknots,
        maxerror,
        parameters.θs,
        parameters.Ss,
    )
end
