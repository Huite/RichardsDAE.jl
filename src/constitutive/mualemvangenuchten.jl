struct MualemVanGenuchten <: ConstitutiveRelationships
    a::Float64      # van Genuchten a [1/L]
    n::Float64      # pore‑size distribution parameter
    m::Float64      # usually m = 1 − 1/n
    l::Float64      # pore‑connectivity (Mualem τ), default ~= 0.5
    ks::Float64     # saturated hydraulic conductivity [L/T]
    θs::Float64     # saturated water content
    θr::Float64     # residual water content
    Ss::Float64
    function MualemVanGenuchten(; a, n, m = nothing, l, ks, θs, θr, Ss)
        if isnothing(m)
            m = 1 - 1 / n
        end
        new(a, n, m, l, ks, θs, θr, Ss)
    end
end

function effective_saturation(ψ, mvg::MualemVanGenuchten)
    (; a, n, m) = mvg
    Se = (1 + (a * abs(ψ))^n)^(-m)
    return ifelse(ψ > 0, 1, Se)
end

function deffective_saturation(ψ, mvg::MualemVanGenuchten)
    (; a, n, m) = mvg
    absψ = abs(ψ)
    base = 1 + (a * absψ)^n
    dSe = m * n * a^n * absψ^(n - 1) * base^(-m - 1)
    return ifelse(ψ > 0, 0, dSe)
end

function moisture_content(ψ, mvg::MualemVanGenuchten)
    (; θs, θr) = mvg
    Se = effective_saturation(ψ, mvg)
    return θr + Se * (θs - θr)
end

function specific_moisture_capacity(ψ, mvg::MualemVanGenuchten)
    (; θs, θr) = mvg
    return deffective_saturation(ψ, mvg) * (θs - θr)
end

function conductivity(ψ, mvg::MualemVanGenuchten)
    (; ks, l, m) = mvg
    Se = effective_saturation(ψ, mvg)
    return ks * Se^l * (1 - (1 - Se^(1 / m))^m)^2
end

"""dK/dψ for Newton formulation."""
function dconductivity(ψ, mvg::MualemVanGenuchten)
    (; ks, l, m) = mvg
    Se = effective_saturation(ψ, mvg)
    dSe_dψ = deffective_saturation(ψ, mvg)

    # Term 1: derivative of Se^l
    term1 = l * Se^(l - 1) * dSe_dψ

    # Term 2: derivative of (1 - (1 - Se^(1/m))^m)^2
    inner = 1 - Se^(1 / m)
    inner_der = -(1 / m) * Se^(1 / m - 1) * dSe_dψ
    outer = 1 - inner^m
    outer_der = -m * inner^(m - 1) * inner_der

    return ifelse(ψ > 0, 0.0, ks * (term1 * (1 - inner^m)^2 + 2 * Se^l * outer * outer_der))
end
