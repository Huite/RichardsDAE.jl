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

# Modified Mualem–van Genuchten relations (Vogel 2000, Ippisch 2006)

# [core]
struct ModifiedMualemVanGenuchten <: ConstitutiveRelationships
    a::Float64      # van Genuchten a [1/L]
    n::Float64      # pore‑size distribution parameter
    m::Float64      # usually m = 1 − 1/n
    l::Float64      # pore‑connectivity (Mualem τ), default ~= 0.5
    ks::Float64     # saturated hydraulic conductivity [L/T]
    θs::Float64     # saturated water content
    θr::Float64     # residual water content
    ψe::Float64     # air‑entry suction (> 0)  [L]
    Sc::Float64     # cut‑off saturation factor
    Ss::Float64
    function ModifiedMualemVanGenuchten(; a, n, m = nothing, l, ks, θs, θr, ψe, Ss)
        if isnothing(m)
            m = 1 - 1 / n
        end
        Sc = (1 + abs(a * ψe)^n)^(-m)
        new(a, n, m, l, ks, θs, θr, ψe, Sc, Ss)
    end
end

# [core]
function effective_saturation(ψ, mvg::ModifiedMualemVanGenuchten)
    (; a, n, m, ψe, Sc) = mvg
    return ifelse(ψ > ψe, 1.0, (1 / Sc) * (1 + (a * abs(ψ))^n)^(-m))
end

# [jacobian]
function deffective_saturation(ψ, mvg::ModifiedMualemVanGenuchten)
    (; a, n, m, Sc, ψe) = mvg
    absψ = abs(ψ)
    f = 1 + (a * absψ)^n
    derivative = (m * n * a^n / Sc) * absψ^(n - 1) * f^(-m - 1)
    return ifelse(ψ > ψe, 0.0, derivative)
end

# [core]
function moisture_content(ψ, mvg::ModifiedMualemVanGenuchten)
    (; θs, θr) = mvg
    return θr + effective_saturation(ψ, mvg) * (θs - θr)
end

# [core]
"""Specific moisture capacity C = dθ/dψ"""
function specific_moisture_capacity(ψ, mvg::ModifiedMualemVanGenuchten)
    (; θs, θr) = mvg
    return deffective_saturation(ψ, mvg) * (θs - θr)
end


# [core]
# Helper function
@inline function _F(x, m)
    return 1 - (1 - x^(1 / m))^m
end

# [core]
function conductivity(ψ, mvg::ModifiedMualemVanGenuchten)
    (; ks, l, m, Sc) = mvg
    Se = effective_saturation(ψ, mvg)
    Fc = _F(Sc, m)
    Se_ = Se * Sc
    F = _F(Se_, m)
    kr = ifelse(Se >= 1, 1.0, Se^l * (F / Fc)^2)
    return ks * kr
end

# [jacobian]
"""∂K/∂ψ for Jacobian"""
function dconductivity(ψ, mvg::ModifiedMualemVanGenuchten)
    (; ks, l, m, Sc, ψe) = mvg
    Se = effective_saturation(ψ, mvg)
    dSe = deffective_saturation(ψ, mvg)
    Fc = _F(Sc, m)
    Se_ = Se * Sc
    G = 1 - Se_^(1 / m)
    F = 1 - G^m
    dF_dSe = Sc * G^(m - 1) * Se_^(1 / m - 1)
    dkr_dSe = l * Se^(l - 1) * (F / Fc)^2 + 2 * Se^l * (F / Fc) * (dF_dSe / Fc)
    derivative = ks * dkr_dSe * dSe
    return ifelse(ψ > ψe, 0.0, derivative)
end
