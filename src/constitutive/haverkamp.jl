"""Haverkamp constitutive relationship."""
@kwdef struct Haverkamp <: ConstitutiveRelationships
    a::Float64
    B::Float64
    y::Float64
    A::Float64
    ks::Float64
    θs::Float64
    θr::Float64
    Ss::Float64
end

function conductivity(ψ, h::Haverkamp)
    return h.ks * h.A / (h.A + abs(min(ψ, 0.0))^h.y)
end

"""dK/dψ for Newton formulation."""
function dconductivity(ψ, h::Haverkamp)
    return -h.A * h.ks * h.y * min(ψ, 0.0) * abs(min(ψ, 0.0))^(h.y - 2) /
           (h.A + abs(min(ψ, 0.0))^h.y)^2
end

function moisture_content(ψ, h::Haverkamp)
    return h.a * (h.θs - h.θr) / (h.a + abs(ψ)^h.B) + h.θr
end

function specific_moisture_capacity(ψ, h::Haverkamp)
    return h.a * h.B * (h.θs - h.θr) * abs(ψ)^(h.B - 1) / (h.a + abs(ψ)^h.B)^2
end

function dspecific_moisture_capacity(ψ, h::Haverkamp)
    (; a, B, θs, θr) = h
    A = a * (θs - θr)
    absψ = abs(ψ)
    num1 = (B - 1) * absψ^(B - 2) * (a + absψ^B)^2
    num2 = 2 * B * ψ * absψ^(2B - 3) * (a + absψ^B)
    denom = (a + absψ^B)^4
    return -A * B * (num1 - num2) / denom
end

function aqueous_saturation(ψ, h::Haverkamp)
    return moisture_content(ψ, h) / h.θs
end

function daqueous_saturation_dψ(ψ, h::Haverkamp)
    return specific_moisture_capacity(ψ, h) / h.θs
end
