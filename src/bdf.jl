mutable struct BDFWorkSpace
    max_order::Int  # Maximum BDF order
    u_prev::Matrix{Float64}  # previous states
    Δt_prev::Vector{Float64}  # previous time steps
    a::Vector{Float64}  # BDF coefficients
    order::Int
end

function BDFWorkSpace(parameters::RichardsParameters)
    n = nunknown(parameters)
    order = max_order(parameters.bdforder)
    return BDFWorkSpace(
        order,
        zeros(order, n),
        zeros(order),
        zeros(order + 1),
        1,  # always initialize with first order
    )
end


function set_bdf_coefficients!(bdf::BDFWorkSpace, Δt)
    order = bdf.order
    a = bdf.a
    Δt_prev = bdf.Δt_prev

    if order == 1
        # t[1] = 0, t[2] = -Δt
        a[1] = 1 / Δt
        a[2] = -1 / Δt
    elseif order == 2
        τ₁ = Δt
        τ₂ = Δt + Δt_prev[1]
        a[1] = (τ₁ + τ₂) / (τ₁ * τ₂)
        a[2] = -τ₂ / (τ₁ * (τ₂ - τ₁))
        a[3] = τ₁ / (τ₂ * (τ₂ - τ₁))
    elseif order == 3
        τ₁ = Δt
        τ₂ = Δt + Δt_prev[1]
        τ₃ = τ₂ + Δt_prev[2]
        a[1] = (τ₁*τ₂ + τ₁*τ₃ + τ₂*τ₃) / (τ₁ * τ₂ * τ₃)
        a[2] = -τ₂ * τ₃ / (τ₁ * (τ₂ - τ₁) * (τ₃ - τ₁))
        a[3] = τ₁ * τ₃ / (τ₂ * (τ₂ - τ₁) * (τ₃ - τ₂))
        a[4] = -τ₁ * τ₂ / (τ₃ * (τ₃ - τ₁) * (τ₃ - τ₂))
    elseif order == 4
        τ₁ = Δt
        τ₂ = Δt + Δt_prev[1]
        τ₃ = τ₂ + Δt_prev[2]
        τ₄ = τ₃ + Δt_prev[3]
        a[1] = (τ₁*τ₂*τ₃ + τ₁*τ₂*τ₄ + τ₁*τ₃*τ₄ + τ₂*τ₃*τ₄) / (τ₁ * τ₂ * τ₃ * τ₄)
        a[2] = -τ₂ * τ₃ * τ₄ / (τ₁ * (τ₂ - τ₁) * (τ₃ - τ₁) * (τ₄ - τ₁))
        a[3] = τ₁ * τ₃ * τ₄ / (τ₂ * (τ₂ - τ₁) * (τ₃ - τ₂) * (τ₄ - τ₂))
        a[4] = -τ₁ * τ₂ * τ₄ / (τ₃ * (τ₃ - τ₁) * (τ₃ - τ₂) * (τ₄ - τ₃))
        a[5] = τ₁ * τ₂ * τ₃ / (τ₄ * (τ₄ - τ₁) * (τ₄ - τ₂) * (τ₄ - τ₃))
    elseif order == 5
        τ₁ = Δt
        τ₂ = Δt + Δt_prev[1]
        τ₃ = τ₂ + Δt_prev[2]
        τ₄ = τ₃ + Δt_prev[3]
        τ₅ = τ₄ + Δt_prev[4]
        s₁₂₃₄ = τ₁*τ₂*τ₃*τ₄
        s₁₂₃₅ = τ₁*τ₂*τ₃*τ₅
        s₁₂₄₅ = τ₁*τ₂*τ₄*τ₅
        s₁₃₄₅ = τ₁*τ₃*τ₄*τ₅
        s₂₃₄₅ = τ₂*τ₃*τ₄*τ₅
        a[1] = (s₁₂₃₄ + s₁₂₃₅ + s₁₂₄₅ + s₁₃₄₅ + s₂₃₄₅) / (τ₁ * τ₂ * τ₃ * τ₄ * τ₅)
        a[2] = -s₂₃₄₅ / (τ₁ * (τ₂ - τ₁) * (τ₃ - τ₁) * (τ₄ - τ₁) * (τ₅ - τ₁))
        a[3] = s₁₃₄₅ / (τ₂ * (τ₂ - τ₁) * (τ₃ - τ₂) * (τ₄ - τ₂) * (τ₅ - τ₂))
        a[4] = -s₁₂₄₅ / (τ₃ * (τ₃ - τ₁) * (τ₃ - τ₂) * (τ₄ - τ₃) * (τ₅ - τ₃))
        a[5] = s₁₂₃₅ / (τ₄ * (τ₄ - τ₁) * (τ₄ - τ₂) * (τ₄ - τ₃) * (τ₅ - τ₄))
        a[6] = -s₁₂₃₄ / (τ₅ * (τ₅ - τ₁) * (τ₅ - τ₂) * (τ₅ - τ₃) * (τ₅ - τ₄))
    else
        error("BDF order $order not implemented")
    end
    return
end

function rotate_history!(bdf::BDFWorkSpace, u, Δt)
    # Shift columns right
    for j = bdf.max_order:-1:2
        bdf.u_prev[j, :] .= bdf.u_prev[j-1, :]
        bdf.Δt_prev[j] = bdf.Δt_prev[j-1]
    end
    bdf.u_prev[1, :] .= u
    bdf.Δt_prev[1] = Δt
    # Increase order if it hasn't reached maximum order yet.
    bdf.order = min(bdf.order + 1, bdf.max_order)
    return
end
