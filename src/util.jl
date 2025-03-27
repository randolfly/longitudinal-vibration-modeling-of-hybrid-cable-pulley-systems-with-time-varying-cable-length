using .HybridCableModel

# =============== COMMON FUNCTIONS ====================

function Φ(i::Integer, x̃::Float64)
    return sin(i * π * x̃)
end

function ∂ₓΦ(i::Integer, x̃::Float64)
    return i * π * cos(i * π * x̃)
end

function ∂ₓₓΦ(i::Integer, x̃::Float64)
    return -(i * π)^2 * sin(i * π * x̃)
end

function ∫Φ(i::Integer)
    return (1 + (-1)^(1 + i)) / (i * pi)
end

function ∫ΦΦ(i::Integer, j::Integer)
    if i == j
        return 1 / 2
    else
        return 0
    end
end

function ∫∂ₓΦ∂ₓΦ(i::Integer, j::Integer)
    if i == j
        return 1 / 2 * (i * pi)^2
    else
        return 0
    end
end

function ∫x∂ₓΦ(i::Integer)
    return (-1 + (-1)^(i)) / (i * pi)
end

function ∫x²∂ₓₓΦ(i::Integer)
    return (2 - 2 * (-1)^i + (-1)^i * i^2 * pi^2) / (i * pi)
end

function ∫x∂ₓΦΦ(i::Integer, j::Integer)
    if i == j
        return -1 / 4
    else
        return ((-1)^(i + j) * i * j) / (i^2 - j^2)
    end
end

function ∫x²∂ₓₓΦΦ(i::Integer, j::Integer)
    if i == j
        return 1 / 4 - i^2 * pi^2 / 6
    else
        return -(4 * (-1)^(i + j) * i^3 * j) / ((i^2 - j^2)^2)
    end
end

# get friction force
function f(∂ₜθ::Float64, T::Float64, C::Float64, tanh_ratio::Float64=1000.0)
    # _f = T * sign(∂ₜθ) + C * ∂ₜθ
    _f = T * tanh(tanh_ratio * ∂ₜθ) + C * ∂ₜθ
    return _f
end

function fm(dxe::Float64, Tm::Float64, Cm::Float64, tanh_ratio::Float64=1000.0)
    # _f = -T * sign(dxe) - C * dxe
    _f = -Tm * tanh(tanh_ratio * dxe) - Cm * dxe
    return _f
end


function numerical_derivative(sol, t)
    ForwardDiff.derivative(t -> sol(t), t)
end


function get_cable_force(mp::ModelParam, xe::Float64, ddxe::Float64)
    cable_force = mp.k * (mp.L - xe) - mp.m * ddxe
    return cable_force
end

# =============== SIMPLE MODEL ====================

## u function

# u(l,t) = ∑ Φᵢ(l) xᵢ
function u(mt::SimpleModel, mp::ModelParam, x::Vector{Float64}, l::Float64)
    xe = x[end]
    _u = 0
    for i = 1:mp.N
        _u += Φ(i, l / xe) * x[i]
    end
    return _u
end

function dθ(mt::SimpleModel, mp::ModelParam, dxe::Float64)
    _dθ = (-dxe) / mp.rd
    return _dθ
end

# =============== PARTIAL MODEL ====================

## u function

# u(l,t) = ∑ Φᵢ(l) xᵢ
function u(mt::PartialModel, mp::ModelParam, x::Vector{Float64}, l::Float64)
    xe = x[end]
    _u = 0
    for i = 1:mp.N
        _u += Φ(i, l / xe) * x[i]
    end
    return _u
end

function ∂ₓu(mt::PartialModel, mp::ModelParam, x::Vector{Float64}, l::Float64)
    xe = x[end]
    _∂ₓu = 0
    x̃ = l / xe

    for i = 1:mp.N
        _∂ₓu += ∂ₓΦ(i, x̃) * x[i] / xe
    end

    return _∂ₓu
end

function ∂ₜu(mt::PartialModel, mp::ModelParam, x::Vector{Float64}, dx::Vector{Float64}, l::Float64)
    _∂ₜu = 0
    xe = x[end]
    x̃ = l / xe

    for i in 1:mp.N
        _∂ₜu += Φ(i, x̃) * dx[i]
    end

    return _∂ₜu
end

function ∂ₓₓu(mt::PartialModel, mp::ModelParam, x::Vector{Float64}, l::Float64)
    xe = x[end]
    x̃ = l / xe

    _∂ₓₓu = 0
    for i in 1:mp.N
        _∂ₓₓu += ∂ₓₓΦ(i, x̃) * x[i] / xe^2
    end

    return _∂ₓₓu
end

function ∂ₓₜu(mt::PartialModel, mp::ModelParam, x::Vector{Float64}, dx::Vector{Float64}, l::Float64)
    _∂ₓₜu = 0
    xe = x[end]
    dxe = dx[end]
    x̃ = l / xe

    for i in 1:mp.N
        _∂ₓₜu += ∂ₓΦ(i, x̃) * dx[i] / xe
    end
    return _∂ₓₜu
end

function ∂ₓₑuₚ(mt::PartialModel, mp::ModelParam, pulley_id::Integer, x::Vector{Float64})
    _∂ₓₑuₚ = 0
    xe = x[end]
    for i in 1:mp.N
        _∂ₓₑuₚ += x[i] * ∂ₓΦ(i, 1 - mp.lp[pulley_id] / xe)
    end
    _∂ₓₑuₚ *= (mp.lp[pulley_id] / xe^2)
    return _∂ₓₑuₚ
end

## θ function

function θ(mt::PartialModel, mp::ModelParam, xe::Float64)
    _θ = (-xe + mp.L) / mp.rd
    return _θ
end

function θp(mt::PartialModel, mp::ModelParam, pulley_id::Integer, x::Vector{Float64})
    xe = x[end]
    _θp = (mp.L - xe + u(mt, mp, x, xe - mp.lp[pulley_id])) / mp.rp[pulley_id]
    return _θp
end

function dθ(mt::PartialModel, mp::ModelParam, dxe::Float64)
    _dθ = (-dxe) / mp.rd
    return _dθ
end

function dθp(mt::PartialModel, mp::ModelParam, pulley_id::Integer, x::Vector{Float64}, dx::Vector{Float64})
    xe = x[end]
    dxe = dx[end]
    _∂ₜθ = (-dxe + ∂ₜu(mt, mp, x, dx, xe - mp.lp[pulley_id]) + ∂ₓu(mt, mp, x, xe - mp.lp[pulley_id]) * dxe) / mp.rp[pulley_id]
    return _∂ₜθ
end

function s(mt::PartialModel, mp::ModelParam, x::Vector{Float64})
    xe = x[end]
    _θ = θ(mt, mp, xe)
    _θₚ = Vector{Float64}(undef, mp.pulley_num)
    for i in 1:mp.pulley_num
        _θₚ[i] = θp(mt, mp, i, x)
    end
    _s = [x; _θ; _θₚ]
    return _s
end

# construct ds = [dq; dθ; dθ₁, dθ₂,...,dθq]
function ds(mt::PartialModel, mp::ModelParam, x::Vector{Float64}, dx::Vector{Float64})
    xe = x[end]
    dxe = dx[end]

    _dθ = dθ(mt, mp, dxe)
    _dθₚ = Vector{Float64}(undef, mp.pulley_num)
    for i in 1:mp.pulley_num
        _dθₚ[i] = dθp(mt, mp, i, x, dx)
    end
    _ds = [dx; _dθ; _dθₚ]
    return _ds
end

# =============== FULL MODEL ====================

## u function

# u(l,t) = ∑ Φᵢ(l) xᵢ
function u(mt::FullModel, mp::ModelParam, x::Vector{Float64}, l::Float64)
    xe = x[end]
    _u = 0
    for i = 1:mp.N
        _u += Φ(i, l / xe) * x[i]
    end
    return _u
end

function ∂ₓu(mt::FullModel, mp::ModelParam, x::Vector{Float64}, l::Float64)
    xe = x[end]
    _∂ₓu = 0
    x̃ = l / xe

    for i = 1:mp.N
        _∂ₓu += ∂ₓΦ(i, x̃) * x[i] / xe
    end

    return _∂ₓu
end

function ∂ₜu(mt::FullModel, mp::ModelParam, x::Vector{Float64}, dx::Vector{Float64}, l::Float64)
    _∂ₜu = 0
    xe = x[end]
    dxe = dx[end]
    x̃ = l / xe

    for i in 1:mp.N
        _∂ₜu += Φ(i, x̃) * dx[i] - x̃ * dxe / xe * ∂ₓΦ(i, x̃) * x[i]
    end

    return _∂ₜu
end

function ∂ₓₓu(mt::FullModel, mp::ModelParam, x::Vector{Float64}, l::Float64)
    xe = x[end]
    x̃ = l / xe

    _∂ₓₓu = 0
    for i in 1:mp.N
        _∂ₓₓu += ∂ₓₓΦ(i, x̃) * x[i] / xe^2
    end

    return _∂ₓₓu
end

function ∂ₓₜu(mt::FullModel, mp::ModelParam, x::Vector{Float64}, dx::Vector{Float64}, l::Float64)
    _∂ₓₜu = 0
    xe = x[end]
    dxe = dx[end]
    x̃ = l / xe

    for i in 1:mp.N
        _∂ₓₜu += xe * ∂ₓΦ(i, x̃) * dx[i]
        _∂ₓₜu += -dxe * ∂ₓΦ(i, x̃) * x[i]
        _∂ₓₜu += -x̃ * dxe * ∂ₓₓΦ(i, x̃) * x[i]
    end

    _∂ₓₜu = _∂ₓₜu / xe^2

    return _∂ₓₜu
end

function ∂ₓₑuₚ(mt::FullModel, mp::ModelParam, pulley_id::Integer, x::Vector{Float64})
    _∂ₓₑuₚ = 0
    xe = x[end]
    for i in 1:mp.N
        _∂ₓₑuₚ += x[i] * ∂ₓΦ(i, 1 - mp.lp[pulley_id] / xe)
    end
    _∂ₓₑuₚ *= (mp.lp[pulley_id] / xe^2)
    return _∂ₓₑuₚ
end

## θ function

function θ(mt::FullModel, mp::ModelParam, xe::Float64)
    _θ = (-xe + mp.L) / mp.rd
    return _θ
end

function θp(mt::FullModel, mp::ModelParam, pulley_id::Integer, x::Vector{Float64})
    xe = x[end]
    _θp = (mp.L - xe + u(mt, mp, x, xe - mp.lp[pulley_id])) / mp.rp[pulley_id]
    return _θp
end

function dθ(mt::FullModel, mp::ModelParam, dxe::Float64)
    _dθ = (-dxe) / mp.rd
    return _dθ
end

function dθp(mt::FullModel, mp::ModelParam, pulley_id::Integer, x::Vector{Float64}, dx::Vector{Float64})
    xe = x[end]
    dxe = dx[end]
    _∂ₜθ = (-dxe + ∂ₜu(mt, mp, x, dx, xe - mp.lp[pulley_id]) + ∂ₓu(mt, mp, x, xe - mp.lp[pulley_id]) * dxe) / mp.rp[pulley_id]
    return _∂ₜθ
end

function s(mt::FullModel, mp::ModelParam, x::Vector{Float64})
    xe = x[end]
    _θ = θ(mt, mp, xe)
    _θₚ = Vector{Float64}(undef, mp.pulley_num)
    for i in 1:mp.pulley_num
        _θₚ[i] = θp(mt, mp, i, x)
    end
    _s = [x; _θ; _θₚ]
    return _s
end

# construct ds = [dq; dθ; dθ₁, dθ₂,...,dθq]
function ds(mt::FullModel, mp::ModelParam, x::Vector{Float64}, dx::Vector{Float64})
    xe = x[end]
    dxe = dx[end]

    _dθ = dθ(mt, mp, dxe)
    _dθₚ = Vector{Float64}(undef, mp.pulley_num)
    for i in 1:mp.pulley_num
        _dθₚ[i] = dθp(mt, mp, i, x, dx)
    end
    _ds = [dx; _dθ; _dθₚ]
    return _ds
end

## ζ function
function ζ(mt::FullModel, mp::ModelParam, x::Vector{Float64}, dx::Vector{Float64}, x̃::Float64)
    _ζ = 0
    xe = x[end]
    dxe = dx[end]

    for i in 1:mp.N
        _ζ += -2 * dxe / xe * x̃ * ∂ₓΦ(i, x̃) * dx[i]
        _ζ += 2 * dxe^2 / xe^2 * x̃ * ∂ₓΦ(i, x̃) * x[i]
        _ζ += dxe^2 / xe^2 * x̃^2 * ∂ₓₓΦ(i, x̃) * x[i]
    end
    return _ζ
end

function ∫ζ(mt::FullModel, mp::ModelParam, x::Vector{Float64}, dx::Vector{Float64})
    _∫ζ = 0
    xe = x[end]
    dxe = dx[end]

    for i in 1:mp.N
        _∫ζ += -2 * dxe / xe * ∫x∂ₓΦ(i) * dx[i]
        _∫ζ += 2 * dxe^2 / xe^2 * ∫x∂ₓΦ(i) * x[i]
        _∫ζ += dxe^2 / xe^2 * ∫x²∂ₓₓΦ(i) * x[i]
    end
    return _∫ζ
end

function ∫ζΦ(mt::FullModel, mp::ModelParam, j::Integer, x::Vector{Float64}, dx::Vector{Float64})
    _∫ζΦ = 0
    xe = x[end]
    dxe = dx[end]

    for i in 1:mp.N
        _∫ζΦ += -2 * dxe / xe * ∫x∂ₓΦΦ(i, j) * dx[i]
        _∫ζΦ += 2 * dxe^2 / xe^2 * ∫x∂ₓΦΦ(i, j) * x[i]
        _∫ζΦ += dxe^2 / xe^2 * ∫x²∂ₓₓΦΦ(i, j) * x[i]
    end
    return _∫ζΦ
end