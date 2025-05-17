using .HybridCableModel

# =============== COMMON FUNCTIONS ====================

# extended model functions, used to accelerate the calculation of the model and avoid the discontinous stress condition
# i: [1, pulley_num]

function aᵢ(li::Float64, xe::Float64)
    ai = sqrt(3 / 2) * (xe) / (xe - li)
    return ai
end

function Φₑ(li::Float64, xe::Float64, x̃::Float64)
    ai = aᵢ(li, xe)
    if (x̃ <= 1 - li / xe)
        return ai * x̃
    else
        return ai * (xe / li - 1) * (1 - x̃)
    end
end

function ∂ₓΦₑ(li::Float64, xe::Float64, x̃::Float64)
    ai = aᵢ(li, xe)
    if (x̃ <= 1 - li / xe)
        return ai
    else
        return -ai * (xe / li - 1)
    end
end

function ∫Φₑ(li::Float64, xe::Float64)
    ai = aᵢ(li, xe)
    return (-li * ai + xe * ai) / (2 * xe)
end

# ∫ΦᵢΦₑ(lⱼ), obviously i ≠ j
function ∫ΦΦₑ(i::Integer, lj::Float64, xe::Float64)
    aj = aᵢ(lj, xe)

    return ((-1)^(1 + i) * xe * sin(i * pi * lj / xe) * aj) / (i^2 * pi^2 * lj)
end

# ∫Φₑ(lᵢ)Φₑ(lⱼ) where i ≠ j, otherwise result is 1/2
function ∫ΦₑΦₑ(li::Float64, lj::Float64, xe::Float64)
    ai = aᵢ(li, xe)
    aj = aᵢ(lj, xe)

    # li < lj case
    function ∫ΦₑΦₑ_tmp(li::Float64, lj::Float64, xe::Float64, ai::Float64, aj::Float64)
        return (xe - lj) * (-li^2 + (2 * xe - lj) * lj) * ai * aj / (6 * xe^2 * lj)
    end

    if li < lj
        return ∫ΦₑΦₑ_tmp(li, lj, xe, ai, aj)
    else
        return ∫ΦₑΦₑ_tmp(lj, li, xe, aj, ai)
    end
end

# ∫∂ₓΦᵢ∂ₓΦₑ(lⱼ), obviously i ≠ j
function ∫∂ₓΦ∂ₓΦₑ(i::Integer, lj::Float64, xe::Float64)
    aj = aᵢ(lj, xe)

    return ((-1)^(1 + i) * xe * sin(i * pi * lj / xe) * aj) / lj
end

# ∫∂ₓΦₑ(lᵢ)∂ₓΦₑ(lⱼ), if i==j this eqn still holds
function ∫∂ₓΦₑ∂ₓΦₑ(li::Float64, lj::Float64, xe::Float64)
    ai = aᵢ(li, xe)
    aj = aᵢ(lj, xe)

    # li < lj case
    function ∫∂ₓΦₑ∂ₓΦₑ_tmp(li::Float64, lj::Float64, xe::Float64, ai::Float64, aj::Float64)
        return (xe - lj) * ai * aj / lj
    end

    if li < lj
        return ∫∂ₓΦₑ∂ₓΦₑ_tmp(li, lj, xe, ai, aj)
    else
        return ∫∂ₓΦₑ∂ₓΦₑ_tmp(lj, li, xe, aj, ai)
    end
end

function ∫x∂ₓΦₑ(li::Float64, xe::Float64)
    ai = aᵢ(li, xe)
    return -(xe - li) * ai / (2 * xe)
end

# ∫x∂ₓΦᵢΦₑ(lj)
function ∫x∂ₓΦΦₑ(i::Integer, lj::Float64, xe::Float64)
    aj = aᵢ(lj, xe)
    return (-1)^i * (i * pi * (-1 + cos(i * pi * lj / xe)) * (xe - lj) + 2 * xe * sin(i * pi * lj / xe)) * aj / (i^2 * pi^2 * lj)
end

function ∫x∂ₓΦₑΦ(li::Float64, j::Integer, xe::Float64)
    ai = aᵢ(li, xe)
    return (-1)^(1 + j) * (j * pi * (-1 + cos(j * pi * li / xe)) * (xe - li) + xe * sin(j * pi * li / xe)) * ai / (j^2 * pi^2 * li)
end

function ∫x∂ₓΦₑΦₑ(li::Float64, lj::Float64, xe::Float64)
    ai = aᵢ(li, xe)
    aj = aᵢ(lj, xe)
    # li < lj case
    function ∫x∂ₓΦₑΦₑ_tmp(li::Float64, lj::Float64, xe::Float64, ai::Float64, aj::Float64)
        return (xe - lj) * (-3 * xe * li + 2 * li^2 + (2 * xe - lj) * lj) * ai * aj / (6 * xe^2 * lj)
    end

    if li < lj
        return ∫x∂ₓΦₑΦₑ_tmp(li, lj, xe, ai, aj)
    else
        return ∫x∂ₓΦₑΦₑ_tmp(lj, li, xe, aj, ai)
    end
end

function ∫x²∂ₓₓΦΦₑ(i::Integer, lj::Float64, xe::Float64)
    aj = aᵢ(lj, xe)
    denominator = (i^2 * pi^2 * xe * lj)
    part1 = -4 * i * pi * xe * (-1 + cos(i * pi * lj / xe)) * (xe - lj)
    part2 = ((-6 + i^2 * pi^2) * xe^2 + i^2 * pi^2 * lj * (-2 * xe + lj)) * sin(i * pi * lj / xe)
    return (-1)^i * (part1 + part2) * aj / denominator
end

function ∫∂ₓₓΦΦₑ(i::Integer, lj::Float64, xe::Float64)
    aj = aᵢ(lj, xe)
    return (-1)^i * xe * sin(i * pi * lj / xe) * aj / lj
end


# --------------------------------------------------------

function Φ(mp::ModelParam, xe::Float64, i::Integer, x̃::Float64)
    if i <= mp.N
        return sin(i * π * x̃)
    else
        return Φₑ(mp.lp[i-mp.N], xe, x̃)
    end
end

function ∂ₓΦ(mp::ModelParam, xe::Float64, i::Integer, x̃::Float64)
    if i <= mp.N
        return i * π * cos(i * π * x̃)
    else
        return ∂ₓΦₑ(mp.lp[i-mp.N], xe, x̃)
    end
end

function ∂ₓₓΦ(mp::ModelParam, i::Integer, x̃::Float64)
    if i <= mp.N
        return -(i * π)^2 * sin(i * π * x̃)
    else
        return 0
    end
end

function ∫Φ(mp::ModelParam, xe::Float64, i::Integer)
    if i <= mp.N
        return (1 + (-1)^(1 + i)) / (i * pi)
    else
        return ∫Φₑ(mp.lp[i-mp.N], xe)
    end
end

function ∫ΦΦ(mp::ModelParam, xe::Float64, i::Integer, j::Integer)
    # all trial function satisfy this eqn
    if i == j
        return 1 / 2
    end
    if (i <= mp.N)
        if (j <= mp.N)
            return 0
        end
        # i is Φ and j is Φₑ
        return ∫ΦΦₑ(i, mp.lp[j-mp.N], xe)
    end
    if j <= mp.N
        return ∫ΦΦₑ(j, mp.lp[i-mp.N], xe)
    end
    # i > N and j > N 
    return ∫ΦₑΦₑ(mp.lp[i-mp.N], mp.lp[j-mp.N], xe)
end

function ∫∂ₓΦ∂ₓΦ(mp::ModelParam, xe::Float64, i::Integer, j::Integer)
    if (i <= mp.N)
        if (j <= mp.N)
            if i == j
                return 1 / 2 * (i * pi)^2
            end
            return 0
        end
        return ∫∂ₓΦ∂ₓΦₑ(i, mp.lp[j-mp.N], xe)
    end
    if (j <= mp.N)
        return ∫∂ₓΦ∂ₓΦₑ(j, mp.lp[i-mp.N], xe)
    end
    # i > N and j > N
    return ∫∂ₓΦₑ∂ₓΦₑ(mp.lp[i-mp.N], mp.lp[j-mp.N], xe)
end

function ∫x∂ₓΦ(mp::ModelParam, xe::Float64, i::Integer)
    if (i <= mp.N)
        return (-1 + (-1)^(i)) / (i * pi)
    end
    return ∫x∂ₓΦₑ(mp.lp[i-mp.N], xe)
end

function ∫x²∂ₓₓΦ(mp::ModelParam, i::Integer)
    if (i <= mp.N)
        return (2 - 2 * (-1)^i + (-1)^i * i^2 * pi^2) / (i * pi)
    end
    return 0
end

function ∫x∂ₓΦΦ(mp::ModelParam, xe::Float64, i::Integer, j::Integer)
    if (i <= mp.N)
        if (j <= mp.N)
            if i == j
                return -1 / 4
            end
            return ((-1)^(i + j) * i * j) / (i^2 - j^2)
        end
        return ∫x∂ₓΦΦₑ(i, mp.lp[j-mp.N], xe)
    end
    if (j <= mp.N)
        return ∫x∂ₓΦₑΦ(mp.lp[i-mp.N], j, xe)
    end
    # i > N and j > N
    return ∫x∂ₓΦₑΦₑ(mp.lp[i-mp.N], mp.lp[j-mp.N], xe)
end

function ∫x²∂ₓₓΦΦ(mp::ModelParam, xe::Float64, i::Integer, j::Integer)
    if (i <= mp.N)
        if (j <= mp.N)
            if i == j
                return 1 / 4 - i^2 * pi^2 / 6
            end
            return -(4 * (-1)^(i + j) * i^3 * j) / ((i^2 - j^2)^2)
        end
        return ∫x²∂ₓₓΦΦₑ(i, mp.lp[j-mp.N], xe)
    end
    return 0
end

function ∫∂ₓₓΦΦ(mp::ModelParam, xe::Float64, i::Integer, j::Integer)
    if (i <= mp.N)
        if (j <= mp.N)
            if i == j
                return -1 / 2 * (i * pi)^2
            end
            return 0
        end
        return ∫∂ₓₓΦΦₑ(i, mp.lp[j-mp.N], xe)
    end
    return 0
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

function get_cable_force(mt::ModelType, mp::ModelParam, x::Vector{Float64}, l::Float64)
    uₓ = ∂ₓu(mt, mp, x, l)
    return uₓ * mp.E * mp.A
end

function get_cable_force(mp::ModelParam, xe::Float64, ddxe::Float64)
    cable_force = mp.k * (mp.L - xe) - mp.m * ddxe
    return cable_force
end

# =============== FULL MODEL ====================

## u function

# u(l,t) = ∑ Φᵢ(l) xᵢ
function u(mt::FullModel, mp::ModelParam, x::Vector{Float64}, l::Float64)
    xe = x[end]
    _u = 0
    for i = 1:(mp.N+mp.pulley_num)
        _u += Φ(mp, xe, i, l / xe) * x[i]
    end
    return _u
end

function ∂ₓu(mt::FullModel, mp::ModelParam, x::Vector{Float64}, l::Float64)
    xe = x[end]
    _∂ₓu = 0
    x̃ = l / xe

    for i = 1:(mp.N+mp.pulley_num)
        _∂ₓu += ∂ₓΦ(mp, xe, i, x̃) * x[i] / xe
    end

    return _∂ₓu
end

function ∂ₜu(mt::FullModel, mp::ModelParam, x::Vector{Float64}, dx::Vector{Float64}, l::Float64)
    _∂ₜu = 0
    xe = x[end]
    dxe = dx[end]
    x̃ = l / xe

    for i in 1:(mp.N+mp.pulley_num)
        _∂ₜu += Φ(mp, xe, i, x̃) * dx[i] - x̃ * dxe / xe * ∂ₓΦ(mp, xe, i, x̃) * x[i]
    end

    return _∂ₜu
end

function ∂ₓₓu(mt::FullModel, mp::ModelParam, x::Vector{Float64}, l::Float64)
    xe = x[end]
    x̃ = l / xe

    _∂ₓₓu = 0
    for i in 1:(mp.N+mp.pulley_num)
        _∂ₓₓu += ∂ₓₓΦ(mp, i, x̃) * x[i] / xe^2
    end

    return _∂ₓₓu
end

function ∂ₓₜu(mt::FullModel, mp::ModelParam, x::Vector{Float64}, dx::Vector{Float64}, l::Float64)
    _∂ₓₜu = 0
    xe = x[end]
    dxe = dx[end]
    x̃ = l / xe

    for i in 1:(mp.N+mp.pulley_num)
        _∂ₓₜu += xe * ∂ₓΦ(mp, xe, i, x̃) * dx[i]
        _∂ₓₜu += -dxe * ∂ₓΦ(mp, xe, i, x̃) * x[i]
        _∂ₓₜu += -x̃ * dxe * ∂ₓₓΦ(mp, i, x̃) * x[i]
    end

    _∂ₓₜu = _∂ₓₜu / xe^2

    return _∂ₓₜu
end

function ∂ₓₑuₚ(mt::FullModel, mp::ModelParam, pulley_id::Integer, x::Vector{Float64})
    _∂ₓₑuₚ = 0
    xe = x[end]
    for i in 1:(mp.N+mp.pulley_num)
        _∂ₓₑuₚ += x[i] * ∂ₓΦ(mp, xe, i, 1 - mp.lp[pulley_id] / xe)
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

    for i in 1:(mp.N+mp.pulley_num)
        _ζ += -2 * dxe / xe * x̃ * ∂ₓΦ(mp, xe, i, x̃) * dx[i]
        _ζ += 2 * dxe^2 / xe^2 * x̃ * ∂ₓΦ(mp, xe, i, x̃) * x[i]
        _ζ += dxe^2 / xe^2 * x̃^2 * ∂ₓₓΦ(mp, i, x̃) * x[i]
    end
    return _ζ
end

function ∫ζ(mt::FullModel, mp::ModelParam, x::Vector{Float64}, dx::Vector{Float64})
    _∫ζ = 0
    xe = x[end]
    dxe = dx[end]

    for i in 1:(mp.N+mp.pulley_num)
        _∫ζ += -2 * dxe / xe * ∫x∂ₓΦ(mp, xe, i) * dx[i]
        _∫ζ += 2 * dxe^2 / xe^2 * ∫x∂ₓΦ(mp, xe, i) * x[i]
        _∫ζ += dxe^2 / xe^2 * ∫x²∂ₓₓΦ(mp, i) * x[i]
    end
    return _∫ζ
end

function ∫ζΦ(mt::FullModel, mp::ModelParam, j::Integer, x::Vector{Float64}, dx::Vector{Float64})
    _∫ζΦ = 0
    xe = x[end]
    dxe = dx[end]

    for i in 1:(mp.N+mp.pulley_num)
        _∫ζΦ += -2 * dxe / xe * ∫x∂ₓΦΦ(mp, xe, i, j) * dx[i]
        _∫ζΦ += 2 * dxe^2 / xe^2 * ∫x∂ₓΦΦ(mp, xe, i, j) * x[i]
        _∫ζΦ += dxe^2 / xe^2 * ∫x²∂ₓₓΦΦ(mp, xe, i, j) * x[i]
    end
    return _∫ζΦ
end

# =============== PARTIAL MODEL ====================

## u function

# u(l,t) = ∑ Φᵢ(l) xᵢ
function u(mt::PartialModel, mp::ModelParam, x::Vector{Float64}, l::Float64)
    xe = x[end]
    _u = 0
    for i = 1:(mp.N+mp.pulley_num)
        _u += Φ(mp, xe, i, l / xe) * x[i]
    end
    return _u
end

function ∂ₓu(mt::PartialModel, mp::ModelParam, x::Vector{Float64}, l::Float64)
    xe = x[end]
    _∂ₓu = 0
    x̃ = l / xe

    for i = 1:(mp.N+mp.pulley_num)
        _∂ₓu += ∂ₓΦ(mp, xe, i, x̃) * x[i] / xe
    end

    return _∂ₓu
end

function ∂ₜu(mt::PartialModel, mp::ModelParam, x::Vector{Float64}, dx::Vector{Float64}, l::Float64)
    _∂ₜu = 0
    xe = x[end]
    dxe = dx[end]
    x̃ = l / xe

    for i in 1:(mp.N+mp.pulley_num)
        _∂ₜu += Φ(mp, xe, i, x̃) * dx[i]
    end

    return _∂ₜu
end

function ∂ₓₓu(mt::PartialModel, mp::ModelParam, x::Vector{Float64}, l::Float64)
    xe = x[end]
    x̃ = l / xe

    _∂ₓₓu = 0
    for i in 1:(mp.N+mp.pulley_num)
        _∂ₓₓu += ∂ₓₓΦ(mp, i, x̃) * x[i] / xe^2
    end

    return _∂ₓₓu
end

function ∂ₓₜu(mt::PartialModel, mp::ModelParam, x::Vector{Float64}, dx::Vector{Float64}, l::Float64)
    _∂ₓₜu = 0
    xe = x[end]
    dxe = dx[end]
    x̃ = l / xe

    for i in 1:(mp.N+mp.pulley_num)
        _∂ₓₜu += ∂ₓΦ(mp, xe, i, x̃) * dx[i] / xe
    end

    return _∂ₓₜu
end

function ∂ₓₑuₚ(mt::PartialModel, mp::ModelParam, pulley_id::Integer, x::Vector{Float64})
    _∂ₓₑuₚ = 0
    xe = x[end]
    for i in 1:(mp.N+mp.pulley_num)
        _∂ₓₑuₚ += x[i] * ∂ₓΦ(mp, xe, i, 1 - mp.lp[pulley_id] / xe)
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

# =============== SIMPLE MODEL ====================

## u function

# u(l,t) = ∑ Φᵢ(l) xᵢ
function u(mt::SimpleModel, mp::ModelParam, x::Vector{Float64}, l::Float64)
    xe = x[end]
    _u = 0
    for i = 1:(mp.N)
        _u += Φ(mp, xe, i, l / xe) * x[i]
    end
    return _u
end

function ∂ₓu(mt::SimpleModel, mp::ModelParam, x::Vector{Float64}, l::Float64)
    xe = x[end]
    _∂ₓu = 0
    x̃ = l / xe

    for i = 1:(mp.N)
        _∂ₓu += ∂ₓΦ(mp, xe, i, x̃) * x[i] / xe
    end

    return _∂ₓu
end

function dθ(mt::SimpleModel, mp::ModelParam, dxe::Float64)
    _dθ = (-dxe) / mp.rd
    return _dθ
end