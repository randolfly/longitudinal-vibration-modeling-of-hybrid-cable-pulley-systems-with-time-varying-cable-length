using .HybridCableModel

function ode_function!(dX, X, mp, t)
    x = X[1:(mp.N+1)]
    dx = X[(mp.N+2):end]

    _s = s(mt, mp, x)
    _ds = ds(mt, mp, x, dx)
    _M = M(mt, mp, _s, _ds)
    _F = F(mt, mp, _s, _ds)
    _G = G(mt, mp, x, dx)
    _C = C(mt, mp, x, dx)

    id_M = I(mp.N + 1)
    zero_c = zeros(mp.N + 1, 1)

    _A = _M * [id_M; _G]
    _B = _F - _M * [zero_c; _C]
    ddx = inv(_A) * _B
    dX[1:mp.N+1] = dx
    dX[mp.N+2:end] = ddx
    nothing
end

function M(mt::PartialModel, mp::ModelParam, s::Vector{Float64}, ds::Vector{Float64})
    # s = [q; θ; θ₁, θ₂,...,θq]
    x = s[1:mp.N+1]
    xe = x[end]
    θ = s[mp.N+2]
    θₚ = s[mp.N+3:end]

    dx = ds[1:mp.N+1]
    dxe = dx[end]
    dθ = ds[mp.N+2]
    dθₚ = ds[mp.N+3:end]

    M11 = (-mp.ρ * xe) .* [∫Φ(i) for i in 1:mp.N]'
    M12 = mp.m + mp.ρ * xe
    M13 = -(mp.Id + mp.ρ * mp.rd^3 * θ) / mp.rd
    M14 = zeros(1, mp.pulley_num)
    for i in 1:mp.pulley_num
        M14[i] = mp.Ip[i] / mp.rp[i] * (∂ₓₑuₚ(mt, mp, i, x) - 1)
    end

    M21 = zeros(mp.N, mp.N)
    for i in 1:mp.N
        for j in 1:mp.N
            M21[i, j] = ∫ΦΦ(i, j)
        end
    end

    M22 = zeros(mp.N, 1)
    M22_2 = zeros(mp.N, 1)
    for i in 1:mp.N
        M22_2[i] = -∫Φ(i)
    end
    M22 = M22_2

    M23 = zeros(mp.N, 1)

    M24 = zeros(mp.N, mp.pulley_num)
    for i in 1:mp.N
        for j in 1:mp.pulley_num
            M24[i, j] = mp.Ip[j] / (mp.ρ * mp.rp[j]) * Φ(i, 1 - mp.lp[j] / xe)
        end
    end

    _M = [M11 M12 M13 M14;
        M21 M22 M23 M24
    ]
    return _M
end

function F(mt::PartialModel, mp::ModelParam, s::Vector{Float64}, ds::Vector{Float64})
    # s = [q; θ; θ₁, θ₂,...,θq]
    x = s[1:mp.N+1]
    xe = x[end]
    θ = s[mp.N+2]
    θₚ = s[mp.N+3:end]

    dx = ds[1:mp.N+1]
    dxe = dx[end]
    dθ = ds[mp.N+2]
    dθₚ = ds[mp.N+3:end]

    a = sqrt(mp.A * mp.E / mp.ρ)

    Fnon = mp.ρ * dxe^2 + mp.k * (xe - mp.L) - fm(dxe, mp.Tm, mp.Cm)
    Fnon += -1 / mp.rd * (1 / 2 * mp.ρ * mp.rd^3 * dθ^2 + f(dθ, mp.Td, mp.Cd) - mp.T)
    for i in 1:mp.pulley_num
        Fnon += (∂ₓₑuₚ(mt, mp, i, x) - 1) / mp.rp[i] * f(dθₚ[i], mp.Tp[i], mp.Cp[i])
    end

    Fs1 = zeros(mp.N + 1, 1)
    Fs1[1] = -Fnon
    for i in 1:mp.N
        _tmp = 0
        for j in 1:mp.N
            _tmp += ∫∂ₓΦ∂ₓΦ(i, j) * x[j]
        end
        Fs1[i+1] = -a^2 / xe^2 * _tmp
    end

    Fs2 = zeros(mp.N + 1, 1)
    for i in 1:mp.N
        _tmp = 0
        for j in 1:mp.pulley_num
            _tmp += f(dθₚ[j], mp.Tp[j], mp.Cp[j]) / (mp.ρ * mp.rp[j]) * Φ(i, 1 - mp.lp[j] / xe)
        end
        Fs2[i+1] = -_tmp
    end

    _F = Fs1 + Fs2
    return _F
end

function G(mt::PartialModel, mp::ModelParam, x::Vector{Float64}, dx::Vector{Float64})
    xe = x[end]
    dxe = dx[end]

    _G = zeros(mp.pulley_num + 1, mp.N + 1)
    for i in 1:mp.N
        _G[1, i] = 0
    end
    _G[1, mp.N+1] = -1 / mp.rd

    for i in 1:mp.pulley_num
        for j in 1:mp.N
            _G[i+1, j] = Φ(j, 1 - mp.lp[i] / xe) / mp.rp[i]
        end
        _G[i+1, mp.N+1] = (∂ₓu(mt, mp, x, xe - mp.lp[i]) - 1) / mp.rp[i]
    end

    return _G
end

function C(mt::PartialModel, mp::ModelParam, x::Vector{Float64}, dx::Vector{Float64})
    xe = x[end]
    dxe = dx[end]

    _C = zeros(mp.pulley_num + 1, 1)
    _C[1] = 0
    for i in 1:mp.pulley_num
        _C[i+1] = 2 / mp.rp[i] * ∂ₓₜu(mt, mp, x, dx, xe - mp.lp[i]) * dxe
        _C[i+1] += 1 / mp.rp[i] * ∂ₓₓu(mt, mp, x, xe - mp.lp[i]) * dxe^2
    end

    return _C
end


function M(mt::FullModel, mp::ModelParam, s::Vector{Float64}, ds::Vector{Float64})
    # s = [q; θ; θ₁, θ₂,...,θq]
    x = s[1:mp.N+1]
    xe = x[end]
    θ = s[mp.N+2]
    θₚ = s[mp.N+3:end]

    dx = ds[1:mp.N+1]
    dxe = dx[end]
    dθ = ds[mp.N+2]
    dθₚ = ds[mp.N+3:end]

    M11 = (-mp.ρ * xe) .* [∫Φ(i) for i in 1:mp.N]'
    M12 = mp.m + mp.ρ * xe
    for i in 1:mp.N
        M12 += mp.ρ * x[i] * ∫x∂ₓΦ(i)
    end
    M13 = -(mp.Id + mp.ρ * mp.rd^3 * θ) / mp.rd
    M14 = zeros(1, mp.pulley_num)
    for i in 1:mp.pulley_num
        M14[i] = mp.Ip[i] / mp.rp[i] * (∂ₓₑuₚ(mt, mp, i, x) - 1)
    end

    M21 = zeros(mp.N, mp.N)
    for i in 1:mp.N
        for j in 1:mp.N
            M21[i, j] = ∫ΦΦ(i, j)
        end
    end

    M22 = zeros(mp.N, 1)
    M22_1 = zeros(mp.N, 1)
    M22_2 = zeros(mp.N, 1)
    for i in 1:mp.N
        M22_1[i] = 0
        for j in 1:mp.N
            M22_1[i] += ∫x∂ₓΦΦ(j, i) * x[j]
        end
    end
    M22_1 = -(1 / xe) * M22_1

    for i in 1:mp.N
        M22_2[i] = -∫Φ(i)
    end

    M22 = M22_1 + M22_2

    M23 = zeros(mp.N, 1)

    M24 = zeros(mp.N, mp.pulley_num)
    for i in 1:mp.N
        for j in 1:mp.pulley_num
            M24[i, j] = mp.Ip[j] / (mp.ρ * mp.rp[j]) * Φ(i, 1 - mp.lp[j] / xe)
        end
    end

    _M = [M11 M12 M13 M14;
        M21 M22 M23 M24
    ]
    return _M
end

function F(mt::FullModel, mp::ModelParam, s::Vector{Float64}, ds::Vector{Float64})
    # s = [q; θ; θ₁, θ₂,...,θq]
    x = s[1:mp.N+1]
    xe = x[end]
    θ = s[mp.N+2]
    θₚ = s[mp.N+3:end]

    dx = ds[1:mp.N+1]
    dxe = dx[end]
    dθ = ds[mp.N+2]
    dθₚ = ds[mp.N+3:end]
    a = sqrt(mp.A * mp.E / mp.ρ)

    Fnon = mp.ρ * dxe^2 + mp.k * (xe - mp.L) - fm(dxe, mp.Tm, mp.Cm)
    Fnon += -1 / mp.rd * (1 / 2 * mp.ρ * mp.rd^3 * dθ^2 + f(dθ, mp.Td, mp.Cd) - mp.T)
    for i in 1:mp.pulley_num
        Fnon += (∂ₓₑuₚ(mt, mp, i, x) - 1) / mp.rp[i] * f(dθₚ[i], mp.Tp[i], mp.Cp[i])
    end

    Fs1 = zeros(mp.N + 1, 1)
    Fs1[1] = -Fnon + mp.ρ * xe * ∫ζ(mt, mp, x, dx)
    for i in 1:mp.N
        _tmp = 0
        for j in 1:mp.N
            _tmp += ∫∂ₓΦ∂ₓΦ(i, j) * x[j]
        end
        Fs1[i+1] = -a^2 / xe^2 * _tmp
    end

    Fs2 = zeros(mp.N + 1, 1)
    for i in 1:mp.N
        _tmp = 0
        for j in 1:mp.pulley_num
            _tmp += f(dθₚ[j], mp.Tp[j], mp.Cp[j]) / (mp.ρ * mp.rp[j]) * Φ(i, 1 - mp.lp[j] / xe)
        end
        Fs2[i+1] = -(_tmp + ∫ζΦ(mt, mp, i, x, dx))
    end

    _F = Fs1 + Fs2
    return _F
end

function G(mt::FullModel, mp::ModelParam, x::Vector{Float64}, dx::Vector{Float64})
    xe = x[end]
    dxe = dx[end]

    _G = zeros(mp.pulley_num + 1, mp.N + 1)
    for i in 1:mp.N
        _G[1, i] = 0
    end
    _G[1, mp.N+1] = -1 / mp.rd

    for i in 1:mp.pulley_num
        for j in 1:mp.N
            _G[i+1, j] = Φ(j, 1 - mp.lp[i] / xe) / mp.rp[i]
        end
        _G[i+1, mp.N+1] = (∂ₓu(mt, mp, x, xe - mp.lp[i]) - 1) / mp.rp[i]
        _tmp = 0
        for j in 1:mp.N
            _tmp += ∂ₓΦ(j, 1 - mp.lp[i] / xe) * x[j]
        end
        _G[i+1, mp.N+1] += -(1 - mp.lp[i] / xe) / (mp.rp[i] * xe) * _tmp
    end

    return _G
end

function C(mt::FullModel, mp::ModelParam, x::Vector{Float64}, dx::Vector{Float64})
    xe = x[end]
    dxe = dx[end]

    _C = zeros(mp.pulley_num + 1, 1)
    _C[1] = 0
    for i in 1:mp.pulley_num
        _C[i+1] = 2 / mp.rp[i] * ∂ₓₜu(mt, mp, x, dx, xe - mp.lp[i]) * dxe
        _C[i+1] += 1 / mp.rp[i] * ∂ₓₓu(mt, mp, x, xe - mp.lp[i]) * dxe^2
        _C[i+1] += 1 / mp.rp[i] * ζ(mt, mp, x, dx, 1 - mp.lp[i] / xe)
    end

    return _C
end
