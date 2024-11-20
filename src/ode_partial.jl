function M(s::Vector{Float64}, ds::Vector{Float64})
    # s = [q; θ; θ₁, θ₂,...,θq]
    x = s[1:N+1]
    xe = x[end]
    θ = s[N+2]
    θₚ = s[N+3:end]

    dx = ds[1:N+1]
    dxe = dx[end]
    dθ = ds[N+2]
    dθₚ = ds[N+3:end]

    M11 = (-ρ * xe) .* [∫Φ(i) for i in 1:N]'
    M12 = m + ρ * xe
    M13 = -(Id + ρ * rd^3 * θ) / rd
    M14 = zeros(1, pulley_num)
    for i in 1:pulley_num
        M14[i] = Ip[i] / rp[i] * (∂ₓₑuₚ(i, x) - 1)
    end

    M21 = zeros(N, N)
    for i in 1:N
        for j in 1:N
            M21[i, j] = ∫ΦΦ(i, j)
        end
    end

    M22 = zeros(N, 1)
    M22_2 = zeros(N, 1)
    for i in 1:N
        M22_2[i] = -∫Φ(i)
    end
    M22 = M22_2

    M23 = zeros(N, 1)

    M24 = zeros(N, pulley_num)
    for i in 1:N
        for j in 1:pulley_num
            M24[i, j] = Ip[j] / (ρ * rp[j]) * Φ(i, 1 - lp[j] / xe)
        end
    end

    _M = [M11 M12 M13 M14;
        M21 M22 M23 M24
    ]
    return _M
end

function F(s::Vector{Float64}, ds::Vector{Float64})
    # s = [q; θ; θ₁, θ₂,...,θq]
    x = s[1:N+1]
    xe = x[end]
    θ = s[N+2]
    θₚ = s[N+3:end]

    dx = ds[1:N+1]
    dxe = dx[end]
    dθ = ds[N+2]
    dθₚ = ds[N+3:end]

    Fnon = ρ * dxe^2 + k * (xe - L) - fm(dxe)
    Fnon += -1 / rd * (1 / 2 * ρ * rd^3 * dθ^2 + f(dθ, Td, Cd) - T())
    for i in 1:pulley_num
        Fnon += (∂ₓₑuₚ(i, x) - 1) / rp[i] * f(dθₚ[i], Tp[i], Cp[i])
    end

    Fs1 = zeros(N + 1, 1)
    Fs1[1] = -Fnon
    for i in 1:N
        _tmp = 0
        for j in 1:N
            _tmp += ∫∂ₓΦ∂ₓΦ(i, j) * x[j]
        end
        Fs1[i+1] = -a^2 / xe^2 * _tmp
    end

    Fs2 = zeros(N + 1, 1)
    for i in 1:N
        _tmp = 0
        for j in 1:pulley_num
            _tmp += f(dθₚ[j], Tp[j], Cp[j]) / (ρ * rp[j]) * Φ(i, 1 - lp[j] / xe)
        end
        Fs2[i+1] = -_tmp
    end

    _F = Fs1 + Fs2
    return _F
end

function G(x::Vector{Float64}, dx::Vector{Float64})
    xe = x[end]
    dxe = dx[end]

    _G = zeros(pulley_num + 1, N + 1)
    for i in 1:N
        _G[1, i] = 0
    end
    _G[1, N+1] = -1 / rd

    for i in 1:pulley_num
        for j in 1:N
            _G[i+1, j] = Φ(j, 1 - lp[i] / xe) / rp[i]
        end
        _G[i+1, N+1] = (∂ₓu(x, xe - lp[i]) - 1) / rp[i]
    end

    return _G
end

function C(x::Vector{Float64}, dx::Vector{Float64})
    xe = x[end]
    dxe = dx[end]

    _C = zeros(pulley_num + 1, 1)
    _C[1] = 0
    for i in 1:pulley_num
        _C[i+1] = 2 / rp[i] * ∂ₓₜu(x, dx, xe - lp[i]) * dxe
        _C[i+1] += 1 / rp[i] * ∂ₓₓu(x, xe - lp[i]) * dxe^2
    end

    return _C
end
