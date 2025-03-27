module HybridCableModel
export ModelType, FullModel, SimpleModel, PartialModel, SolverParam, ModelParam, HybridCableODE

abstract type ModelType end
struct FullModel <: ModelType end
struct PartialModel <: ModelType end
struct SimpleModel <: ModelType end

@kwdef struct ModelParam
    E = 78.0e9
    A = π * (1.2e-3 / 2)^2
    ρ = 7.5e-3

    L = 5.0

    # winch params
    rd = 1.5e-2
    ld = 0.0  # residual length to xe
    Id = 1.7485e-5

    # pulley params
    Ip = [1.7485e-5;]
    lp = [0.42;]    # residual length to xe
    rp = [rd;]
    pulley_num = 1

    # mass params
    m = 0.65
    #  k = 100.0
    k = 100

    # friction params
    Td = 0.0 * rd
    Cd = 0.0 * rd
    Tp = [Td;]
    Cp = [Cd;]
    Tm = 0.0
    Cm = 0.0

    T = 48.37 * rd
    N = 3
    tspan = (0.0, 1.0)
    tol = 1e-10
end
end