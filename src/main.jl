begin
    using DifferentialEquations
    using GLMakie
    using LinearAlgebra
    # using DelimitedFiles
    using ForwardDiff
    using PrettyTables
    # using BenchmarkTools, Profile
    # using MAT
end

begin
    include("hybrid_cable_model.jl")
    using .HybridCableModel
    include("util.jl")
    include("ode.jl")
    include("bc.jl")
    include("post.jl")

    function solve_hybrid_ode(mp::ModelParam)
        X0, dX0 = load_init_conditions(mp)
        prob = ODEProblem(ode_function!, X0, mp.tspan, mp)
        @time sol = solve(prob, TRBDF2(autodiff=false), reltol=mp.tol, abstol=mp.tol)
        return sol
    end
end

begin
    N = 3
    tspan = (0.0, 1.0)
    tol = 1e-10

    pulley_num = 3
    rd = 1.5e-2
    T = 48.37 * rd
    Td = 0.155 * rd
    Cd = 0.020 * rd

    mt = FullModel()
    # mt = PartialModel()
    mp = ModelParam(
        N=N,
        tspan=tspan,
        tol=tol,
        L=0.42 + 1.04 * pulley_num,
        rd=rd,
        T=T,
        Ip=fill(1.7485e-5, pulley_num),
        lp=[0.42 + 1.04 * (i - 1) for i = 1:pulley_num],
        rp=fill(rd, pulley_num),
        pulley_num=pulley_num,
        Tp=fill(Td, pulley_num),
        Cp=fill(Cd, pulley_num)
    )
    show_model_param(mp)
end

sol = solve_hybrid_ode(mp);
post_sol(mp, sol);