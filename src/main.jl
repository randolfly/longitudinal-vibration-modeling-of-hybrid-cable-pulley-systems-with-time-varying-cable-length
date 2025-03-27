begin
    using DifferentialEquations
    using GLMakie
    using LinearAlgebra
    # using DelimitedFiles
    using ForwardDiff
    using DataFrames
    using PrettyTables
    using BenchmarkTools, Profile
    # using MAT
    using .HybridCableModel
end
begin
    include("hybrid_cable_model.jl")
    include("util.jl")
    include("ode.jl")
    include("bc.jl")
    include("post.jl")
end

begin
    N = 3
    tspan = (0.0, 0.3)
    tol = 1e-10

    L = 5.0
    pulley_num = 1
    T = 0.75

    rd = 1.5e-2
    Id = 1.7485e-5 + 1.2e-3
    Td = 2.0 * rd
    Cd = 0.0 * rd

    rp = 1.5e-2
    Ip = 1.7485e-5
    Tp = 1.0 * rd
    Cp = 0.0 * rd

    k = 100.0
    m = 0.65

    # mt = FullModel()
    # mt = PartialModel()
    mt = SimpleModel()

    mp = ModelParam(
        N=N,
        tspan=tspan,
        tol=tol,
        L=L,
        rd=rd,
        Id=Id,
        k=k,
        T=T,
        Td=Td,
        Cd=Cd,
        pulley_num=pulley_num,
        rp=fill(rp, pulley_num),
        lp=[0.42 + 1.04 * (i - 1) for i = 1:pulley_num],
        Ip=fill(Ip, pulley_num),
        Tp=fill(Tp, pulley_num),
        Cp=fill(Cp, pulley_num)
    )
    display(mt)
    show_model_param(mp)
end

begin
    X0, dX0 = load_init_conditions(mp)
    prob = ODEProblem(ode_function!, X0, mp.tspan, mp)
    @time sol = solve(prob, TRBDF2(autodiff=false), reltol=mp.tol, abstol=mp.tol)
    display("solver finished!")
    # @benchmark sol = solve(prob, TRBDF2(autodiff=false), reltol=mp.tol, abstol=mp.tol) samples = 40 seconds = 80
end

post_sol(mp, sol)