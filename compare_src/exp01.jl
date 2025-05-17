# EXPERIMENT 01: comparision of the extended trial space results and the normal test space results
# COMPARISION MODELS:
# 1. FullModel: full model with all the dynamic parts, discretization N is N=3,10,50
# EXPERIMENTS RESULTS: the extended trial space results are better than the normal test space results, but the calculation time is longer.

begin
    using DifferentialEquations
    using Sundials
    using GLMakie
    using LinearAlgebra
    # using DelimitedFiles
    using ForwardDiff
    using DataFrames
    using PrettyTables
    using BenchmarkTools, Profile
    using MAT
    include("hybrid_cable_model.jl")
    import .HybridCableModel
end
begin
    include("hybrid_cable_model.jl")
    include("util.jl")
    include("ode.jl")
    include("bc.jl")
    include("post.jl")
end

begin
    N = 50
    tspan = (0.0, 0.04)
    tol = 1e-8

    L = 5.0
    pulley_num = 3
    rd = 1.5e-2
    Id = 1.7485e-5
    Td = 0.0 * rd
    Cd = 0.0 * rd

    rp = 1.5e-2
    Ip = 1.7485e-5
    Tp = Td
    Cp = Cd

    T = 100.0 * rd

    k = 300.0
    m = 0.65

    mt = FullModel()
    # mt = PartialModel()

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
    # io = open("solver_stat.txt", "w")
    X0, dX0 = load_init_conditions(mp)
    prob = ODEProblem(ode_function!, X0, mp.tspan, mp)

    # @time sol = solve(prob, TRBDF2(autodiff=false), reltol=mp.tol, abstol=mp.tol)

    display("solver finished!")
    # close(io)
    @benchmark sol = solve(prob, TRBDF2(autodiff=false), reltol=mp.tol, abstol=mp.tol) samples = 30 seconds = 1000
end

# title = "Exp01_Full_P3_N3"
# title = "Exp01_Full_P3_N10"
# title = "Exp01_Full_P3_N50"
# title = "Exp01_Full_P3_N100"

# title = "Exp01_Full_P1_N3"
# title = "Exp01_Full_P1_N10"
# title = "Exp01_Full_P1_N50"
# title = "Exp01_Full_P1_N100"

# post_sol(mt, mp, sol, true, title)

# generate 3D deformation data
begin
    frame_rate = 24
    sol_t_range = LinRange(sol.t[1], sol.t[end], frame_rate * 5)
    length_size = 1000

    function get_all_deformation(sol_t, sol, mp, mt, length_size=1000)
        sol_u_all = zeros(length_size)
        x = sol(sol_t)[1:mp.N+1]
        ux = LinRange(0, x[end], length_size)
        for i in 1:length_size
            sol_u_all[i] = u(mt, mp, x, ux[i]) * 1000 # unit: mm
        end
        return sol_u_all
    end

    function get_all_strains(sol_t, sol, mp, mt, length_size=1000)
        sol_uₓ_all = zeros(length_size)
        x = sol(sol_t)[1:mp.N+1]
        ux = LinRange(0, x[end], length_size)
        for i in 1:length_size
            sol_uₓ_all[i] = ∂ₓu(mt, mp, x, ux[i]) * 1000 # unit: mm
        end
        return sol_uₓ_all
    end

    fig = Figure()
    ax1 = Axis(fig[1, 1:2], xlabel="coordinate(x/xe)", ylabel="deformation(mm)", title="t=" * string(0))
    ax2 = Axis(fig[2, 1:2], xlabel="coordinate(x/xe)", ylabel="strain")
    ux = LinRange(0, 1, length_size)
    # record(fig, "ideal_1pulley_deformation_with_fixlimit_strain.mp4", sol_t_range; framerate=frame_rate) do t
    #     display(t)
    #     sol_u_all = get_all_deformation(t, sol, mp, mt)
    #     sol_uₓ_all = get_all_strains(t, sol, mp, mt)
    #     empty!(ax1)
    #     empty!(ax2)

    #     lines!(ax1, ux, sol_u_all, label="u_all", color=:tomato)
    #     lines!(ax2, ux, sol_uₓ_all, label="u_all", color=:tomato)

    #     ax1.title = "t=" * string(t)

    #     ylims!(ax1, -4e-1, 1e-5)

    #     display(fig)
    # end

    record(fig, title * "_deformation.mp4", sol_t_range; framerate=frame_rate) do t
        display(t)
        sol_u_all = get_all_deformation(t, sol, mp, mt)
        sol_uₓ_all = get_all_strains(t, sol, mp, mt)
        empty!(ax1)
        empty!(ax2)

        lines!(ax1, ux, sol_u_all, label="u_all", color=:tomato)
        lines!(ax2, ux, sol_uₓ_all, label="u_all", color=:tomato)
        autolimits!(ax1)
        autolimits!(ax2)
        ax1.title = "t=" * string(t)

        display(fig)
    end
end