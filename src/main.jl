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
    include("hybrid_cable_model.jl")
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
    N = 10
    tspan = (0.0, 0.001)
    tol = 1e-10

    L = 5.0
    pulley_num = 1

    rd = 1.5e-2
    Id = 1.7485e-5
    Td = 0.155 * rd
    Cd = 0.020 * rd

    rp = 1.5e-2
    Ip = 1.7485e-5
    Tp = 0.0 * rd
    Cp = 0.0 * rd

    T = 100 * rd

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

# begin
#     # shows the cable deformation along the cable at t=0.2
#     # plot_deformation(tspan[2] * 1e-2, sol, mp, mt)
#     # plot_deformation(tspan[2] * 1e-4, sol, mp, mt)
#     # plot_deformation(tspan[2] * 1e-6, sol, mp, mt)
# end

begin
    sol_t_range = LinRange(tspan[1], tspan[2], 30)
    length_size = 1000
    ux = LinRange(0, x[end], length_size)

    function get_all_deformation(sol_t, sol, mp, mt, length_size=1000)
        sol_u_all = zeros(length_size)
        x = sol(sol_t)[1:mp.N+1]
        ux = LinRange(0, x[end], length_size)
        for i in 1:length_size
            sol_u_all[i] = u(mt, mp, x, ux[i]) * 1000 # unit: mm
        end
        return sol_u_all
    end

    fig = Figure()
    ax = Axis(fig[1, 1:2], xlabel="coordinate(m)", ylabel="deformation(mm)", title="t=" * string(0))

    record(fig, "all_deformation.mp4", sol_t_range; framerate=2) do t
        display(t)
        sol_u_all = get_all_deformation(t, sol, mp, mt)
        empty!(ax)
        lines!(ax, ux, sol_u_all, label="u_all", color=:tomato)
        ax.title = "t=" * string(t)
        ylims!(ax, -8e-2, 1e-5)
        display(fig)
    end

    record(fig, "all_deformation_autoaxis.mp4", sol_t_range; framerate=2) do t
        display(t)
        sol_u_all = get_all_deformation(t, sol, mp, mt)
        empty!(ax)
        lines!(ax, ux, sol_u_all, label="u_all", color=:tomato)
        ax.title = "t=" * string(t)
        autolimits!(ax)
        # ylims!(ax, -8e-2, 1e-5)
        display(fig)
    end
end


function plot_deformation(sol_t, sol, mp, mt, length_size=1000)
    sol_u_all = zeros(length_size)
    x = sol(sol_t)[1:mp.N+1]
    ux = LinRange(0, x[end], length_size)
    for i in 1:length_size
        sol_u_all[i] = u(mt, mp, x, ux[i]) * 1000 # unit: mm
    end

    fig4 = Figure()
    ax_dcu = Axis(fig4[1, 1:2], xlabel="coordinate(m)", ylabel="deformation(mm)", title="t=" * string(sol_t))
    lines!(ax_dcu, ux, sol_u_all, label="u_all")
    display(GLMakie.Screen(), fig4)
end