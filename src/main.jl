## init program

# include the necessary packages
begin
    using DifferentialEquations
    using GLMakie
    using LinearAlgebra
    using DelimitedFiles
    using ForwardDiff
    using BenchmarkTools, Profile
end

model_type = "full"
# model_type = "partial"
# param_type = "ideal_1p"
param_type = "ideal"
# param_type = "3p"
# param_type = "1p"

# init config params
begin
    include("param_" * param_type * ".jl")
end

# simulation params
begin
    const N = 3
    # tspan = (0.0, 1.0)
    # tspan = (0.0, 8.0)
    tspan = (0.0, 0.1)  # benchmark time span
    tol = 1e-10
end

# define ode

begin
    include("util_" * model_type * ".jl")
    include("ode_" * model_type * ".jl")
    function ode_function!(dX, X, p, t)
        x = X[1:(N+1)]
        dx = X[(N+2):end]

        _s = s(x)
        _ds = ds(x, dx)
        _M = M(_s, _ds)
        _F = F(_s, _ds)
        _G = G(x, dx)
        _C = C(x, dx)

        id_M = I(N + 1)
        zero_c = zeros(N + 1, 1)

        _A = _M * [id_M; _G]
        _B = _F - _M * [zero_c; _C]

        # id_X = I(N + 1)
        # zero_X = zeros(N + 1, N + 1)

        # LHS = [id_X zero_X; zero_X _A]
        # RHS = [dx; _B]
        # dX_tmp = LHS \ RHS
        ddx = inv(_A) * _B
        dX[1:N+1] = dx
        dX[N+2:end] = ddx
        # dX_tmp = [dx; ddx]
        # copyto!(dX, dX_tmp)
        nothing
    end
end

begin
    include("bc.jl")
    X0, dX0 = load_init_conditions()
    prob = ODEProblem(ode_function!, X0, tspan)
end

# show important param

N
T() / rd
tspan

if param_type == "ideal" || param_type == "ideal_1p"
    title = model_type * "_N" * string(N) * "_T" * string(Int(T() / rd)) * "_" * param_type
else
    title = model_type * "_N" * string(N) * "_" * param_type
end

# solve the ode

# @time sol = solve(prob, TRBDF2(autodiff=false), reltol=tol, abstol=tol);
@benchmark sol = solve(prob, TRBDF2(autodiff=false), reltol=tol, abstol=tol) samples = 40 seconds = 80

## post process

# generate data points from sol
begin
    plt_tspan = tspan
    # plt_tspan = (0.0, 2.0)
    plot_size = Int(20000 * plt_tspan[2] - 20000 * plt_tspan[1])
    include("post.jl")
    display("post data generated!")
end

# plot1
begin
    # GLMakie.activate!()
    dispnew(figure) = display(GLMakie.Screen(), figure)
    # plot xe
    fig = Figure()
    ax_xe = Axis(fig[1, 1:2], ylabel="xe")
    ax_dxe = Axis(fig[2, 1:2], ylabel="dxe")
    ax_ddxe = Axis(fig[3, 1:2], ylabel="ddxe")

    lines!(ax_xe, t, sol_xe)
    lines!(ax_dxe, t, sol_dxe)
    lines!(ax_ddxe, t, sol_ddxe)

    # dispnew(fig)

    # function eta
    fig1 = Figure()
    ax_eta = Axis(fig1[1, 1:2], ylabel="eta")
    ax_deta = Axis(fig1[2, 1:2], ylabel="deta")

    for eta_id in 1:3
        lines!(ax_eta, t, sol_eta[:, eta_id], label="eta_$eta_id")
        lines!(ax_deta, t, sol_deta[:, eta_id], label="deta_$eta_id")
    end
    axislegend(ax_deta, "Mode Shapes", position=:rt)
    # dispnew(fig1)

    # deformation plot
    fig2 = Figure()

    ax_u1 = Axis(fig2[1, 1:2], ylabel="u(xe-rp[1])")
    ax_u2 = Axis(fig2[2, 1:2], ylabel="u(L/3)")
    lines!(ax_u1, t, sol_u1)
    lines!(ax_u2, t, sol_u2)
    # dispnew(fig2)

    # force plot
    fig3 = Figure()

    ax_cf = Axis(fig3[1, 1:2], ylabel="cable force")
    lines!(ax_cf, t, sol_cable_force)
    dispnew(fig3)
end

# export data

begin
    using MAT
    file = matopen("./data/" * title * ".mat", "w")
    write(file, title * "_xe", sol_xe)
    write(file, title * "_dxe", sol_dxe)
    write(file, title * "_force", sol_cable_force)

    close(file)
    display("post data exported!")
end