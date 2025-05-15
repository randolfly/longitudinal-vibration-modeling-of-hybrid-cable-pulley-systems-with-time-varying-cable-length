# include the necessary packages
begin
    using DifferentialEquations
    using GLMakie
    using LinearAlgebra
    using DelimitedFiles
    using ForwardDiff
    using BenchmarkTools, Profile
end

begin
    E = 1.0
    A = 1.0
    L = 1.0
    rho = 1.0
    a = sqrt(E / rho)
    f = 1.0
    l = 0.3
end

begin
    const N = 10
    tspan = (0.0, 3)
    tol = 1e-15
end

begin
    function ode_function!(dX, X, p, t)
        x = X[1:N]
        dx = X[N+1:end]

        Kx = zeros(N, N)
        for i in 1:N
            Kx[i, i] = i^2 * pi^2 / 2
        end

        Fnon = zeros(N, 1)
        for i in 1:N
            Fnon[i] = sin(i * pi * l / L)
        end

        ddx = 2 * (f * Fnon - (a^2 / L^2) * (Kx * x))
        dX[1:N] .= dx
        dX[N+1:end] .= ddx
        nothing
    end
end

begin
    x0 = zeros(N, 1)
    dx0 = zeros(N, 1)
    X0 = [x0; dx0]

    prob = ODEProblem(ode_function!, X0, tspan)
end

@time sol = solve(prob, TRBDF2(autodiff=false), reltol=tol, abstol=tol);

begin
    # resample plot points
    plot_size = 10000
    t = LinRange(tspan..., plot_size)

    sol_eta = zeros(plot_size, N)
    sol_deta = zeros(plot_size, N)
    sol_ddeta = zeros(plot_size, N)

    sol_u1 = zeros(plot_size)
    sol_u2 = zeros(plot_size)
    sol_u3 = zeros(plot_size)
    sol_u4 = zeros(plot_size)
    sol_u5 = zeros(plot_size)

    sol_ux1 = zeros(plot_size)
    sol_ux2 = zeros(plot_size)
    sol_ux3 = zeros(plot_size)
    sol_ux4 = zeros(plot_size)
    sol_ux5 = zeros(plot_size)


    function numerical_derivative(sol, t)
        ForwardDiff.derivative(t -> sol(t), t)
    end

    function u(l::Real, eta::Vector{Float64})
        _u = 0
        for i in 1:N
            _u += eta[i] * sin(i * pi * l / L)
        end
        return _u
    end

    function ux(l::Real, eta::Vector{Float64})
        _ux = 0.0
        for i in 1:N
            _ux += eta[i] * i * pi * cos(i * pi * l / L) / L
        end
        return _ux
    end

    for i in 1:plot_size
        dX = numerical_derivative(sol, t[i])
        x = sol(t[i])[1:N]
        dx = sol(t[i])[N+1:end]
        sol_eta[i, :] = x[1:N]
        sol_deta[i, :] = dx[1:N]

        sol_u1[i] = u(0.1, x)
        sol_u2[i] = u(0.2, x)
        sol_u3[i] = u(0.3, x)
        sol_u4[i] = u(0.4, x)
        sol_u5[i] = u(0.5, x)


        sol_ux1[i] = ux(0.1, x)
        sol_ux2[i] = ux(0.2, x)
        sol_ux3[i] = ux(0.3, x)
        sol_ux4[i] = ux(0.4, x)
        sol_ux5[i] = ux(0.5, x)
    end
end

begin
    # GLMakie.activate!()
    dispnew(figure) = display(GLMakie.Screen(), figure)
    fig1 = Figure()
    ax_u = Axis(fig1[1, 1:2], ylabel="u",xlabel="t")
    lines!(ax_u, t, sol_u1, label="u(0.1)")
    lines!(ax_u, t, sol_u2, label="u(0.2)")
    lines!(ax_u, t, sol_u3, label="u(0.3)")
    lines!(ax_u, t, sol_u4, label="u(0.4)")
    lines!(ax_u, t, sol_u5, label="u(0.5)")
    axislegend(ax_u, "deformation", position=:rt)
    dispnew(fig1)

    fig2 = Figure()
    ax_ux = Axis(fig2[1, 1:2], ylabel="strain")
    lines!(ax_ux, t, sol_ux1, label="u(0.1)")
    lines!(ax_ux, t, sol_ux2, label="u(0.2)")
    lines!(ax_ux, t, sol_ux3, label="u(0.3)")
    lines!(ax_ux, t, sol_ux4, label="u(0.4)")
    lines!(ax_ux, t, sol_ux5, label="u(0.5)")
    axislegend(ax_ux, "strain", position=:rt)
    dispnew(fig2)
end