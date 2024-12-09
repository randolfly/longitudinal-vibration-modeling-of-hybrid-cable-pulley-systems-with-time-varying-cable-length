begin
    # resample plot points

    t = LinRange(plt_tspan..., plot_size)
    sol_xe = zeros(plot_size)
    sol_dxe = zeros(plot_size)
    sol_ddxe = zeros(plot_size)

    sol_eta = zeros(plot_size, N)
    sol_deta = zeros(plot_size, N)
    # sol_ddeta = zeros(plot_size, N)

    sol_u1 = zeros(plot_size)
    sol_u2 = zeros(plot_size)

    sol_cable_force = zeros(plot_size)

    # X = [dq; q]; q=[eta(1:N); xe]
    for i in 1:plot_size
        dX = numerical_derivative(sol, t[i])
        # sol_ddeta[i, :] = dX[N+2:2*N+1]
        sol_ddxe[i] = dX[end]

        x = sol(t[i])[1:N+1]
        dx = sol(t[i])[N+2:end]

        sol_eta[i, :] = x[1:N]
        sol_xe[i] = x[end]
        sol_deta[i, :] = dx[1:N]
        sol_dxe[i] = dx[end]

        sol_u1[i] = u(x, sol_xe[i])
        sol_u2[i] = u(x, L / 3)

        sol_cable_force[i] = get_cable_force(sol_xe[i], sol_ddxe[i])
    end

end