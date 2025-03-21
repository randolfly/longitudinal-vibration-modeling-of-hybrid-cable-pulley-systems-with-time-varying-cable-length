
function export_post_data(title, sol_xe, sol_dxe, sol_cable_force)
    file = matopen("./data/" * title * ".mat", "w")
    write(file, title * "_xe", sol_xe)
    write(file, title * "_dxe", sol_dxe)
    write(file, title * "_force", sol_cable_force)

    close(file)
    display("post data exported!")
end

function show_model_param(mp::ModelParam)
    data = ["N" "tspan" "tol" "pulley_num" "T" "T/rd";
        mp.N mp.tspan[2] mp.tol mp.pulley_num mp.T mp.T/mp.rd]
    pretty_table(data)
end

function post_sol(mp::ModelParam, sol, export_data=false)
    # extract sol data into 20000 points
    plt_tspan = mp.tspan
    plot_size = Int(20000 * plt_tspan[2] - 20000 * plt_tspan[1])
    begin
        # resample plot points
        t = LinRange(plt_tspan..., plot_size)
        sol_xe = zeros(plot_size)
        sol_dxe = zeros(plot_size)
        sol_ddxe = zeros(plot_size)

        sol_eta = zeros(plot_size, mp.N)
        sol_deta = zeros(plot_size, mp.N)
        # sol_ddeta = zeros(plot_size, N)

        sol_u1 = zeros(plot_size)
        sol_u2 = zeros(plot_size)

        sol_cable_force = zeros(plot_size)

        # X = [dq; q]; q=[eta(1:N); xe]
        for i in 1:plot_size
            dX = numerical_derivative(sol, t[i])
            # sol_ddeta[i, :] = dX[N+2:2*N+1]
            sol_ddxe[i] = dX[end]

            x = sol(t[i])[1:mp.N+1]
            dx = sol(t[i])[mp.N+2:end]

            sol_eta[i, :] = x[1:mp.N]
            sol_xe[i] = x[end]
            sol_deta[i, :] = dx[1:mp.N]
            sol_dxe[i] = dx[end]

            sol_u1[i] = u(mt, mp, x, sol_xe[i])
            sol_u2[i] = u(mt, mp, x, mp.L / 3)

            sol_cable_force[i] = get_cable_force(mp, sol_xe[i], sol_ddxe[i])
        end
    end
    display("post data generated!")

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
        dispnew(fig)

        # function eta
        # fig1 = Figure()
        # ax_eta = Axis(fig1[1, 1:2], ylabel="eta")
        # ax_deta = Axis(fig1[2, 1:2], ylabel="deta")

        # for eta_id in 1:3
        #     lines!(ax_eta, t, sol_eta[:, eta_id], label="eta_$eta_id")
        #     lines!(ax_deta, t, sol_deta[:, eta_id], label="deta_$eta_id")
        # end
        # axislegend(ax_deta, "Mode Shapes", position=:rt)
        # dispnew(fig1)

        # # deformation plot
        # fig2 = Figure()

        # ax_u1 = Axis(fig2[1, 1:2], ylabel="u(xe-rp[1])")
        # ax_u2 = Axis(fig2[2, 1:2], ylabel="u(L/3)")
        # lines!(ax_u1, t, sol_u1)
        # lines!(ax_u2, t, sol_u2)
        # dispnew(fig2)

        # force plot

        fig3 = Figure()
        ax_cf = Axis(fig3[1, 1:2], ylabel="cable force")
        lines!(ax_cf, t, sol_cable_force)
        dispnew(fig3)
    end

    # export data
    title = string(mp.pulley_num) * "_pulley_" * string(mp.N) * "_N" * string(mp.T) * "_T"
    if export_data
        export_post_data(title, sol_xe, sol_dxe, sol_cable_force)
    end
end


