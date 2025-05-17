
function export_post_data(title, sol_xe, sol_dxe, sol_cable_force)
    file = matopen("./data/" * title * ".mat", "w")
    write(file, title * "_xe", sol_xe)
    write(file, title * "_dxe", sol_dxe)
    write(file, title * "_force", sol_cable_force)

    close(file)
    display("post data exported!")
end


function show_model_param(mp::ModelParam)
    df1 = DataFrame(
        N=[mp.N],
        tspan=[mp.tspan[2]],
        tol=[mp.tol],
        L=[mp.L],
        rd=[mp.rd],
        Id=[mp.Id],
        m=[mp.m],
        k=[mp.k],
        T_rd=[mp.T / mp.rd],
    )
    df2 = DataFrame(
        Td_rd=[mp.Td / mp.rd],
        Cd=[mp.Cd],
        pulley_num=[mp.pulley_num],
        lp=[mp.lp],
        rp=[mp.rp[1]],
        Ip=[mp.Ip[1]],
        Tp_rp=[mp.Tp[1] / mp.rp[1]],
        Cp=[mp.Cp[1]],
    )
    pretty_table(df1)
    pretty_table(df2)
end

function post_sol(mt::ModelType, mp::ModelParam, sol, export_data=false, title="test")
    # extract sol data into 20000 points
    plt_tspan = (sol.t[1], sol.t[end])
    # plot_size = round(20000 * plt_tspan[2] - 20000 * plt_tspan[1])
    plot_size = 5000
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

        sol_cable_force1 = zeros(plot_size)
        sol_cable_force2 = zeros(plot_size)
        sol_cable_force_tmp = zeros(plot_size)
        # sol_cable_equivalent_stiffness = zeros(plot_size)

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

            sol_cable_force1[i] = get_cable_force(mp, sol_xe[i], sol_ddxe[i])
            # sol_cable_force2[i] = get_cable_force(mt, mp, x, x[end])
            # display(sol_cable_force2[i])
            sol_cable_force_tmp[i] = sol_cable_force1[i] - sol_cable_force2[i]
            # sol_cable_equivalent_stiffness[i] = sol_cable_force[i] / (mp.L - sol_xe[i])
            # if (sol_cable_equivalent_stiffness[i] > 3e3)
            #     # remove nan
            #     sol_cable_equivalent_stiffness[i] = 3e3
            # end
        end

    end
    display("post data generated!")

    # Plots plot
    begin

    end

    # GLMakie plot
    begin
        GLMakie.activate!()
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
        fig1 = Figure()
        ax_eta = Axis(fig1[1, 1:2], ylabel="eta")
        ax_deta = Axis(fig1[2, 1:2], ylabel="deta")

        for eta_id in 1:3
            lines!(ax_eta, t, sol_eta[:, eta_id], label="eta_$eta_id")
            lines!(ax_deta, t, sol_deta[:, eta_id], label="deta_$eta_id")
        end
        axislegend(ax_deta, "Mode Shapes", position=:rt)
        dispnew(fig1)

        # # deformation plot
        fig2 = Figure()

        ax_u1 = Axis(fig2[1, 1:2], ylabel="u(xe-rp[1])")
        ax_u2 = Axis(fig2[2, 1:2], ylabel="u(L/3)")
        lines!(ax_u1, t, sol_u1)
        lines!(ax_u2, t, sol_u2)
        # dispnew(fig2)

        # force plot

        fig3 = Figure()
        ax_cf = Axis(fig3[1, 1:2], ylabel="cable force1")
        # ax_ck = Axis(fig3[2, 1:2], ylabel="cable force2")
        # ax_tmp = Axis(fig3[3, 1:2], ylabel="cable force1+2")
        lines!(ax_cf, t, sol_cable_force1)
        # lines!(ax_ck, t, sol_cable_force2)
        # lines!(ax_tmp, t, sol_cable_force_tmp)
        dispnew(fig3)
    end

    # FullModel or PartialModel 3D deformation data
    begin
        deformation_length_size = 300

        deformation_sol_t_range = LinRange(sol.t[1], sol.t[end], deformation_length_size)
        deformation_sol_x_range = LinRange(0, 1, deformation_length_size)
        # x axis is the coordinate of the cable, y axis is the time, z axis is the deformation/strain
        deformation_sol_u_all = zeros(deformation_length_size, deformation_length_size)
        deformation_sol_uₓ_all = zeros(deformation_length_size, deformation_length_size)

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
                sol_uₓ_all[i] = ∂ₓu(mt, mp, x, ux[i])
            end
            return sol_uₓ_all
        end

        for id in 1:deformation_length_size
            t = deformation_sol_t_range[id]
            deformation_sol_u_all[:, id] = get_all_deformation(t, sol, mp, mt, deformation_length_size)
            deformation_sol_uₓ_all[:, id] = get_all_strains(t, sol, mp, mt, deformation_length_size)
        end

        GLMakie.activate!()
        fig = surface(deformation_sol_x_range, deformation_sol_t_range, deformation_sol_u_all,
            axis=(type=Axis3, azimuth=pi / 4,
                zlabel="deformation", xlabel="x̃", ylabel="t"))

        display(GLMakie.Screen(), fig)

        GLMakie.activate!()
        fig = surface(deformation_sol_x_range, deformation_sol_t_range, deformation_sol_uₓ_all,
            axis=(type=Axis3, azimuth=pi / 4,
                zlabel="strain", xlabel="x̃", ylabel="t"))
        display(GLMakie.Screen(), fig)
        display("3D deformation data generated!")
    end

    # export data
    if export_data
        file = matopen("./data/" * title * ".mat", "w")
        write(file, title * "_xe", sol_xe)
        write(file, title * "_dxe", sol_dxe)
        write(file, title * "_force", sol_cable_force1)

        write(file, title * "_eta1", sol_eta[:, 1])
        write(file, title * "_eta2", sol_eta[:, 2])
        write(file, title * "_eta3", sol_eta[:, 3])
        write(file, title * "_deta1", sol_deta[:, 1])
        write(file, title * "_deta2", sol_deta[:, 2])
        write(file, title * "_deta3", sol_deta[:, 3])

        # 3D deformation
        write(file, title * "_sol_u_all", deformation_sol_u_all)
        write(file, title * "_sol_ux_all", deformation_sol_uₓ_all)
        write(file, title * "_sol_t_range", Array(deformation_sol_t_range))
        write(file, title * "_sol_x_range", Array(deformation_sol_x_range))

        # 2D deformation at t = 0.001 and t = 0.02
        sol_u_t1 = get_all_deformation(0.001, sol, mp, mt)
        sol_uₓ_t1 = get_all_strains(0.001, sol, mp, mt)
        sol_u_t2 = get_all_deformation(0.02, sol, mp, mt)
        sol_uₓ_t2 = get_all_strains(0.02, sol, mp, mt)
        write(file, title * "_sol_u_t1", sol_u_t1)
        write(file, title * "_sol_ux_t1", sol_uₓ_t1)
        write(file, title * "_sol_u_t2", sol_u_t2)
        write(file, title * "_sol_ux_t2", sol_uₓ_t2)

        close(file)
        display("post data exported!")
    end
end