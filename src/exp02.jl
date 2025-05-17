# EXPERIMENT 02: comparision of the extended trial space results and the simplified model results
# COMPARISION MODELS:
# 1. FullModel with extended trial space: full model with all the dynamic parts, discretization N is N=3
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
    N = 3
    tol = 1e-8

    pulley_num = 1
    L = 0.42 + 1.04 * pulley_num

    rd = 1.5e-2
    # Id = 1.7485e-5 + 1.11e-3
    Id = 1.7485e-5
    Td = 0.155 * rd
    Cd = 0.020 * rd

    rp = 1.5e-2
    Ip = 1.7485e-5
    Tp = Td
    Cp = Cd

    T = 48.37 * rd

    k = 1290.44
    m = 0.65

    mt = FullModel()
end

# full model with extended trial function, pulley 1, drum small, decritzation N=3
title = "Exp02_FullEx_P1_DS_N3"
begin
    N = 3
    tspan = (0.0, 0.3)
    pulley_num = 1
    L = 0.42 + 1.04 * pulley_num
    Id = 1.7485e-5
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
    mt = FullModel()
    show_model_param(mp)
    # io = open("solver_stat.txt", "w")
    X0, dX0 = load_init_conditions(mt, mp)
    prob = ODEProblem(ode_function!, X0, mp.tspan, mp)
    @time sol = solve(prob, TRBDF2(autodiff=false), reltol=mp.tol, abstol=mp.tol)
    display("solver finished!")

    post_sol(mt, mp, sol, true, title)
    display("post data exported!")

    export_deformation_video(title)
    display("deformation video exported!")
    # close(io)
    # @benchmark sol = solve(prob, TRBDF2(autodiff=false), reltol=mp.tol, abstol=mp.tol) samples = 30 seconds = 1200
end

# full model with extended trial function, pulley 1, drum big, decritzation N=3(with actual experiment results)
title = "Exp02_FullEx_P1_DB_N3"
begin
    N = 3
    tspan = (0.0, 6.0)
    pulley_num = 1
    L = 0.42 + 1.04 * pulley_num
    Id = 1.7485e-5 + 1.11e-3
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
    mt = FullModel()
    show_model_param(mp)
    # io = open("solver_stat.txt", "w")
    X0, dX0 = load_init_conditions(mt, mp)
    prob = ODEProblem(ode_function!, X0, mp.tspan, mp)
    @time sol = solve(prob, TRBDF2(autodiff=false), reltol=mp.tol, abstol=mp.tol)
    display("solver finished!")

    post_sol(mt, mp, sol, true, title)
    display("post data exported!")

    export_deformation_video(title)
    display("deformation video exported!")
    # close(io)
    # @benchmark sol = solve(prob, TRBDF2(autodiff=false), reltol=mp.tol, abstol=mp.tol) samples = 30 seconds = 1200
end

# full model with extended trial function, pulley 3, drum small, decritzation N=3
title = "Exp02_FullEx_P3_DS_N3"
begin
    N = 3
    tspan = (0.0, 0.3)
    pulley_num = 3
    L = 0.42 + 1.04 * pulley_num
    Id = 1.7485e-5
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
    mt = FullModel()
    show_model_param(mp)
    # io = open("solver_stat.txt", "w")
    X0, dX0 = load_init_conditions(mt, mp)
    prob = ODEProblem(ode_function!, X0, mp.tspan, mp)
    @time sol = solve(prob, TRBDF2(autodiff=false), reltol=mp.tol, abstol=mp.tol)
    display("solver finished!")

    post_sol(mt, mp, sol, true, title)
    display("post data exported!")

    export_deformation_video(title)
    display("deformation video exported!")
    # close(io)
    # @benchmark sol = solve(prob, TRBDF2(autodiff=false), reltol=mp.tol, abstol=mp.tol) samples = 30 seconds = 1200
end

# full model with extended trial function, pulley 3, drum big, decritzation N=3(with actual experiment results)
title = "Exp02_FullEx_P3_DB_N3"
begin
    N = 3
    tspan = (0.0, 3.0)
    pulley_num = 3
    L = 0.42 + 1.04 * pulley_num
    Id = 1.7485e-5 + 1.11e-3
    mt = FullModel()
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
    show_model_param(mp)
    # io = open("solver_stat.txt", "w")
    X0, dX0 = load_init_conditions(mt, mp)
    prob = ODEProblem(ode_function!, X0, mp.tspan, mp)
    @time sol = solve(prob, TRBDF2(autodiff=false), reltol=mp.tol, abstol=mp.tol)
    display("solver finished!")

    post_sol(mt, mp, sol, true, title)
    display("post data exported!")

    export_deformation_video(title)
    display("deformation video exported!")
    # close(io)
    # @benchmark sol = solve(prob, TRBDF2(autodiff=false), reltol=mp.tol, abstol=mp.tol) samples = 30 seconds = 1200
end

# ======================= SIMPLE MODEL RESULTS ==================================
# simple model with extended trial function, pulley 1, drum small, decritzation N=3, concertrate pulley-induced forces
title = "Exp02_Simple_P1_DS_N3_PC"
begin
    N = 3
    tspan = (0.0, 0.3)
    pulley_num = 1
    L = 0.42 + 1.04 * pulley_num
    Td = 0.155 * rd * (pulley_num + 1)
    Cd = 0.020 * rd * (pulley_num + 1)
    Id = 1.7485e-5 * (pulley_num + 1)
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
    mt = SimpleModel()
    show_model_param(mp)
    # io = open("solver_stat.txt", "w")
    X0, dX0 = load_init_conditions(mt, mp)
    prob = ODEProblem(ode_function!, X0, mp.tspan, mp)
    @time sol = solve(prob, TRBDF2(autodiff=false), reltol=mp.tol, abstol=mp.tol)
    display("solver finished!")

    post_sol(mt, mp, sol, true, title)
    display("post data exported!")

    export_deformation_video(title)
    display("deformation video exported!")
    # close(io)
    # @benchmark sol = solve(prob, TRBDF2(autodiff=false), reltol=mp.tol, abstol=mp.tol) samples = 30 seconds = 1200
end

# simple model with extended trial function, pulley 1, drum big, decritzation N=3, concertrate pulley-induced forces
title = "Exp02_Simple_P1_DB_N3_PC"
begin
    N = 3
    tspan = (0.0, 6.0)
    pulley_num = 1
    L = 0.42 + 1.04 * pulley_num
    Td = 0.155 * rd * (pulley_num + 1)
    Cd = 0.020 * rd * (pulley_num + 1)
    Id = + 1.11e-3 + 1.7485e-5 * (pulley_num + 1)
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
    mt = SimpleModel()
    show_model_param(mp)
    # io = open("solver_stat.txt", "w")
    X0, dX0 = load_init_conditions(mt, mp)
    prob = ODEProblem(ode_function!, X0, mp.tspan, mp)
    @time sol = solve(prob, TRBDF2(autodiff=false), reltol=mp.tol, abstol=mp.tol)
    display("solver finished!")

    post_sol(mt, mp, sol, true, title)
    display("post data exported!")

    export_deformation_video(title)
    display("deformation video exported!")
    # close(io)
    # @benchmark sol = solve(prob, TRBDF2(autodiff=false), reltol=mp.tol, abstol=mp.tol) samples = 30 seconds = 1200
end

# simple model with extended trial function, pulley 3, drum small, decritzation N=3, concertrate pulley-induced forces
title = "Exp02_Simple_P3_DS_N3_PC"
begin
    N = 3
    tspan = (0.0, 0.3)
    pulley_num = 3
    L = 0.42 + 1.04 * pulley_num
    Td = 0.155 * rd * (pulley_num + 1)
    Cd = 0.020 * rd * (pulley_num + 1)
    Id = 1.7485e-5 * (pulley_num + 1)
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
    mt = SimpleModel()
    show_model_param(mp)
    # io = open("solver_stat.txt", "w")
    X0, dX0 = load_init_conditions(mt, mp)
    prob = ODEProblem(ode_function!, X0, mp.tspan, mp)
    @time sol = solve(prob, TRBDF2(autodiff=false), reltol=mp.tol, abstol=mp.tol)
    display("solver finished!")

    post_sol(mt, mp, sol, true, title)
    display("post data exported!")

    export_deformation_video(title)
    display("deformation video exported!")
    # close(io)
    # @benchmark sol = solve(prob, TRBDF2(autodiff=false), reltol=mp.tol, abstol=mp.tol) samples = 30 seconds = 1200
end

# simple model with extended trial function, pulley 3, drum big, decritzation N=3, concertrate pulley-induced forces
title = "Exp02_Simple_P3_DB_N3_PC"
begin
    N = 3
    tspan = (0.0, 3.0)
    pulley_num = 3
    L = 0.42 + 1.04 * pulley_num
    Td = 0.155 * rd * (pulley_num + 1)
    Cd = 0.020 * rd * (pulley_num + 1)
    Id = 1.11e-3 +  1.7485e-5 * (pulley_num + 1)
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
    mt = SimpleModel()
    show_model_param(mp)
    # io = open("solver_stat.txt", "w")
    X0, dX0 = load_init_conditions(mt, mp)
    prob = ODEProblem(ode_function!, X0, mp.tspan, mp)
    @time sol = solve(prob, TRBDF2(autodiff=false), reltol=mp.tol, abstol=mp.tol)
    display("solver finished!")

    post_sol(mt, mp, sol, true, title)
    display("post data exported!")

    export_deformation_video(title)
    display("deformation video exported!")
    # close(io)
    # @benchmark sol = solve(prob, TRBDF2(autodiff=false), reltol=mp.tol, abstol=mp.tol) samples = 30 seconds = 1200
end

# ----- IGNORE PULLEY INDUCED FORCES -----

# simple model with extended trial function, pulley 1, drum small, decritzation N=3, ignore pulley-induced forces
title = "Exp02_Simple_P1_DS_N3_PI"
begin
    N = 3
    tspan = (0.0, 0.3)
    pulley_num = 1
    L = 0.42 + 1.04 * pulley_num
    Td = 0.155 * rd * (1)
    Cd = 0.020 * rd * (1)
    Id = 1.7485e-5 * (1)
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
    mt = SimpleModel()
    show_model_param(mp)
    # io = open("solver_stat.txt", "w")
    X0, dX0 = load_init_conditions(mt, mp)
    prob = ODEProblem(ode_function!, X0, mp.tspan, mp)
    @time sol = solve(prob, TRBDF2(autodiff=false), reltol=mp.tol, abstol=mp.tol)
    display("solver finished!")

    post_sol(mt, mp, sol, true, title)
    display("post data exported!")

    export_deformation_video(title)
    display("deformation video exported!")
    # close(io)
    # @benchmark sol = solve(prob, TRBDF2(autodiff=false), reltol=mp.tol, abstol=mp.tol) samples = 30 seconds = 1200
end

# simple model with extended trial function, pulley 1, drum big, decritzation N=3, ignore pulley-induced forces
title = "Exp02_Simple_P1_DB_N3_PI"
begin
    N = 3
    tspan = (0.0, 6.0)
    pulley_num = 1
    L = 0.42 + 1.04 * pulley_num
    Td = 0.155 * rd * (1)
    Cd = 0.020 * rd * (1)
    Id = + 1.11e-3 + 1.7485e-5 * (1)
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
    mt = SimpleModel()
    show_model_param(mp)
    # io = open("solver_stat.txt", "w")
    X0, dX0 = load_init_conditions(mt, mp)
    prob = ODEProblem(ode_function!, X0, mp.tspan, mp)
    @time sol = solve(prob, TRBDF2(autodiff=false), reltol=mp.tol, abstol=mp.tol)
    display("solver finished!")

    post_sol(mt, mp, sol, true, title)
    display("post data exported!")

    export_deformation_video(title)
    display("deformation video exported!")
    # close(io)
    # @benchmark sol = solve(prob, TRBDF2(autodiff=false), reltol=mp.tol, abstol=mp.tol) samples = 30 seconds = 1200
end

# simple model with extended trial function, pulley 3, drum small, decritzation N=3, ignore pulley-induced forces
title = "Exp02_Simple_P3_DS_N3_PI"
begin
    N = 3
    tspan = (0.0, 0.3)
    pulley_num = 3
    L = 0.42 + 1.04 * pulley_num
    Td = 0.155 * rd * (1)
    Cd = 0.020 * rd * (1)
    Id = 1.7485e-5 * (1)
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
    mt = SimpleModel()
    show_model_param(mp)
    # io = open("solver_stat.txt", "w")
    X0, dX0 = load_init_conditions(mt, mp)
    prob = ODEProblem(ode_function!, X0, mp.tspan, mp)
    @time sol = solve(prob, TRBDF2(autodiff=false), reltol=mp.tol, abstol=mp.tol)
    display("solver finished!")

    post_sol(mt, mp, sol, true, title)
    display("post data exported!")

    export_deformation_video(title)
    display("deformation video exported!")
    # close(io)
    # @benchmark sol = solve(prob, TRBDF2(autodiff=false), reltol=mp.tol, abstol=mp.tol) samples = 30 seconds = 1200
end

# simple model with extended trial function, pulley 3, drum big, decritzation N=3, ignore pulley-induced forces
title = "Exp02_Simple_P3_DB_N3_PI"
begin
    N = 3
    tspan = (0.0, 3.0)
    pulley_num = 3
    L = 0.42 + 1.04 * pulley_num
    Td = 0.155 * rd * (1)
    Cd = 0.020 * rd * (1)
    Id = 1.11e-3 +  1.7485e-5 * (1)
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
    mt = SimpleModel()
    show_model_param(mp)
    # io = open("solver_stat.txt", "w")
    X0, dX0 = load_init_conditions(mt, mp)
    prob = ODEProblem(ode_function!, X0, mp.tspan, mp)
    @time sol = solve(prob, TRBDF2(autodiff=false), reltol=mp.tol, abstol=mp.tol)
    display("solver finished!")

    post_sol(mt, mp, sol, true, title)
    display("post data exported!")

    export_deformation_video(title)
    display("deformation video exported!")
    # close(io)
    # @benchmark sol = solve(prob, TRBDF2(autodiff=false), reltol=mp.tol, abstol=mp.tol) samples = 30 seconds = 1200
end

# ----------- N=10 -------------
# simple model with extended trial function, pulley 1, drum small, decritzation N=10, concertrate pulley-induced forces
title = "Exp02_Simple_P1_DS_N10_PC"
begin
    N = 10
    tspan = (0.0, 0.3)
    pulley_num = 1
    L = 0.42 + 1.04 * pulley_num
    Td = 0.155 * rd * (pulley_num + 1)
    Cd = 0.020 * rd * (pulley_num + 1)
    Id = 1.7485e-5 * (pulley_num + 1)
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
    mt = SimpleModel()
    show_model_param(mp)
    # io = open("solver_stat.txt", "w")
    X0, dX0 = load_init_conditions(mt, mp)
    prob = ODEProblem(ode_function!, X0, mp.tspan, mp)
    @time sol = solve(prob, TRBDF2(autodiff=false), reltol=mp.tol, abstol=mp.tol)
    display("solver finished!")

    post_sol(mt, mp, sol, true, title)
    display("post data exported!")

    export_deformation_video(title)
    display("deformation video exported!")
    # close(io)
    # @benchmark sol = solve(prob, TRBDF2(autodiff=false), reltol=mp.tol, abstol=mp.tol) samples = 30 seconds = 1200
end

# simple model with extended trial function, pulley 1, drum big, decritzation N=10, concertrate pulley-induced forces
title = "Exp02_Simple_P1_DB_N10_PC"
begin
    N = 10
    tspan = (0.0, 6.0)
    pulley_num = 1
    L = 0.42 + 1.04 * pulley_num
    Td = 0.155 * rd * (pulley_num + 1)
    Cd = 0.020 * rd * (pulley_num + 1)
    Id = + 1.11e-3 + 1.7485e-5 * (pulley_num + 1)
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
    mt = SimpleModel()
    show_model_param(mp)
    # io = open("solver_stat.txt", "w")
    X0, dX0 = load_init_conditions(mt, mp)
    prob = ODEProblem(ode_function!, X0, mp.tspan, mp)
    @time sol = solve(prob, TRBDF2(autodiff=false), reltol=mp.tol, abstol=mp.tol)
    display("solver finished!")

    post_sol(mt, mp, sol, true, title)
    display("post data exported!")

    export_deformation_video(title)
    display("deformation video exported!")
    # close(io)
    # @benchmark sol = solve(prob, TRBDF2(autodiff=false), reltol=mp.tol, abstol=mp.tol) samples = 30 seconds = 1200
end

# simple model with extended trial function, pulley 3, drum small, decritzation N=10, concertrate pulley-induced forces
title = "Exp02_Simple_P3_DS_N10_PC"
begin
    N = 10
    tspan = (0.0, 0.3)
    pulley_num = 3
    L = 0.42 + 1.04 * pulley_num
    Td = 0.155 * rd * (pulley_num + 1)
    Cd = 0.020 * rd * (pulley_num + 1)
    Id = 1.7485e-5 * (pulley_num + 1)
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
    mt = SimpleModel()
    show_model_param(mp)
    # io = open("solver_stat.txt", "w")
    X0, dX0 = load_init_conditions(mt, mp)
    prob = ODEProblem(ode_function!, X0, mp.tspan, mp)
    @time sol = solve(prob, TRBDF2(autodiff=false), reltol=mp.tol, abstol=mp.tol)
    display("solver finished!")

    post_sol(mt, mp, sol, true, title)
    display("post data exported!")

    export_deformation_video(title)
    display("deformation video exported!")
    # close(io)
    # @benchmark sol = solve(prob, TRBDF2(autodiff=false), reltol=mp.tol, abstol=mp.tol) samples = 30 seconds = 1200
end

# simple model with extended trial function, pulley 3, drum big, decritzation N=10, concertrate pulley-induced forces
title = "Exp02_Simple_P3_DB_N10_PC"
begin
    N = 10
    tspan = (0.0, 3.0)
    pulley_num = 3
    L = 0.42 + 1.04 * pulley_num
    Td = 0.155 * rd * (pulley_num + 1)
    Cd = 0.020 * rd * (pulley_num + 1)
    Id = 1.11e-3 +  1.7485e-5 * (pulley_num + 1)
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
    mt = SimpleModel()
    show_model_param(mp)
    # io = open("solver_stat.txt", "w")
    X0, dX0 = load_init_conditions(mt, mp)
    prob = ODEProblem(ode_function!, X0, mp.tspan, mp)
    @time sol = solve(prob, TRBDF2(autodiff=false), reltol=mp.tol, abstol=mp.tol)
    display("solver finished!")

    post_sol(mt, mp, sol, true, title)
    display("post data exported!")

    export_deformation_video(title)
    display("deformation video exported!")
    # close(io)
    # @benchmark sol = solve(prob, TRBDF2(autodiff=false), reltol=mp.tol, abstol=mp.tol) samples = 30 seconds = 1200
end

# ----- IGNORE PULLEY INDUCED FORCES -----

# simple model with extended trial function, pulley 1, drum small, decritzation N=10, ignore pulley-induced forces
title = "Exp02_Simple_P1_DS_N10_PI"
begin
    N = 10
    tspan = (0.0, 0.3)
    pulley_num = 1
    L = 0.42 + 1.04 * pulley_num
    Td = 0.155 * rd * (1)
    Cd = 0.020 * rd * (1)
    Id = 1.7485e-5 * (1)
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
    mt = SimpleModel()
    show_model_param(mp)
    # io = open("solver_stat.txt", "w")
    X0, dX0 = load_init_conditions(mt, mp)
    prob = ODEProblem(ode_function!, X0, mp.tspan, mp)
    @time sol = solve(prob, TRBDF2(autodiff=false), reltol=mp.tol, abstol=mp.tol)
    display("solver finished!")

    post_sol(mt, mp, sol, true, title)
    display("post data exported!")

    export_deformation_video(title)
    display("deformation video exported!")
    # close(io)
    # @benchmark sol = solve(prob, TRBDF2(autodiff=false), reltol=mp.tol, abstol=mp.tol) samples = 30 seconds = 1200
end

# simple model with extended trial function, pulley 1, drum big, decritzation N=10, ignore pulley-induced forces
title = "Exp02_Simple_P1_DB_N10_PI"
begin
    N = 10
    tspan = (0.0, 6.0)
    pulley_num = 1
    L = 0.42 + 1.04 * pulley_num
    Td = 0.155 * rd * (1)
    Cd = 0.020 * rd * (1)
    Id = + 1.11e-3 + 1.7485e-5 * (1)
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
    mt = SimpleModel()
    show_model_param(mp)
    # io = open("solver_stat.txt", "w")
    X0, dX0 = load_init_conditions(mt, mp)
    prob = ODEProblem(ode_function!, X0, mp.tspan, mp)
    @time sol = solve(prob, TRBDF2(autodiff=false), reltol=mp.tol, abstol=mp.tol)
    display("solver finished!")

    post_sol(mt, mp, sol, true, title)
    display("post data exported!")

    export_deformation_video(title)
    display("deformation video exported!")
    # close(io)
    # @benchmark sol = solve(prob, TRBDF2(autodiff=false), reltol=mp.tol, abstol=mp.tol) samples = 30 seconds = 1200
end

# simple model with extended trial function, pulley 3, drum small, decritzation N=10, ignore pulley-induced forces
title = "Exp02_Simple_P3_DS_N10_PI"
begin
    N = 10
    tspan = (0.0, 0.3)
    pulley_num = 3
    L = 0.42 + 1.04 * pulley_num
    Td = 0.155 * rd * (1)
    Cd = 0.020 * rd * (1)
    Id = 1.7485e-5 * (1)
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
    mt = SimpleModel()
    show_model_param(mp)
    # io = open("solver_stat.txt", "w")
    X0, dX0 = load_init_conditions(mt, mp)
    prob = ODEProblem(ode_function!, X0, mp.tspan, mp)
    @time sol = solve(prob, TRBDF2(autodiff=false), reltol=mp.tol, abstol=mp.tol)
    display("solver finished!")

    post_sol(mt, mp, sol, true, title)
    display("post data exported!")

    export_deformation_video(title)
    display("deformation video exported!")
    # close(io)
    # @benchmark sol = solve(prob, TRBDF2(autodiff=false), reltol=mp.tol, abstol=mp.tol) samples = 30 seconds = 1200
end

# simple model with extended trial function, pulley 3, drum big, decritzation N=10, ignore pulley-induced forces
title = "Exp02_Simple_P3_DB_N10_PI"
begin
    N = 10
    tspan = (0.0, 3.0)
    pulley_num = 3
    L = 0.42 + 1.04 * pulley_num
    Td = 0.155 * rd * (1)
    Cd = 0.020 * rd * (1)
    Id = 1.11e-3 +  1.7485e-5 * (1)
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
    mt = SimpleModel()
    show_model_param(mp)
    # io = open("solver_stat.txt", "w")
    X0, dX0 = load_init_conditions(mt, mp)
    prob = ODEProblem(ode_function!, X0, mp.tspan, mp)
    @time sol = solve(prob, TRBDF2(autodiff=false), reltol=mp.tol, abstol=mp.tol)
    display("solver finished!")

    post_sol(mt, mp, sol, true, title)
    display("post data exported!")

    export_deformation_video(title)
    display("deformation video exported!")
    # close(io)
    # @benchmark sol = solve(prob, TRBDF2(autodiff=false), reltol=mp.tol, abstol=mp.tol) samples = 30 seconds = 1200
end

# generate dispalcement video for all time steps
function export_deformation_video(title="test")
    frame_rate = 24
    sol_t_range = LinRange(sol.t[1], sol.t[end], frame_rate * 5)
    length_size = 1000

    function get_all_deformation(sol_t, sol, mp, mt, length_size=1000)
        sol_u_all = zeros(length_size)
        x = sol(sol_t)[1:mp.N+mp.pulley_num+1]
        ux = LinRange(0, x[end], length_size)
        for i in 1:length_size
            sol_u_all[i] = u(mt, mp, x, ux[i]) * 1000 # unit: mm
        end
        return sol_u_all
    end

    function get_all_strains(sol_t, sol, mp, mt, length_size=1000)
        sol_uₓ_all = zeros(length_size)
        x = sol(sol_t)[1:mp.N+mp.pulley_num+1]
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
    record(fig, title * "_deformation.mp4", sol_t_range; framerate=frame_rate) do t
        # display(t)
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