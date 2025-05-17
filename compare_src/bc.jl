using .HybridCableModel

function load_init_conditions(mp::ModelParam)
    q0 = zeros(mp.N + 1, 1)
    # xe
    q0[end] = mp.L
    dq0 = zeros(mp.N + 1, 1)
    # dxe
    dq0[end] = 0.0
    X0 = [q0; dq0]

    dX0 = zeros(2 * mp.N + 2, 1)
    return X0, dX0
end