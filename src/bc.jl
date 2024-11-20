# functions
function load_init_conditions()
    q0 = zeros(N + 1, 1)
    # xe
    q0[end] = L
    dq0 = zeros(N + 1, 1)
    # dxe
    dq0[end] = 0.0
    X0 = [q0; dq0]

    dX0 = zeros(2 * N + 2, 1)
    return X0, dX0
end