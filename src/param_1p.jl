# cable params
# const E = 137.0e7
const E = 78.0e9
# const E = 10.0

cable_diameter = 1.2e-3
const A = π * (cable_diameter / 2)^2
# const A = 1.0
const ρ = 7.5e-3
# const ρ = 1.0

const L = 1.40
const a = sqrt(A * E / ρ)

# winch params
const rd = 1.5e-2
const ld = 0.0  # residual length to xe
const Id = 1.7485e-5 + 1.11e-3

# pulley params
const Ip = [1.7485e-5; 1.7485e-5]
const lp = [0.26; 1.11]    # residual length to xe
const rp = [rd; rd; rd; rd]
const pulley_num = length(Ip)

# mass params
const m = 0.65
# const k = 100.0
const k = 1290.44


# friction params
const Td = 0.4 * rd
const Cd = 0.0 * rd
const Tp_i = Td
const Cp_i = Cd
const Tp = [Tp_i; Tp_i; Tp_i]
const Cp = [Cp_i; Cp_i; Cp_i]
const Tm = 0.0
const Cm = 0.0

# driven force
function T()
    return 0.72555  # 48.37*rd
end