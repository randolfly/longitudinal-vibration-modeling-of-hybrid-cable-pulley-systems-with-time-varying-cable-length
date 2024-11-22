# cable params
#  E = 137.0e7
E = 78.0e9
#  E = 10.0

cable_diameter = 1.2e-3
A = π * (cable_diameter / 2)^2
#  A = 1.0
ρ = 7.5e-3
#  ρ = 1.0

L = 0.42 + 1.04 * 3
a = sqrt(A * E / ρ)

# winch params
rd = 1.5e-2
ld = 0.0  # residual length to xe
Id = 1.7485e-5 + 1.11e-3

# pulley params
Ip = [1.7485e-5; 1.7485e-5; 1.7485e-5]
lp = [0.42; 0.42 + 1.04; 0.42 + 1.04 * 2]    # residual length to xe
rp = [rd; rd; rd]
pulley_num = length(Ip)

# mass params
m = 0.65
#  k = 100.0
k = 1290.44

# friction params
#  Td = 0.433 * rd
#  Cd = 0.000 * rd
Td = 0.155 * rd
Cd = 0.020 * rd
Tp_i = Td
Cp_i = Cd
Tp = [Tp_i; Tp_i; Tp_i]
Cp = [Cp_i; Cp_i; Cp_i]
Tm = 0.0
Cm = 0.0

# driven force
function T()
    return 48.37 * rd
end