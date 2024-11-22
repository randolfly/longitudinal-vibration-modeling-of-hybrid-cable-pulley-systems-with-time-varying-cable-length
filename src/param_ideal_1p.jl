# cable params
E = 78.0e9

cable_diameter = 1.2e-3
A = π * (cable_diameter / 2)^2
#  A = 1.0
ρ = 7.5e-3
#  ρ = 1.0

L = 5.0
a = sqrt(A * E / ρ)

# winch params
rd = 1.5e-2
ld = 0.0  # residual length to xe
Id = 1.7485e-5

# pulley params
pulley_num = 1
Ipi = 1.7485e-5;
Ip = ones(pulley_num) * Ipi
lp = 0.42 .+ 1.04 .* [i for i in 0:(pulley_num-1)]'
rp = ones(pulley_num) * rd

# mass params
m = 0.65
k = 100

# friction params
Td = 0.0 * rd
Cd = 0.0 * rd
Tp_i = Td
Cp_i = Cd
Tp = ones(pulley_num) * Tp_i
Cp = ones(pulley_num) * Cp_i
Tm = 0.0
Cm = 0.0

# driven force
function T()
    return 200 * rd
end