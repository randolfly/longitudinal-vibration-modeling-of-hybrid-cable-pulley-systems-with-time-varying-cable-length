file = matopen("./ideal_0pulley_deformation_N3_3D.mat", "w")
write(file, "zp_N3_tall_u", sol_u_all)
write(file, "zp_N3_tall_ux", sol_uₓ_all)
write(file, "zp_N3_tall_t", Array(sol_t_range))
write(file, "zp_N3_tall_x", Array(sol_x_range))
close(file)