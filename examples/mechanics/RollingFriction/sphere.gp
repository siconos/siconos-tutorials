
result_velocity ="sphere-velocity-body_1.dat"
result_position ="sphere-position-body_1.dat"
result_velocity ="sphere_transverse-velocity-body_1.dat"
result_velocity_absolute ="sphere_transverse-velocity-absolute-body_1.dat"
result_position ="sphere_transverse-position-body_1.dat"

set term GNUTERM 1
plot result_velocity u 1:2 w l, result_velocity u 1:3 w l, result_velocity u 1:4 w l
set term GNUTERM 2
plot result_velocity u 1:5 w l, result_velocity u 1:6 w l, result_velocity u 1:7 w l
set term GNUTERM 3
plot result_position u 1:2 w l, result_position u 1:3 w l, result_position u 1:4 w l
set term GNUTERM 4
plot result_position u 1:5 w l, result_position u 1:6 w l, result_position u 1:7 w l,  result_position u 1:7 w l
set term GNUTERM 5
plot result_velocity_absolute u 1:5 w l, result_velocity_absolute u 1:6 w l, result_velocity_absolute u 1:7 w l
