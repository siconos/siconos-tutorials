terminal_type='aqua'
terminal_type='X11'


set  term terminal_type 0
resultfile = 'BouncingBallNETS.dat'
plot \
resultfile every ::2 u 1:2 t "Ball  x position" w l,\
resultfile every ::2 u 1:3 t "Ball  y position" w l,\
resultfile every ::2 u 1:4 t "Ball  z position" w l

set  term terminal_type 1
plot resultfile every ::2 u 1:26 t "Ball  angle rotation " w l

set  term terminal_type 2
set yrange [-1.2:1.2]
plot \
resultfile every ::2 u 1:5 t "Ball  quat_1 " w l,\
resultfile every ::2 u 1:6 t "Ball  quat_2" w l,\
resultfile every ::2 u 1:7 t "Ball  quat_3" w l,\
resultfile every ::2 u 1:8 t "Ball  quat_4" w l

set  term terminal_type 2
set yrange restore
set yrange [-5.0:5.0]
plot \
resultfile every ::2 u 1:9 t "Ball velocity x " w l,\
resultfile every ::2 u 1:10 t "Ball  velocity y" w l,\
resultfile every ::2 u 1:11 t "Ball  velocity z" w l

set  term terminal_type 3
set yrange restore
set yrange [-1.0:1.0]
plot \
resultfile every ::2 u 1:12 t "Ball  rot velocity x" w l,\
resultfile every ::2 u 1:13 t "Ball  rot velocity y" w l,\
resultfile every ::2 u 1:14 t "Ball  rot velocity z" w l

