#set  term wx
#!tail -n 100000 result.dat > result-gp.dat
resultfile = "ImpactingStick-FremondNSL.dat"
#resultfile = "ImpactingStick.dat"


basheight = 0.2
heightoff = 0.2
winratio = 1.0
winheight = basheight*winratio



set term aqua 0
set grid
plot \
resultfile every ::2 u 2:3 t "position y vs x " w l


set term aqua 1
set size 1,1	
set origin	0,0
set multiplot layout 3,1 columnsfirst scale 1.0,1.0
plot \
resultfile every ::2 u 1:2 t "position x vs t" w l, \
resultfile every ::2 u 1:3 t "position y vs t" w l, \
resultfile every ::2 u 1:4 t "position theta vs t" w l
plot \
resultfile every ::2 u 1:5 t "velocity x vs t" w l, \
resultfile every ::2 u 1:6 t "velocity y vs t" w l,\
resultfile every ::2 u 1:7 t "velocity theta vs t" w l
plot \
resultfile every ::2 u 1:11 t "velocity u_n vs t " w l,\
resultfile every ::2 u 1:12 t "velocity u_t vs t " w l


unset multiplot


set term aqua 4
plot \
resultfile every ::2 u 13:14 t "cone " w p

set term aqua 5

set size 1,1	
set origin	0,0
set multiplot layout 3,1 columnsfirst scale 1.0,1.0
#     [ up to 6 plot commands here ]
plot \
resultfile every ::2 u 1:15 t "normal work" w l
plot \
resultfile every ::2 u 1:16 t "tangent work" w l
plot \
resultfile every ::2 u 1:17 t "dissipation" w l
unset multiplot


set term aqua 6

set size 1,1	
set origin	0,0
set multiplot layout 3,1 columnsfirst scale 1.0,1.0
#     [ up to 6 plot commands here ]
plot \
resultfile every ::2 u 1:18 t "normal dissipation" w l
plot \
resultfile every ::2 u 1:19 t "tangent dissipation" w l
plot \
resultfile every ::2 u 1:20 t "kinetic energy" w l,\
resultfile every ::2 u 1:21 t "potential energy" w l,\
resultfile every ::2 u 1:($20+$21) t "mechanical energy" w l,\
resultfile every ::2 u 1:22 t "work forces" w l

unset multiplot


