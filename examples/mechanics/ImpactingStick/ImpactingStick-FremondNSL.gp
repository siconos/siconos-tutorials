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
plot \
resultfile every ::2 u 1:2 t "position x vs t" w l, \
resultfile every ::2 u 1:3 t "position y vs t" w l, \
resultfile every ::2 u 1:4 t "position theta vs t" w l


set term aqua 2
plot \
resultfile every ::2 u 1:5 t "velocity x vs t" w l, \
resultfile every ::2 u 1:6 t "velocity y vs t" w l,\
resultfile every ::2 u 1:7 t "velocity theta vs t" w l




#set origin 0.0,winheight*0.0+heightoff
set term aqua 3
plot \
resultfile every ::2 u 1:11 t "velocity u_n vs t " w l,\
resultfile every ::2 u 1:12 t "velocity u_t vs t " w l,\

set term aqua 4
plot \
resultfile every ::2 u 13:14 t "cone " w p


#unset multiplot
     


