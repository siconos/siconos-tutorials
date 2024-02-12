#set  term wx
#!tail -n 100000 result.dat > result-gp.dat
resultfile = "BouncingBallTS-FremondNSL.dat"


basheight = 0.2
heightoff = 0.2
winratio = 1.0
winheight = basheight*winratio




set multiplot
set size winratio,winheight


set origin 0.0,winheight*3.0+heightoff
plot \
resultfile every ::2 u 2:3 t "Ball position y vs x " w l

set origin 0.0,winheight*2.0+heightoff
plot \
resultfile every ::2 u 1:3 t "Ball position y vs t" w l


set origin 0.0,winheight*1.0+heightoff
plot \
resultfile every ::2 u 1:7 t "normal gap  vs t " w l




set origin 0.0,winheight*0.0+heightoff
plot \
resultfile every ::2 u 10:12 t "Ball position u_t vs r_t " w l 0

unset multiplot
     


