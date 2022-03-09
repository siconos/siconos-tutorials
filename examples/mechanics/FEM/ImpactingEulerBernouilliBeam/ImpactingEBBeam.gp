set  term X11
set term pdf
set term pdf monochrome dashed size 5in,3.5in
set term pdf color size 5in,3.5in
extension= ".pdf"
#set term tikz standalone color solid size 5in,3in
#set term tikz standalone monochrome  size 5in,3in font '\small\sf'
#extension=".tex"


resultfile='ImpactingEBBeam.dat'
method="MoreauJean-theta05-"

model="ImpactingEBBeam"
timestep = "1e-4"
h=1e-4
ndof = "2-"


file = model.method.ndof.timestep


set output file."-Beam-Energy".extension
set auto
plot \
     resultfile u 1:6 t "Beam Potential energy" w l,\
     resultfile u 1:7 t "Beam Kinetic energy" w l,\
     resultfile u 1:($7+$6) t "Beam Total energy" w l

set output file."-Total-Energy".extension
set auto
plot \
     resultfile u 1:6 t "Beam Potential energy" w l,\
     resultfile u 1:($7+$17) t "Kinetic energy" w l,\
     resultfile u 1:($7+$6+$17) t "Total energy" w l



set output file."-ImpactEnergy".extension
set auto
plot \
     resultfile u 1:12 t "Impact energy" w l


set output file."-Caplambda".extension
plot \
     resultfile u 1:4 t "Reaction impulse" w l

set output file."-ball-position".extension
plot \
      resultfile u 1:($13) t "Ball position" w l ,\
      resultfile u 1:($2) t " Beam end position" w l ,\
      
set output file."-lambda".extension
plot \
     resultfile u 1:($4) t "lambda" w l
set output file."-u".extension
plot \
     resultfile u 1:($14) t "u" w l ,\
     resultfile u 1:($15) t "u_{k+theta}" w l




set output file."-qvr".extension
offset=0.05
xoff=0.05
set format y "%.1e"

basheight = 0.3
heightoff = 0.07
winratio = 0.95
winheight = basheight*winratio

set lmargin 8 
set bmargin 0
set tmargin 0

set size 1.0 , 1.2

set xzeroaxis
set format x ""

set mxtics 
set grid xtics
set grid mxtics

set multiplot
set size winratio,winheight

set ylabel "m" 
set origin xoff,winheight*2.0+heightoff+offset+offset
plot \
     resultfile u (1000*$1):2 t "q" w l
set ylabel "m/s"
set origin xoff,winheight+heightoff+offset
plot \
     resultfile u (1000*$1):3 t "v" w l
set format x "%.2f"
#set bmargin 0
set xtics border
set xlabel "time [ms]"
set ylabel "N"

set origin xoff,0.0+heightoff
plot \
     resultfile u (1000*$1):($4) t "r" w l


