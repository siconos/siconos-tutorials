set  term X11
set term pdf
set output "SimpleExampleRelay.pdf"
resultfile="SimpleExampleRelay_py.dat"
resultfile="SimpleExampleRelay.dat"
set xrange [0:1.2]
set yrange [-1.2:1.2]
plot \
resultfile u 1:2 t "Siconos Platform -- INRIA                                             Process x(1)" w l,\
resultfile u 1:3 t "Process x(2)" w l,\
resultfile u 1:6 t "lambda(1)" w l,\
resultfile u 1:7 t "lambda(2)" w l
