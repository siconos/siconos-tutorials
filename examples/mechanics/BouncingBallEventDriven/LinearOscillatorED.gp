#set  term X11
resultfile= 'LinearOscillatorED.dat'
#resultfile= 'BouncingBallED.ref'
#resultfile = 'BouncingBallEDwithRestingContact.ref'
set term aqua 0
plot \
resultfile every ::2  u 1:2 t "Ball position -- q x " w l,\
resultfile every ::2  u 1:3 t "Ball Velocity -- v x" w l,\
resultfile every ::2  u 1:4 t "Impulse -- p(1)" w l,\
resultfile every ::2  u 1:5 t "Reaction force -- p(2)" w l,\
resultfile every ::2  u 1:6 t "Impulse multiplier -- lambda(1)" w l,\
resultfile every ::2  u 1:7 t "Reaction multiplier -- lambda(2)" w l

set term aqua 1
plot \
resultfile every ::2  u 1:8 t "Ball position -- q y " w l,\
resultfile every ::2  u 1:9 t "Ball Velocity -- v y" w l

set term aqua 2
plot \
resultfile every ::2  u 1:2 t "Ball position -- q x " w l