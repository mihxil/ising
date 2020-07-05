set nologscale 
set grid
set autoscale
#set title 'aantal spins up en down als functie van het magneetveld'
set title 'aantal spins up en down als functie van 1/temperatuur'
#set xlabel 'B'
set xlabel 'J'
set ylabel 'aantal up/down'
set yrange [0:400]
#plot 'result3' using 2:3 t 'upb' with linespoints, 'result3' using 2:4 t 'downb' with linespoints, 'result4' using 2:3 t 'upm' with linespoints, 'result4' using 2:4 t 'downm' with linespoint
           
plot 'result' using 1:3 t 'upb' with linespoints, 'result' using 1:4 t 'downb' with linespoints
#, 'result1' using 1:3 t 'upm' with linespoints, 'result1' using 1:4 t 'downm' with linespoints