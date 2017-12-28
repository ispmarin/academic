#!/usr/bin/gnuplot
#Gnuplot script to plot comparison data between phreeqc and bionapl
#Ivan Marin, Université Laval, 18/12/2012
#ispmarin@gmail.com


set style data lines
#plot Benzene concentration in kg
set xlabel "Time, days"
set ylabel "Concentration, kg/m^3 (g/L)"
set title "Benzene Comparison for different cells"
plot "./concentration1.data" u 1:2, "./concentration2.data" u 1:2, "./concentration3.data" u 1:2,"./concentration4.data" u 1:2
#pause -1
set terminal wxt 2
set title "Electron acceptor (O2) Comparison"
plot "./concentration1.data" u 1:3, "./concentration2.data" u 1:3, "./concentration3.data" u 1:3,"./concentration4.data" u 1:3

set terminal wxt 3
set title "Biomass"
plot "./concentration1.data" u 1:4, "./concentration2.data" u 1:4, "./concentration3.data" u 1:4,"./concentration4.data" u 1:4
pause -1
