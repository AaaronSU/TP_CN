set terminal png
set output 'Courbe de la convergence.png'

set title "La courbe de convergence des différents méthodes"
set xlabel "Le nombre d'itération"
set ylabel "Erreur résiduelle"

plot 'RESVEC_alpha.dat' with lines title 'Richardson alpha', 'RESVEC_GS.dat' with lines title 'Gauss-Seidel', 'RESVEC_jacobi.dat' with lines title 'Jacobi'
