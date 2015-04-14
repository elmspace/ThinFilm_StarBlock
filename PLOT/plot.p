reset

elmspace = "/home/ash/Desktop/SCFT_CODES/ThinFilm/ThinFilm_StarBlock/PHI/"


set autoscale
unset key
plot   elmspace."phi_x.dat" using 1:($2+$3+$4+$5) with linespoints lt 1,\
    elmspace."phi_x.dat" using 1:($6+$7+$8+$9) with linespoints lt 3

