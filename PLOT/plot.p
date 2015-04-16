reset

# File path
elmspace = "/home/ash/"
landua= "/1/home/dehgha2/"


path= landua."Desktop/SCFT_CODES/ThinFilm/ThinFilm_StarBlock/PHI/"


set autoscale
unset key
plot   path."phi_x.dat" using 1:($2+$3+$4+$5) with linespoints lt 1,\
    path."phi_x.dat" using 1:($6+$7+$8+$9) with linespoints lt 3

