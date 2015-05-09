reset

# File path
elmspace = "/home/ash/"
landua= "/1/home/dehgha2/"
mac="/Users/ashkandehghan/"

path= mac."Desktop/SCFT_CODES/ThinFilm/ThinFilm_StarBlock/PHI/"


set autoscale
set key t c
plot   path."phi_z.dat" using 1:($2+$3+$4+$5) title "A" with linespoints lt 1,\
    path."phi_z.dat" using 1:($6+$7+$8+$9) title "B" with linespoints lt 3,\
       path."phi_z.dat" using 1:10 title "HA" with linespoints lt 4,\
	  path."phi_z.dat" using 1:11 title "HS" with linespoints lt 5

