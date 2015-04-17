reset

# File path
elmspace = "/home/ash/"
landua= "/1/home/dehgha2/"


path= landua."Desktop/SCFT_CODES/ThinFilm/ThinFilm_StarBlock/DATA/"


set autoscale
set key b r

plot   path."1Arm/HEX/HOR_VER.dat" using 1:2  title "HOR"with linespoints lt 1,\
 path."1Arm/HEX/HOR_VER.dat" using 1:3 title "VER" with linespoints lt 3

