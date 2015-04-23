reset

# File path
elmspace = "/home/ash/"
landua= "/1/home/dehgha2/"
mac = "/Users/ashkandehghan/"


path= mac."Desktop/SCFT_CODES/ThinFilm/ThinFilm_StarBlock/RESULTS/"


set autoscale
set key b r

plot   path."MOD1_Hor.dat" using 6:1  title "HOR"with linespoints lt 1,\
 path."MOD1_Ver.dat" using 6:1 title "VER" with linespoints lt 3
    
    pause(-1)

