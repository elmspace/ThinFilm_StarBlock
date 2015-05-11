reset

# File path
elmspace = "/home/ash/"
landua= "/1/home/dehgha2/"
mac="/Users/ashkandehghan/"

path= mac."Desktop/SCFT_CODES/ThinFilm/ThinFilm_StarBlock/RESULTS/"


set autoscale
unset key
plot   path."Run.dat" using 1:3 with linespoints lt 1

   pause(-1)


   reset

# File path
elmspace = "/home/ash/"
landua= "/1/home/dehgha2/"
mac="/Users/ashkandehghan/"

path= mac."Desktop/SCFT_CODES/ThinFilm/ThinFilm_StarBlock/RESULTS/"


set autoscale
unset key
plot   path."fE.dat" using 5:1 with linespoints lt 1

   pause(-1)
