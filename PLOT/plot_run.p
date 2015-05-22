
reset

# File path
elmspace = "/home/ash/"
landua= "/1/home/dehgha2/"
mac="/Users/ashkandehghan/"



x1=4.83609
y1=9.78680

x2=5.84060
y2=9.54263

#set xr [x1 : x2]

m=(y2-y1)/(x2-x1)
b=y2-m*x2

f(x)=m*x+b




path= mac."Desktop/SCFT_CODES/ThinFilm/ThinFilm_StarBlock/RESULTS/"


set autoscale
unset key
#plot   path."fE.dat" using 5:1 with linespoints lt 1
plot   path."fE.dat" using 5:($1-f($5)) with linespoints lt 1


   pause(-1)
