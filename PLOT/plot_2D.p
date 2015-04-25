reset

# File path
elmspace = "/home/ash/"
landua= "/1/home/dehgha2/"
mac = "/Users/ashkandehghan/"

path= mac."Desktop/SCFT_CODES/ThinFilm/ThinFilm_StarBlock/PHI/"


set pm3d
set iso 100
set samp 100
set palette model RGB
set dgrid3d 20,20,1
set pm3d flush begin ftriangles scansforward interpolate 10,5

set xlabel 'x'  font ",22"
set ylabel 'y'  font ",22"
   
unset key
unset sur
set hidden3d
set view map 
set autoscale
set size square
# Format: index1 index2 A1 A2 A3 A4 B1 B2 B3 B4
splot path."phi_xy.dat" using 1:2:($3+$4+$5+$6)

pause(-1)
   
#####################################################################################################
reset

# File path
elmspace = "/home/ash/"
landua= "/1/home/dehgha2/"
mac = "/Users/ashkandehghan/"
   
path= mac."Desktop/SCFT_CODES/ThinFilm/ThinFilm_StarBlock/PHI/"


set pm3d
set iso 100
set samp 100
set palette model RGB
set dgrid3d 20,20,1
set pm3d flush begin ftriangles scansforward interpolate 10,5

set xlabel 'x'  font ",22"
set ylabel 'z'  font ",22"
unset key
unset sur
set hidden3d
set view map 
set autoscale
set size square
# Format: index1 index2 A1 A2 A3 A4 B1 B2 B3 B4
splot path."phi_xz.dat" using 1:2:($3+$4+$5+$6)

pause(-1)

#####################################################################################################

reset

# File path
elmspace = "/home/ash/"
landua= "/1/home/dehgha2/"
mac = "/Users/ashkandehghan/"
   
path= mac."Desktop/SCFT_CODES/ThinFilm/ThinFilm_StarBlock/PHI/"


set pm3d
set iso 100
set samp 100
set palette model RGB
set dgrid3d 20,20,1
set pm3d flush begin ftriangles scansforward interpolate 10,5

set xlabel 'y'  font ",22"
set ylabel 'z'  font ",22"
unset key
unset sur
set hidden3d
set view map 
set autoscale
set size square
# Format: index1 index2 A1 A2 A3 A4 B1 B2 B3 B4
splot path."phi_yz.dat" using 1:2:($3+$4+$5+$6)
