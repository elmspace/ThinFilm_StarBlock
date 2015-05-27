reset

set term postscript enhanced color
set output "~/Desktop/fE.ps"

set xlabel '{/Symbol D}{/Symbol c}'  font ",32" offset 0,-1
set ylabel 'fE'  font ",32"

set bmargin 5
set lmargin 12

# File path
elmspace = "/home/ash/"
landua= "/1/home/dehgha2/"
mac = "/Users/ashkandehghan/"

path= mac."Desktop/SCFT_CODES/ThinFilm/ThinFilm_StarBlock/DATA/LAM_New/"



set key b r font ",18" spacing 3.5


set xr [0.0 : 1.33]
set yr [3.768 : 3.803]
set ytics 0.01
set xtics 0.25
set xtics font ", 20"
set ytics font ", 20"

f(x) = 0.0

plot path."1Arm_Hor_Ver.dat" using ($6-$4):1  title "Parallel" w lp lc 1 pt 4 lw 5,\
path."1Arm_Hor_Ver.dat" using ($6-$4):9  title "Perpendicular" w lp lc 2 pt 6 lw 5

    pause(-1)

