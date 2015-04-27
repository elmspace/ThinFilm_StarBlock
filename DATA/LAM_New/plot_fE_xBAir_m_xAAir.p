reset

#set term postscript enhanced color
#set output "~/Desktop/fE.ps"

set xlabel '{/Symbol D}{/Symbol c}'  font ",32" offset 0,-1
set ylabel 'fE'  font ",32"

set bmargin 5
set lmargin 12

# File path
elmspace = "/home/ash/"
landua= "/1/home/dehgha2/"
mac = "/Users/ashkandehghan/"

path= landua."Desktop/SCFT_CODES/ThinFilm/ThinFilm_StarBlock/DATA/LAM_New/"



set key at 6.5,9.89 font ",24" spacing 3.0


#set xr [-1.6 : 7.7]
#set yr [9.66 : 9.9]
set ytics 0.025
set xtics 1.0
set xtics font ", 18"
set ytics font ", 18"

f(x) = 0.0

plot path."4Arm_Hor.dat" using ($6-$4):1  title "Parallel" w lp lc 1 pt 4 lw 4,\
path."4Arm_Ver.dat" using ($6-$4):1  title "Perpendicular" w lp lc 2 pt 6 lw 4

    pause(-1)

