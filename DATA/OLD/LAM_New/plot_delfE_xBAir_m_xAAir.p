reset

set term postscript enhanced color
set output "~/Desktop/fE.ps"


# File path
elmspace = "/home/ash/"
landua= "/1/home/dehgha2/"
mac = "/Users/ashkandehghan/"

path= mac."Desktop/SCFT_CODES/ThinFilm/ThinFilm_StarBlock/DATA/LAM_New/"


set xlabel '{/Symbol D}{/Symbol c}'  font ",32" offset 0,-1
set ylabel '{/Symbol D}fE'  font ",32"

set autoscale
set key t r font ",24" spacing 3.0

set xr [-0.5 : 2.5]
set ytics 0.02
set xtics 0.4
set xtics font ", 18"
set ytics font ", 18"

f(x) = 0.0

plot f(x) notitle " " w l lc -1,\
path."1Arm_Hor_Ver.dat" using ($6-$4):($1-$9)  title "1 Arm" w lp lc 1 lw 4,\
    path."2Arm_Hor_Ver.dat" using ($6-$4):($1-$9)  title "2 Arm" w lp lc 2 lw 4 ,\
         path."3Arm_Hor_Ver.dat" using ($6-$4):($1-$9)  title "3 Arm" w lp lc 3 lw 4,\
	      path."4Arm_Hor_Ver.dat" using ($6-$4):($1-$9)  title "4 Arm" w lp lc 4 lw 4


#path."1Arm_Hor_Ver.dat" using ($14-$12):9 title "VER" with linespoints lt 3
    
    pause(-1)

