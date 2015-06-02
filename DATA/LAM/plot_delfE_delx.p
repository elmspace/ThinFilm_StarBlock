reset

#set term postscript enhanced color
#set output "~/Desktop/fE.ps"

mac = "/Users/ashkandehghan/"

path= mac."Desktop/SCFT_CODES/ThinFilm/ThinFilm_StarBlock/DATA/LAM/"


set xlabel '{/Symbol D}{/Symbol c}N'  font ",32" offset 0,-1
set ylabel '{/Symbol D}fE'  font ",32"

set xr [0.0 : 100.0]
set yr [-0.15 : 0.1]
set key t r font ",24" spacing 2.5
set pointsize 2

f(x) = 0

plot f(x) notitle " " w l lc -1 lw 4,\
path."Lam_1Arm.dat" using ($1-$2):($3-$4)  title "1Arm" with linespoints lt 1 pt 2 lw 4,\
   path."Lam_2Arm.dat" using (($1-$2)/2):(($3-$4)/2)  title "2Arm" with linespoints lt 2 pt 4 lw 4,\
       path."Lam_3Arm.dat" using (($1-$2)/3):(($3-$4)/3)  title "3Arm" with linespoints lt 3 pt 6 lw 4,\
	  path."Lam_4Arm.dat" using (($1-$2)/4):(($3-$4)/4)  title "4Arm" with linespoints lt 4 pt 8 lw 4
    
    pause(-1)

