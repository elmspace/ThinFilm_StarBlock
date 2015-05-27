reset

set term postscript enhanced color
set output "~/Desktop/fE.ps"

mac = "/Users/ashkandehghan/"

path= mac."Desktop/SCFT_CODES/ThinFilm/ThinFilm_StarBlock/DATA/LAM/"


set xlabel '{/Symbol D}{/Symbol c}N'  font ",32" offset 0,-1
set ylabel 'fE'  font ",32"

set autoscale
set key b r font ",24" spacing 3.0
set pointsize 2

plot   path."Lam_4Arm.dat" using ($1-$2):3  title "Parallel" with linespoints lt 1 pt 4 lw 4,\
   path."Lam_4Arm.dat" using ($1-$2):4  title "Perpendicular" with linespoints lt 2 pt 6 lw 4,\
      path."Lam_4Arm.dat" using ($1-$2):5  title "Mixed" with linespoints lt 3 pt 8 lw 4
    
    pause(-1)

