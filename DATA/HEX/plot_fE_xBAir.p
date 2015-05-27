reset

set term postscript enhanced color
set output "~/Desktop/fE.ps"

numb_arms=4

mac = "/Users/ashkandehghan/"

path= mac."Desktop/SCFT_CODES/ThinFilm/ThinFilm_StarBlock/DATA/HEX/"

path_Hor_2D= mac."Desktop/SCFT_CODES/ThinFilm/ThinFilm_StarBlock/RESULTS/"


set xlabel '{/Symbol D}{/Symbol c}N'  font ",32" offset 0,-1
set ylabel 'fE'  font ",32"

set autoscale
set key b r font ",24" spacing 3.0
set pointsize 2
set xr [0.0 : 200]

	   
plot   path."Hex_".numb_arms."Arm.dat" using (($1-$2)*100*numb_arms):3  title "Parallel" with linespoints lt 1 pt 2 lw 6,\
   path."Hex_".numb_arms."Arm.dat" using (($1-$2)*100*numb_arms):4  title "Perpendicular" with linespoints lt 3 pt 4 lw 6,\
      path."Hex_".numb_arms."Arm.dat" using (($1-$2)*100*numb_arms):5  title "Mixed 1" with linespoints lt 2 pt 6 lw 6,\
	 path."Hex_".numb_arms."Arm.dat" using (($1-$2)*100*numb_arms):6  title "Mixed 2" with linespoints lt 4 pt 8 lw 6


#	 path_Hor_2D."MOD1.dat" using ($6-$4):1  title "HOR" with linespoints lt 4





	 
    
    pause(-1)

