reset

# File path
elmspace = "/home/ash/"
landua= "/1/home/dehgha2/"
mac = "/Users/ashkandehghan/"

path= elmspace."Desktop/SCFT_CODES/ThinFilm/ThinFilm_StarBlock/DATA/LAM/"


set autoscale
set key t r

#set xr [-0.5 : 2.0]

f(x) = 0.0

plot f(x) notitle " " w l lc -1,\
path."1Arm_Hor_Ver.dat" using ($6-$4):($1-$9)  title "1 Arm" w lp lc 1,\
    path."2Arm_Hor_Ver.dat" using ($6-$4):($1-$9)  title "2 Arm" w lp lc 2,\
         path."3Arm_Hor_Ver.dat" using ($6-$4):($1-$9)  title "3 Arm" w lp lc 3,\
	      path."4Arm_Hor_Ver.dat" using ($6-$4):($1-$9)  title "4 Arm" w lp lc 4


#path."1Arm_Hor_Ver.dat" using ($14-$12):9 title "VER" with linespoints lt 3
    
    pause(-1)

