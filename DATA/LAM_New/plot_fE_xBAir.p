reset


# File path
elmspace = "/home/ash/"
landua= "/1/home/dehgha2/"
mac = "/Users/ashkandehghan/"

path= elmspace."Desktop/SCFT_CODES/ThinFilm/ThinFilm_StarBlock/DATA/LAM/"



f(x) = 0.0

plot f(x) lc -1,\
path."1Arm_Hor_Ver.dat" using ($6/106/0.25):($1-$9)  title "Parallel" w lp lc 1 pt 4 lw 4,\
path."2Arm_Hor_Ver.dat" using ($6/204/0.25):($1-$9)  title "Perpendicular" w lp lc 2 pt 6 lw 4,\
path."3Arm_Hor_Ver.dat" using ($6/302/0.25):($1-$9)  title "Perpendicular" w lp lc 3 pt 6 lw 4,\
path."4Arm_Hor_Ver.dat" using ($6/400/0.25):($1-$9)  title "Perpendicular" w lp lc 4 pt 6 lw 4

    pause(-1)

