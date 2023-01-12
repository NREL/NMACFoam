set terminal pngcairo dashed size 800,500
set output "fig1.png"
set termoption font "Helvetica,18"
set key spacing 1.2
set lmargin 9
set rmargin 5
set bmargin 3.5
set xlabel "x/h"
set ylabel "pressure drop (Pa)"
set yrange [0:3000]
plot 'pres.dat' u ($1*2/0.001):(2012.13-$2) w l lw 2 title "NMACFoam",\
'../literature_data/two_porous_walls_Berman_JApplPhys_24_1232_1953.dat' u 1:($2*1000) w p ps 2 pt 5 title "Berman J Appl. Phys. 24 1953"
