set terminal pngcairo dashed size 800,500
set output "fig2.png"
set termoption font "Helvetica,18"
set key spacing 1.2
set lmargin 9
set rmargin 5
set bmargin 3.5
set xlabel "axial distance (non-dim)"
set ylabel "CP factor"

set yrange [1:1.3]
plot 'conc.dat' u ($1/0.0762):($2/32) w l lw 4 dt 1 title "NMACFoam",\
     '../literature_data/permeatevel_4umpersec_Subramani_JMS_277_2006.dat' u 1:2 w lp lc -1 ps 3 pt 5 dt 1 title "Subramani et al., JMS, 277, 2006"
