#ZADANIE2
set term png
set out 'Zadanie2.png'
set logscale x
plot 'Zad2_1.1.dat' u 1:2 w l, 'Zad2_1.2.dat' u 1:2 w l,'Zad2_1.3.dat' u 1:2 w l,'Zad2_1.4.dat' u 1:2 w l,'Zad2_1.5.dat' u 1:2 w l,'Zad2_1.6.dat' u 1:2 w l, 'Zad2_1.7.dat' u 1:2 w l,'Zad2_1.9.dat' u 1:2 w l

#ZADANIE2
set term png
set out 'Zadanie2_iteracje.png'
unset logscale x
plot 'zad2_count.dat' u 1:2 w l,

#ZADANIE1
set term png
set out 'Zadanie1.png'
set logscale x
plot 'Zad1_0.1.dat' u 1:2 w l, 'Zad1_0.2.dat' u 1:2 w l,'Zad1_0.3.dat' u 1:2 w l,'Zad1_0.4.dat' u 1:2 w l,'Zad1_0.5.dat' u 1:2 w l,'Zad1_0.6.dat' u 1:2 w l, 'Zad1_0.7.dat' u 1:2 w l,'Zad1_0.8.dat' u 1:2 w l,'Zad1_0.9.dat' u 1:2 w l

set out 'Zadanie1_gestosc.png'
unset logscale x
set pm3d map
splot 'zad1_gestosc.dat' u 1:2:3

set out 'Zadanie1_C.png'
unset logscale x
set pm3d map
splot 'zad1_C.dat' u 1:2:3

set out 'Zadanie2_C.png'
unset logscale x
set pm3d map
splot 'zad2_C.dat' u 1:2:3

#ZADANIE1
set term png
set out 'Zadanie1_iteracje.png'
unset logscale x
plot 'zad1_count.dat' u 1:2 w l,

set out 'Zadanie2_gestosc.png'
unset logscale x
set pm3d map
splot 'zad2_gestosc.dat' u 1:2:3

#ZADANIE3
set term png
set out 'Zadanie3.png'
set logscale x
plot 'Zad3.dat' u 1:2 w l