set terminal gif animate delay 3
set output "Gifs/Rocoso3.gif"  

rock = 'Datos/Rock.dat'

set xrange[-120:120]
set yrange[-120:120]

do for [bb = 0:1000]{
	set title sprintf('Frame:%03.0f',bb)
	plot rock u 2:3 index bb 
}