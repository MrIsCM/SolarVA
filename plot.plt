set terminal gif animate delay 3
set output "Gifs/R_SunProx.gif"  

rock = 'Datos/Rock.dat'

set xrange[-60:60]
set yrange[-60:60]

do for [bb = 0:800]{
	set title sprintf('Frame:%03.0f',bb)
	plot rock u 2:3 index bb 
}