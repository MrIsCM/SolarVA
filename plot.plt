set terminal gif animate delay 10
set output "Gifs/rocoso.gif"  

rock = 'Datos/Rock.dat'

set xrange[-100:120]
set yrange[-50:50]

do for [bb = 0:200]{
	set title sprintf('Frame:%03.0f',bb)
	plot rock every ::::bb u 2:3
}