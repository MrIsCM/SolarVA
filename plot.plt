set terminal gif animate delay 3
set output "Gifs/R100G500.gif"  

rock = 'Datos/Rock.dat'
gas = 'Datos/Gas.dat'


set xrange[-60:60]
set yrange[-60:60]

stats rock nooutput
stats gas nooutput

do for [bb = 0:int(STATS_blocks)-2]{
	set title sprintf('Frame:%03.0f',bb)
	plot 	rock u 2:3 index bb title 'Masas Rocosas',\
			gas u 2:3 index bb pt 7 title 'Masas Gaseosa'

}