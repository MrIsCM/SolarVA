!==============================================
!
! 		Simulacion Formacion Planetoides
! 
! 				 30/06/2022
!
!==============================================

program planetesimales
	implicit none

	!----------------------------------------------------------
	! 					PLANETAS ROCOSOS
	!----------------------------------------------------------
	
	! Cantidad de planetesimales iniciales
	integer, parameter :: rock_n = 50

	
	! Masa, posicion, velocidad, acelerock_acion, funcion aux w, 
	double precision :: rock_m(1:rock_n), rock_x(1:rock_n), rock_y(1:rock_n), rock_vx(1:rock_n)
	double precision ::  rock_vy(1:rock_n), rock_ax(1:rock_n), rock_ay(1:rock_n), rock_wx(1:rock_n), rock_wy(1:rock_n)


	!----------------------------------------------------------
	! 					PLANETAS GASEOSO
	!----------------------------------------------------------

	! Cantidad de planetesimales iniciales Tipo: Gas
	integer, parameter :: gas_n = 450

	! Masa, posicion, velocidad, aceleracion, funcion aux w, 
	double precision :: gas_m(1:gas_n), gas_x(1:gas_n), gas_y(1:gas_n), gas_vx(1:gas_n), gas_vy(1:gas_n)
	double precision ::  gas_ax(1:gas_n), gas_ay(1:gas_n), gas_wx(1:gas_n), gas_wy(1:gas_n)


	!----------------------------------------------------------
	! 					PARAMETROS GENERALES
	!----------------------------------------------------------

	! Parametros temporales (tiempo e incremento)
	double precision :: t, h 

	! Declaracion de constantes de utilidad (U.Astro, Cte Gravitacion, Masa Sol)
	double precision :: c, G, Ms

	! Factores de escala (distancia, tiempo, velocidad)
	double precision :: f_r, f_t, f_v

	! Energia y momento angular total
	double precision :: E_t, L_t	

	! Parametro iterativo
	integer :: i
	integer, parameter :: iter = 10000



	!###################################################
	! 		Definicion de algunos parametros
	!###################################################
	c = 1.493E11 		! m
	G = 6.67E-11 		! m^3* kg^-1* s^-2 
	Ms = 1.99E30 		! kg

	f_r = 6.684491979E-12 					! 		1/c 			! m
	f_t = 1.991095381E-7					! (G*Ms/c**3)**(1/2)	! s
	f_v = 3.357193253E-5					! f_r/f_t				! m/s

	t=0
	h=0.001

	open(10, file='Datos/INFO_Simul.txt', status='unknown')
		write(10,*) "Datos simulacion formacion S.Solar - Fortran95"
		write(10,*) "h = 		", h, "		Intervalo de tiempo (~h*58 dias)"
		write(10,*) "i = 		", iter, "		Numero de iteraciones"

	open(1, file='Datos/Rock.txt', status='unknown')
	open(2, file='Datos/Gas.txt', status='unknown')
	open(3, file='Datos/Constantes.txt', status='unknown')



	call condiciones_iniciales(rock_x, rock_y, rock_vx, rock_vy, f_v, rock_n)
	call condiciones_iniciales(gas_x, gas_y, gas_vx, gas_vy, f_v, gas_n)

	do i = 1, rock_n
		write(1,*) rock_x(i), rock_y(i), rock_vy(i)
	end do

	do i = 1, gas_n 
		write(2,*) gas_x(i), gas_y(i), gas_vy(i)
	end do


	close(1)
	close(2)
	close(3)



end program planetesimales


subroutine condiciones_iniciales(x, y, vx, vy, f_v, n)
	implicit none

	integer, intent (in) :: n
	double precision, intent(out) :: f_v  
	double precision,intent(inout) :: x(1:n), y(1:n), vx(1:n), vy(1:n)
	
	integer :: i

	call random_seed()
	call random_number(x)
	call random_number(vx)
	call random_number(vy)

	x = 80*x + 30 			! 30 < x < 110  Unidades Astronomicas
	y = 0 					! Parten del eje x 

	vx = 0 					! Al partir del eje x --> vx = 0 
	vy = (30 + 30*vy)*f_v 	! 30 < vy < 60 		Luego lo reescalo con f_v

	do i = 1, n
		if (mod(i,2) == 0) then 
			vy(i) = -vy(i)		! Quiero que algunos giren en sentido opuesto
		end if 
	end do


end subroutine condiciones_iniciales