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
	integer, parameter :: rock_n = 10

	
	! Masa, posicion, velocidad, acelerock_acion, funcion aux w, 
	double precision :: rock_m(1:rock_n), rock_x(1:rock_n), rock_y(1:rock_n), rock_vx(1:rock_n)
	double precision ::  rock_vy(1:rock_n), rock_ax(1:rock_n), rock_ay(1:rock_n), rock_wx(1:rock_n), rock_wy(1:rock_n)
	double precision :: rock_Emec


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
	integer :: i, j
	integer, parameter :: iter = 100000



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
	h=0.01

	open(10, file='Datos/INFO_Simul.txt', status='unknown')
		write(10,*) "Datos simulacion formacion S.Solar - Fortran95"
		write(10,*) "h = 		", h, "		Intervalo de tiempo (~h*58 dias)"
		write(10,*) "i = 		", iter, "		Numero de iteraciones"

	open(1, file='Datos/Rock.dat', status='unknown')
	open(2, file='Datos/Gas.dat', status='unknown')
	open(3, file='Datos/Constantes.dat', status='unknown')



	call condiciones_iniciales(rock_x, rock_y, rock_vx, rock_vy, f_v, rock_n)
	call condiciones_iniciales(gas_x, gas_y, gas_vx, gas_vy, f_v, gas_n)

	call aceleraciones(rock_x, rock_y, rock_m, rock_ax, rock_ay, rock_n)

	do i = 0, iter
		call vervelet_pos(rock_x, rock_y, rock_vx, rock_vy, rock_ax, rock_ay, rock_n, h)
		call w(rock_vx, rock_vy, rock_ax, rock_ay, rock_wx, rock_wy, rock_n, h)
		call vervelet_vel(rock_vx, rock_vy, rock_ax, rock_ay, rock_wx, rock_wy, rock_n, h)
		call aceleraciones(rock_x, rock_y, rock_m, rock_ax, rock_ay, rock_n)

		if (mod(i,5*100)==0) then
			do j = 1, rock_n
				write(1,*) t, rock_x(j), rock_y(j)
			end do
			write(1,*) 
			write(1,*)

			call Ener_mec(rock_x, rock_y, rock_vx, rock_vy, rock_m, rock_Emec, rock_n)
			write(3,*) t, rock_Emec
		end if 

		t = t+h
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
	! call random_number(vx)
	call random_number(vy)

	x = 80*x + 30 			! 30 < x < 110  Unidades Astronomicas
	y = 0 					! Parten del eje x 

	! vx = vx*1.1*f_v				! Al partir del eje x considero velocidades muy pequeñas 
	vy = (50 + 100*vy)*f_v	 	! 30 < vy < 60 		Luego lo reescalo con f_v

	! do i = 1, n
	! 	if (mod(i,4) == 0) then 
	! 		vy(i) = -vy(i)		! Quiero que algunos giren en sentido opuesto
	! 	end if 
	! end do


end subroutine condiciones_iniciales



!-------------------------------------------------------------------
!	Solo me interesa la interaccion de las masas con el Sol
!	NO considero las interacciones gravitatorias entre las masas
!-------------------------------------------------------------------
subroutine aceleraciones(x, y, m, ax, ay, n)
	implicit none

	integer, intent(in) :: n
	double precision, intent(in) :: x(1:n), y(1:n), m(1:n)
	double precision, intent(out) :: ax(1:n), ay(1:n)

	double precision, parameter :: G = 6.67E-11

	integer :: i, j

	ax = 0
	ay = 0
	
	do i = 1, n
		ax(i) = ax(i) - (x(i))/(x(i)**2 + y(i)**2)**(3.0/2.0)
		ay(i) = ay(i) - (y(i))/(x(i)**2 + y(i)**2)**(3.0/2.0)
	end do

end subroutine aceleraciones


subroutine vervelet_pos(x, y, vx, vy, ax, ay, n, h)
	implicit none

	integer, intent(in) :: n
	double precision,intent(in) :: h, ax(1:n), ay(1:n)
	double precision, intent(inout) :: x(1:n), y(1:n), vx(1:n), vy(1:n)

	integer :: i

	do i = 1, n
		x(i) = x(i) + h*vx(i) + 0.5*h**2*ax(i)
		y(i) = y(i) + h*vy(i) + 0.5*h**2*ay(i)
	end do

end subroutine vervelet_pos


subroutine w(vx, vy, ax, ay, wx, wy, n, h)
	implicit none

	integer, intent(in) :: n
	double precision, intent(in) :: h, vx(1:n), vy(1:n), ax(1:n), ay(1:n)
	double precision, intent(inout) :: wx(1:n), wy(1:n)
	
	integer :: i
	
	do i = 1, n
		wx(i) = vx(i) + 0.5*h*ax(i)
		wy(i) = vy(i) + 0.5*h*ay(i)
	end do
end subroutine w



subroutine vervelet_vel(vx, vy, ax, ay, wx, wy, n, h)
	implicit none

	integer, intent(in) :: n
	double precision, intent(in) :: h,wx(1:n), wy(1:n), ax(1:n), ay(1:n)
	double precision, intent(out) :: vx(1:n), vy(1:n)

	integer :: i 

	do i = 1, n
		vx(i) = wx(i) + 0.5*h*ax(i)
		vy(i) = wy(i) + 0.5*h*ay(i) 
	end do

end subroutine vervelet_vel


!-----------------------------------------------------
! 		ENERGIA MECANICA
!
!	No tiene en cuenta las perdidas en las colisiones
! 	inelasticas entre las masas
!
!-----------------------------------------------------

subroutine Ener_mec(x, y, vx, vy, m, E_mec, n)
	implicit none

	integer, intent(in) :: n
	double precision, intent(in) :: x(1:10), y(1:10), vx(1:10), vy(1:10), m(1:10)
	double precision, intent(out) :: E_mec
	double precision :: T, V

	integer :: i, j

	T = 0
	V = 0

	do i = 1, n
		T = T + 0.5*m(i)*(vx(i)**2 + vy(i)**2)
	end do
	do i = 1, n
		V = V - m(i)/(x(i)**2+y(i)**2)**0.5
	end do

	E_mec = T + V

	
end subroutine Ener_mec