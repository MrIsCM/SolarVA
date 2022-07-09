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
	integer, parameter :: rock_n = 100

	
	! Masa, posicion, velocidad, acelerock_acion, funcion aux w, 
	double precision :: rock_m(1:rock_n), rock_r(1:rock_n), rock_x(1:rock_n), rock_y(1:rock_n), rock_vx(1:rock_n)
	double precision ::  rock_vy(1:rock_n), rock_ax(1:rock_n), rock_ay(1:rock_n), rock_wx(1:rock_n), rock_wy(1:rock_n)
	double precision :: rock_Emec, rock_ESun, rock_Qcol, rock_Et


	!----------------------------------------------------------
	! 					PLANETAS GASEOSO
	!----------------------------------------------------------

	! Cantidad de planetesimales iniciales Tipo: Gas
	integer, parameter :: gas_n = 900

	! Masa, posicion, velocidad, aceleracion, funcion aux w, 
	double precision :: gas_m(1:gas_n), gas_r(1:gas_n), gas_x(1:gas_n), gas_y(1:gas_n), gas_vx(1:gas_n), gas_vy(1:gas_n)
	double precision ::  gas_ax(1:gas_n), gas_ay(1:gas_n), gas_wx(1:gas_n), gas_wy(1:gas_n)
	double precision :: gas_Emec, gas_ESun, gas_Qcol, gas_Et 


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
	integer, parameter :: iter = 8*1000000



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
	h=0.007

	open(10, file='Datos/INFO_Simul.txt', status='unknown')
		write(10,*) "Datos simulacion formacion S.Solar - Fortran95"
		write(10,*) "h = 		", h, "		Intervalo de tiempo (~h*58 dias)"
		write(10,*) "i = 		", iter, "		Numero de iteraciones"

	open(1, file='Datos/Rock.dat', status='unknown')
	open(2, file='Datos/Gas.dat', status='unknown')
	open(3, file='Datos/Constantes.dat', status='unknown')



	call condiciones_iniciales_rock(rock_m, rock_r, rock_x, rock_y, rock_vx, rock_vy, rock_n)
	call condiciones_iniciales_gas(gas_m, gas_r, gas_x, gas_y, gas_vx, gas_vy, gas_n)

	call aceleraciones(rock_x, rock_y, rock_m, rock_ax, rock_ay, rock_n)
	call aceleraciones(gas_x, gas_y, gas_m, gas_ax, gas_ay, gas_n)


	! Energia 'destruida' al dejar de considerar las masas que se acercan demasiado a Sol
	rock_ESun = 0 
	rock_Qcol = 0

	do i = 0, iter

		! Deteccion de colisiones entre ROCOSOS
		call colision_rock(rock_r, rock_x, rock_y, rock_vx, rock_vy, rock_m, rock_Qcol, rock_n)

		! Deteccion de colisiones entre GASEOSO
		call colision_gas(gas_r, gas_x, gas_y, gas_vx, gas_vy, gas_m, gas_Qcol, gas_n)

		! 'Eliminacion' de masas por distancia al Sol  (ROCOSO)
		call SunDist(rock_m, rock_x, rock_y, rock_vx, rock_vy, rock_ESun, rock_n)

		! 'Eliminacion' de masas por distancia al Sol (GASEOSO)
		call SunDist(gas_m, gas_x, gas_y, gas_vx, gas_vy, gas_ESun, gas_n)


		! Algoritmo de Vervelet - Calculo de la posicion en t+h - Masas ROCOSAS
		call vervelet_pos(rock_m, rock_x, rock_y, rock_vx, rock_vy, rock_ax, rock_ay, rock_n, h)
		call w(rock_m, rock_vx, rock_vy, rock_ax, rock_ay, rock_wx, rock_wy, rock_n, h)
		call vervelet_vel(rock_m, rock_vx, rock_vy, rock_ax, rock_ay, rock_wx, rock_wy, rock_n, h)
		call aceleraciones(rock_x, rock_y, rock_m, rock_ax, rock_ay, rock_n)

		! Algoritmo de Vervelet - Calculo de la posicion en t+h - Masas GASEOSAS
		call vervelet_pos(gas_m, gas_x, gas_y, gas_vx, gas_vy, gas_ax, gas_ay, gas_n, h)
		call w(gas_m, gas_vx, gas_vy, gas_ax, gas_ay, gas_wx, gas_wy, gas_n, h)
		call vervelet_vel(gas_m, gas_vx, gas_vy, gas_ax, gas_ay, gas_wx, gas_wy, gas_n, h)
		call aceleraciones(gas_x, gas_y, gas_m, gas_ax, gas_ay, gas_n)

		if (mod(i,3000)==0) then
			do j = 1, rock_n
				if (rock_m(j) /= 0) then
					write(1,*) t, rock_x(j), rock_y(j), rock_m(j), rock_r(j)
				end if 
			end do
			do j = 1, gas_n
				if (gas_m(j) /= 0) then
					write(2,*) t, gas_x(j), gas_y(j), gas_m(j), gas_r(j)
				end if 
			end do
			write(1,*) 
			write(1,*)
			write(2,*) 
			write(2,*)

			call Ener_mec(rock_x, rock_y, rock_vx, rock_vy, rock_m, rock_Emec, rock_n)
			call Ener_mec(gas_x, gas_y, gas_vx, gas_vy, gas_m, gas_Emec, gas_n)
			rock_Et = rock_Emec + rock_ESun + rock_Qcol
			gas_Et = gas_Emec + gas_ESun + gas_Qcol

			write(3,*) t, rock_Et + gas_Et, rock_Et, gas_Et
		end if 

		t = t+h
	end do

	close(1)
	close(2)
	close(3)

	j=0
	do i = 1, rock_n
		if (rock_m(i)/=0) then
			j = j+1
		end if
	end do

	write(*,*) 'Quedan: ', j, ' masas ROCOSAS'

	j=0
	do i = 1, gas_n
		if (gas_m(i)/=0) then
			j = j+1
		end if
	end do

	write(*,*) 'Quedan: ', j, ' masas GASEOSAS'


	



end program planetesimales


subroutine condiciones_iniciales_rock(m, r, x, y, vx, vy, n)

	!===========================================================================
	!
	!	Subrutina para adjudicar las condiciones iniciales.
	!
	!			+ El cojunto de posiciones i-esimas de cada vector corresponde
	!			  con el conjunto de caracteristicas del cuerpo i-esimo
	!
	!			+ Primero se adjudica la posicion
	!
	!			+ La velocidad se calcula en funcion de la posicion
	!			  y se le añade una magnitud aleatoria proporcional a v
	!
	!--------------------------------------------------------------------------
	!	Parametros de entrada: 
	!		- m: Vector de Masas
	!		- x, y : Vectores de posiciones
	!		- vx, vy: Vectores de las velocidades
	!		- n : Numero de cuerpos de cada tipo (rocoso o gaseoso)
	!--------------------------------------------------------------------------
	!	Parametros de salida:
	!		- x, y, vx, vy: Adjudica las posiciones y velocidades iniciales 
	!===========================================================================

	implicit none

	integer, intent (in) :: n 
	double precision,intent(inout) :: x(1:n), y(1:n), vx(1:n), vy(1:n), m(1:n), r(1:n)
	
	integer :: i
	double precision :: aux1(1:n), aux2(1:n), v(1:n), alpha 


	!-----------------------------------------
	! 		POSICIONES Y VELOCIDADES
	!-----------------------------------------

	call random_seed()
	call random_number(x)
	call random_number(y)
	call random_number(aux1)
	call random_number(aux2)
	! call random_number(vx)

	x = 15*x + 0.01 		! 0.01 < x < 15.01  Unidades Astronomicas
	y = 2*10*y - 10			! -10 < y < 10  	Unidades Astronomicas

	do i = 1, n
		v(i) = 1/(x(i)**2 + y(i)**2)**0.25 		! La raiz de la norma -- RAIZ((x^2+y^2)^1/2)
		alpha = atan(y(i)/x(i))
		vx(i) = v(i) * sin(alpha)
		if ( alpha > 0 ) then 
			vy(i) = -v(i) * cos(alpha)*(0.2*aux1(i) + 0.9)
		else
			vy(i) = v(i) * cos(alpha)*(0.2*aux2(i) + 0.9)
		end if 
	end do

	
	!------------------------
	!		MASA Y RADIO
	!------------------------

	m = 2.5E-8
	r = 8E-6

end subroutine condiciones_iniciales_rock

subroutine condiciones_iniciales_gas(m, r, x, y, vx, vy, n)

	!===========================================================================
	!
	!	Subrutina para adjudicar las condiciones iniciales.
	!
	!			+ El cojunto de posiciones i-esimas de cada vector corresponde
	!			  con el conjunto de caracteristicas del cuerpo i-esimo
	!
	!			+ Primero se adjudica la posicion
	!
	!			+ La velocidad se calcula en funcion de la posicion
	!			  y se le añade una magnitud aleatoria proporcional a v
	!
	!--------------------------------------------------------------------------
	!	Parametros de entrada: 
	!		- m: Vector de Masas
	!		- x, y : Vectores de posiciones
	!		- vx, vy: Vectores de las velocidades
	!		- n : Numero de cuerpos de cada tipo (rocoso o gaseoso)
	!--------------------------------------------------------------------------
	!	Parametros de salida:
	!		- x, y, vx, vy: Adjudica las posiciones y velocidades iniciales 
	!===========================================================================

	implicit none

	integer, intent (in) :: n 
	double precision,intent(inout) :: x(1:n), y(1:n), vx(1:n), vy(1:n), m(1:n), r(1:n)
	
	integer :: i
	double precision :: aux1(1:n), aux2(1:n), v(1:n), alpha 


	!-----------------------------------------
	! 		POSICIONES Y VELOCIDADES
	!-----------------------------------------

	call random_seed()
	call random_number(x)
	call random_number(y)
	call random_number(aux1)
	call random_number(aux2)
	! call random_number(vx)

	x = -17*x - 3 		! -3 < x < -20  	Unidades Astronomicas
	y = 2*15*y - 15		! -15 < y < 15  	Unidades Astronomicas

	do i = 1, n
		v(i) = 1/(x(i)**2 + y(i)**2)**0.25 		! La raiz de la norma -- RAIZ((x^2+y^2)^1/2)
		alpha = atan(y(i)/x(i))
		if ( alpha > 0 ) then 
			vx(i) = -v(i) * sin(alpha)
			vy(i) = v(i) * cos(alpha)*(0.2*aux1(i) + 0.9)
		else
			vx(i) = -v(i) * sin(alpha)
			vy(i) = v(i) * cos(alpha)*(0.2*aux2(i) + 0.9)
		end if 
	end do

	
	!------------------------
	!		MASA Y RADIO
	!------------------------

	! Refs: 
	! 	Masa Jupiter: 1.898E27 [Kg] --> Mj/Ms = Mj' = 9.537688E-4 kg' 
	! 	Radio Jupiter: 69.911E6 [m] -->Rj/c = Rj' 	= 4.673195E-4 m'
	!	M/r^3 = 9345475.461 

	m = 9.5E-6
	r = 1.0054E-4

end subroutine condiciones_iniciales_gas


subroutine aceleraciones(x, y, m, ax, ay, n)

	!==================================================================
	!	Subrutina que calcula la aceleracion producida por la 
	! 	interaccion gravitatoria de una masa en la posicion (x,y).
	!-----------------------------------------------------------------
	!	SOLO me interesa la interaccion de las masas con el Sol
	!	NO considero las interacciones gravitatorias entre las masas
	!-----------------------------------------------------------------
	!	Devuelve las aceleraciones dividas por componentes
	!==================================================================

	implicit none

	integer, intent(in) :: n
	double precision, intent(in) :: x(1:n), y(1:n), m(1:n)
	double precision, intent(out) :: ax(1:n), ay(1:n)

	double precision, parameter :: G = 6.67E-11

	integer :: i

	ax = 0
	ay = 0
	
	do i = 1, n
		if (m(i) /= 0) then
			ax(i) = ax(i) - (x(i))/(x(i)**2 + y(i)**2)**(3.0/2.0)
			ay(i) = ay(i) - (y(i))/(x(i)**2 + y(i)**2)**(3.0/2.0)
		end if 
	end do

end subroutine aceleraciones


subroutine vervelet_pos(m, x, y, vx, vy, ax, ay, n, h)
	implicit none

	integer, intent(in) :: n
	double precision,intent(in) :: h, ax(1:n), ay(1:n), m(1:n)
	double precision, intent(inout) :: x(1:n), y(1:n), vx(1:n), vy(1:n)

	integer :: i

	do i = 1, n
		if (m(i) /= 0) then	
			x(i) = x(i) + h*vx(i) + 0.5*h**2*ax(i)
			y(i) = y(i) + h*vy(i) + 0.5*h**2*ay(i)
		end if	
	end do

end subroutine vervelet_pos


subroutine w(m, vx, vy, ax, ay, wx, wy, n, h)
	implicit none

	integer, intent(in) :: n
	double precision, intent(in) :: h, vx(1:n), vy(1:n), ax(1:n), ay(1:n), m(1:n)
	double precision, intent(inout) :: wx(1:n), wy(1:n)
	
	integer :: i
	
	do i = 1, n
		if (m(i) /= 0) then
			wx(i) = vx(i) + 0.5*h*ax(i)
			wy(i) = vy(i) + 0.5*h*ay(i)
		end if
		end do
end subroutine w


subroutine vervelet_vel(m, vx, vy, ax, ay, wx, wy, n, h)
	implicit none

	integer, intent(in) :: n
	double precision, intent(in) :: h,wx(1:n), wy(1:n), ax(1:n), ay(1:n), m(1:n)
	double precision, intent(out) :: vx(1:n), vy(1:n)

	integer :: i 

	do i = 1, n
		if (m(i) /= 0) then
			vx(i) = wx(i) + 0.5*h*ax(i)
			vy(i) = wy(i) + 0.5*h*ay(i) 
		end if
	end do

end subroutine vervelet_vel


subroutine SunDist(m, x, y, vx, vy, ESun, n)
	implicit none
	integer, intent(in) :: n 
	double precision, intent(in) :: x(1:n), y(1:n) 
	double precision, intent(inout) :: m(1:n), vx(1:n), vy(1:n), ESun 

	integer :: i 
	real, parameter :: R_Sol2 = 2.162687122E-5 	! U.Astro 

	do i = 1, n
		! Si se acerca demasiado al Sol o se aleja demasiado
		if (((x(i)**2 + y(i)**2) <=  R_Sol2) .or. (x(i)**2 + y(i)**2) >= 80**2) then

			! Calculo la E.Mec que tiene para tenerla en cuenta en la conserv. de la Energia
			ESun = Esun + 0.5*m(i)*(vx(i)**2 + vy(i)**2) - m(i)/(x(i)**2+y(i)**2)**0.5

			! 'Elimino' la masa. Comporbar siempre que m/=0
			vx(i) = 0
			vy(i) = 0
			m(i) = 0

		end if 
	end do

	
end subroutine SunDist

!
! 	Subrutina que detecta cuando se produce una colision entre dos cuerpos.
!
! 	Se considera colision cuando la distancia entre los centros es menor o igual que la suma de los 
!	radios (opcionalmente se puede multiplicar la suma de los radios por un parametros para compensar 
!	el hecho de que no estamos considerando la interracion gravitartoria entre las masas)
!	
!	Antes de hacer la comprobacion, controlar que sea una masa valida (/= 0 ).
!
!	La velocidad de la masa resultante es la suma de las dos que colisionan.
!
subroutine colision_rock(r, x, y, vx, vy, m, Qcol, n)
	
	implicit none 
	
	integer, intent(in) :: n
	double precision, intent(in) :: x(1:n), y(1:n)
	double precision, intent(inout) :: r(1:n), vx(1:n), vy(1:n), m(1:n), Qcol

	integer :: i, j

	do i = 1, n
		if (m(i) /= 0) then 
			do j = i+1, n
				if ( ((x(i)-x(j))**2 + (y(i)-y(j))**2) <= 6*(r(i)+r(j))  ) then
					write(*,*) "Colision  producida"

					! Energia disipada en la colision 
						! Energia cinetica
						Qcol = -0.5*m(i)*m(j)/(m(i)+m(j))*(sqrt(vx(i)**2+vy(i)**2)-sqrt(vx(j)**2+vy(j)**2))**2 

						! + Energia potencial
						Qcol = Qcol + m(j)*(1/sqrt(x(j)**2+y(j)**2) - 1/sqrt(x(i)**2+y(i)**2))

					! Radio (variacion)
					r(i) = r(i)*((m(i)+m(j))/m(i))**(1.0/3.0)
					r(j) = 0 

					! Velocidades (variacion)
					vx(i) = (m(i)*vx(i) + m(j)*vx(j))/(m(i)+m(j))
					vy(i) = (m(i)*vy(i) + m(j)*vy(j))/(m(i)+m(j))

					vx(j) = 0 
					vy(j) = 0
				
					! Masa (variacion)
					m(i) = m(i)+m(j)
					m(j) = 0

				end if
			end do
		end if
	end do

	
end subroutine colision_rock


!
!	COLISIONES MASAS GAS
!
subroutine colision_gas(r, x, y, vx, vy, m, Qcol, n)
	
	implicit none 
	
	integer, intent(in) :: n
	double precision, intent(in) :: x(1:n), y(1:n)
	double precision, intent(inout) :: r(1:n), vx(1:n), vy(1:n), m(1:n), Qcol

	integer :: i, j

	do i = 1, n
		if ((m(i) /= 0) .and. (x(i)**2+y(i)**2) >= 4**2) then 
			do j = i+1, n
				if ( ((x(i)-x(j))**2 + (y(i)-y(j))**2) <= 2*(r(i)+r(j))  ) then
					write(*,*) "Colision  producida"

					! Energia disipada en la colision 
						! Energia cinetica
						Qcol = -0.5*m(i)*m(j)/(m(i)+m(j))*(sqrt(vx(i)**2+vy(i)**2)-sqrt(vx(j)**2+vy(j)**2))**2 

						! + Energia potencial
						Qcol = Qcol + m(j)*(1/sqrt(x(j)**2+y(j)**2) - 1/sqrt(x(i)**2+y(i)**2))

					! Radio (variacion)
					r(i) = r(i)*((m(i)+m(j))/m(i))**(1.0/3.0)
					r(j) = 0 

					! Velocidades (variacion)
					vx(i) = (m(i)*vx(i) + m(j)*vx(j))/(m(i)+m(j))
					vy(i) = (m(i)*vy(i) + m(j)*vy(j))/(m(i)+m(j))

					vx(j) = 0 
					vy(j) = 0
				
					! Masa (variacion)
					m(i) = m(i)+m(j)
					m(j) = 0

				end if
			end do
		end if
	end do

	
end subroutine colision_gas


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
		if (m(i)/=0) then
			T = T + 0.5*m(i)*(vx(i)**2 + vy(i)**2)
		end if 
	end do
	do i = 1, n
		if (m(i)/=0) then
			V = V - m(i)/(x(i)**2+y(i)**2)**0.5
		end if
	end do

	E_mec = T + V

	
end subroutine Ener_mec