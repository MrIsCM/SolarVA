!==============================================
!
! 		Simulacion Formacion Planetoides
! 
! 				 30/06/2022
!
!==============================================
program sistema_solar
	implicit none

	! Posicion y derivadas
	integer, parameter :: n = 10			! Numero de masas orbitando inicialmente

	! Masas de tipo rocoso
	double precision :: rock_x(1:n), rock_y(1:n), rock_vx(1:n), rock_vy(1:n), rock_ax(1:n), rock_ay(1:n), rock_wx(1:n), rock_wy(1:n), rock_m(1:n)

	! Masas de tipo gaseoso
	double precision :: gas_x(1:n), gas_y(1:n), gas_vx(1:n), gas_vy(1:n), gas_ax(1:n), gas_ay(1:n), gas_wx(1:n), gas_wy(1:n), gas_m(1:n) 

	
	! Paramteros temporales
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

	f_r = c 					! m
	f_t = 5022361.106 			! s
	f_v = 29786.7869			! m/s

	t=0
	h=0.001


	open(10, file='Datos/README.txt', status='unknown')
		write(10,*) "Datos simulacion formacion S.Solar - Fortran95"
		write(10,*) "h = 		", h, "		Intervalo de tiempo (~h*58 dias)"
		write(10,*) "i = 		", iter, "		Numero de iteraciones"

	open(1, file='Datos/Rock.txt', status='unknown')
	open(2, file='Datos/Gas.txt', status='unknown')
	open(4, file='Datos/Constantes.txt', status='unknown')

	do i = 1, iter

		call Ener_mec(x, y, vx, vy, m, E_t)
		call momento_angular(x, y, vx, vy, m, L_t)

		if (mod(i,10) == 0) then 
			write(1,*) t, rock_x, rock_y
			!write(2,*) t, vx, vy
			!write(3,*) t, ax, ay
			write(4,*) t, E_t, L_t

		end if 
		
		call aceleraciones(x, y, m, ax, ay)				! Aceleraciones para las posiciones de la iter anterior
		call vervelet_pos(x, y, vx, vy, ax, ay, h)		! Calculo nuevas posiciones
		call w(vx, vy, ax, ay, wx, wy, h)				! Calculo la funcion aux w
		call aceleraciones(x, y, m, ax, ay)				! Calculo las aceleraciones a partir de las nuevas pos
		call vervelet_vel(vx, vy, ax, ay, wx, wy, h)	! Calculo las nuevas velocidades mediante las ac nuevas

		t = t + h
	end do 												! Repite

	close(1)
	close(2)
	close(4)
	
end program sistema_solar


subroutine condiciones_iniciales(rock_x, rock_y, rock_vx, rock_vy, rock_m, gas_x, gas_y, gas_vx, gas_vy, gas_m, f_r, f_v, n)
	implicit none

	integer :: n
	double precision,intent(inout) :: rock_x(1:n), rock_y(1:n), rock_vx(1:n), rock_vy(1:n), rock_m(1:n),  gas_x(1:n), gas_y(1:n), gas_vx(1:n), gas_vy(1:n), gas_m(1:n)
	double precision, intent(in) :: f_r, f_v

	integer :: i

	

end subroutine condiciones_iniciales

subroutine aceleraciones(x, y, m, ax, ay)
	implicit none

	integer, parameter :: n = 10
	double precision, intent(in) :: x(1:n), y(1:n), m(1:n)
	double precision, intent(out) :: ax(1:n), ay(1:n)

	double precision, parameter :: G = 6.67 * 10.0**(-11)

	integer :: i, j

	ax = 0
	ay = 0
	
	do i = 1, n
		do j = 1, n
			if (i /= j) then
				ax(i) = ax(i) - m(j)*(x(i)-x(j))/((x(i)-x(j))**2 + (y(i)-y(j))**2)**(3.0/2.0)
				ay(i) = ay(i) - m(j)*(y(i)-y(j))/((x(i)-x(j))**2 + (y(i)-y(j))**2)**(3.0/2.0)
			end if 
		end do
	end do

end subroutine aceleraciones

subroutine vervelet_pos(x, y, vx, vy, ax, ay, h)
	implicit none

	integer, parameter :: n=10
	double precision,intent(in) :: h, ax(1:n), ay(1:n)
	double precision, intent(inout) :: x(1:n), y(1:n), vx(1:n), vy(1:n)

	integer :: i

	do i = 1, n
		x(i) = x(i) + h*vx(i) + 0.5*h**2*ax(i)
		y(i) = y(i) + h*vy(i) + 0.5*h**2*ay(i)
	end do

end subroutine vervelet_pos


subroutine w(vx, vy, ax, ay, wx, wy, h)
	implicit none

	integer, parameter :: n = 10
	double precision, intent(in) :: h, vx(1:n), vy(1:n), ax(1:n), ay(1:n)
	double precision, intent(inout) :: wx(1:n), wy(1:n)
	
	integer :: i
	
	do i = 1, n
		wx(i) = vx(i) + 0.5*h*ax(i)
		wy(i) = vy(i) + 0.5*h*ay(i)
	end do
end subroutine w

subroutine vervelet_vel(vx, vy, ax, ay, wx, wy, h)
	implicit none

	integer, parameter :: n = 10
	double precision, intent(in) :: h,wx(1:n), wy(1:n), ax(1:n), ay(1:n)
	double precision, intent(out) :: vx(1:n), vy(1:n)

	integer :: i 

	do i = 1, n
		vx(i) = wx(i) + 0.5*h*ax(i)
		vy(i) = wy(i) + 0.5*h*ay(i) 
	end do

end subroutine vervelet_vel

subroutine Ener_mec(x, y, vx, vy, m, E_mec)
	implicit none

	integer, parameter :: n = 10
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
		do j = 1, n
			if (i /= j) then
				V = V - m(i)*m(j)/((x(i)-x(j))**2+(y(i)-y(j))**2)**0.5
			end if
		end do
	end do

	E_mec = T + 0.5*V			! La mitad de E. Potencial ya que la E.Pot de i,j es = a j,i

	
end subroutine Ener_mec

subroutine momento_angular(x, y, vx, vy, m, L)
	implicit none
	double precision,intent(in) :: x(1:10), y(1:10), vx(1:10), vy(1:10), m(1:10)
	double precision,intent(out) ::  L

	integer, parameter :: n = 10
	integer :: i
	L = 0

	do i = 1, n
		L = L + x(i)*m(i)*vy(i) - y(i)*m(i)*vx(i)
	end do

end subroutine momento_angular

subroutine periodo(y, t, period_orbital)
	implicit none
	double precision ,intent(in) :: y, t 
	double precision ,intent(out) :: period_orbital
end subroutine periodo