module boris_push

implicit none

contains

 SUBROUTINE push(particle, time, dt)
  REAL, INTENT(IN) :: time, dt
  REAL, DIMENSION(1:7), INTENT(INOUT) :: particle
  REAL, DIMENSION(1:3) :: x, p, u_minus, E, t, B, s, u_plus
  INTEGER :: i
  ! 
  DO i=1,3
   x(i) = particle(i)
   p(i) = particle(i+3)
  END DO
  CALL e_field(x,time,E)
  u_minus = p + E*dt*0.5
  CALL b_field(x,time,B)
  t = B*dt*0.5*inv_gamma(u_minus)
  s = 2*t/(1 + dot_product(t,t))
  call cross_product(u_minus, s, u_plus)
  u_plus = u_plus + u_minus + t*dot_product(u_minus,s)
  p = u_plus + E*dt*0.5
  x = x + p*dt*inv_gamma(p)
  do i=1,3
   particle(i) = x(i)
   particle(i+3) = p(i)
  end do
 END SUBROUTINE push

 SUBROUTINE cross_product(A,B,C)
  REAL, DIMENSION(1:3), INTENT(IN) :: A,B
  REAL, DIMENSION(1:3), INTENT(OUT) :: C
  !AxB = C
  C(1) =  A(2)*B(3) - A(3)*B(2)
  C(2) = -A(1)*B(3) + A(3)*B(1)
  C(3) =  A(1)*B(2) - A(2)*B(1)
 END SUBROUTINE cross_product

 REAL FUNCTION inv_gamma(p)
  IMPLICIT NONE
  REAL, DIMENSION(1:3), INTENT(IN)  :: p
  REAL :: gamma_factor
  !
  gamma_factor = sqrt(1+p(1)**2+p(2)**2+p(3)**2)
  inv_gamma = 1/gamma_factor
 END FUNCTION inv_gamma

 SUBROUTINE e_field(Y,t,E)
  IMPLICIT NONE
  REAL, DIMENSION(1:3), INTENT(IN) :: Y
  REAL, INTENT(IN) :: t
  REAL, DIMENSION(1:3), INTENT(OUT) :: E
  REAL :: x
  REAL :: E0 = 6.0 
  !
  E = 0
  x = Y(1)
  IF (x <= 0) THEN
   E(2) =  E0 * ( 8./13.*cos(t+x) + sin(t-x) - 1./13.*sin(t+x))
  ELSE
   E(2) =  E0 * ( 8./13.*exp(-x)*cos(t-x/2) + 12./13.*exp(-x)*sin(t-x/2) )
  END IF
 END SUBROUTINE e_field
 
 SUBROUTINE b_field(Y,t,B)
  IMPLICIT NONE
  REAL, DIMENSION(1:3), INTENT(IN) :: Y
  REAL, INTENT(IN) :: t
  REAL, DIMENSION(1:3), INTENT(OUT) :: B
  REAL :: B0 = 6.0 
  REAL :: x
  !
  B = 0
  x = Y(1)
  IF (x <= 0) THEN
   B(3) =  B0 * ( -8./13.*cos(t+x) + sin(t-x) + 1./13.*sin(t+x) )
  ELSE
   B(3) =  B0 * ( -8./13.*exp(-x)*cos(t-x/2) + 14./13.*exp(-x)*sin(t-x/2) )
  END IF
 END SUBROUTINE b_field

end module boris_push