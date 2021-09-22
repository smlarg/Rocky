program just_one_particle
 IMPLICIT NONE
 REAL,DIMENSION(6) :: particle
 REAL :: t, t_0, dt, t_max
 REAL :: p1, p2
 
 dt = 0.0001

 p1 = -1.0
 p2 = 1.5
 t_0 = 3.15
 
 t = t_0
 t_max = t + 50
 particle = (/0.,0.,0.,p1,p2,0./)
 
 DO
  IF ( particle(1) > 0 .and. particle(4) > 0) EXIT
  IF (t > t_max) THEN
   write(*,*) 'stuck in vacuum'
   EXIT
  END IF
  call push(particle, t, dt)
  t = t + dt
 END DO
 
 write(*,*) (/0.,0.,0.,p1,p2,0./)
 write(*,*) t_0
 write(*,*) particle
 write(*,*) t
 
CONTAINS  
  
 SUBROUTINE push(particle, time, dt)
  REAL, INTENT(IN) :: time, dt
  REAL, DIMENSION(1:6), INTENT(INOUT) :: particle
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

 SUBROUTINE b_field(x,t,B)
  IMPLICIT NONE
  REAL, DIMENSION(1:3), INTENT(IN) :: x
  REAL, INTENT(IN) :: t
  REAL, DIMENSION(1:3), INTENT(OUT) :: B
  REAL :: B0 = 6.0, R = 1.0
  !
  B = 0
  IF (x(1) <= 0) THEN
   B(3) = B0 * (sin(x(1)-t) + R *sin(-x(1)-t))
  END IF
 END SUBROUTINE b_field

 SUBROUTINE e_field(x,t,E)
  IMPLICIT NONE
  REAL, DIMENSION(1:3), INTENT(IN) :: x
  REAL, INTENT(IN) :: t
  REAL, DIMENSION(1:3), INTENT(OUT) :: E
  REAL :: E0 = 6.0, R = 1.0
  !
  E = 0
  IF (x(1) <= 0) THEN
   E(2) = E0 * (sin(x(1)-t) - R *sin(-x(1)-t))
  END IF
 END SUBROUTINE e_field

end program just_one_particle