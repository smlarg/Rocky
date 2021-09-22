program just_one_particle
 IMPLICIT NONE
 REAL,DIMENSION(3) :: p,x,v,f
 REAL :: t, t_0, dt, t_max
 REAL :: p1, p2
 INTEGER :: i
 
 dt = 0.000001

 p1 = -1.0
 p2 = 1.5
 t_0 = 3.15
 
 t = t_0
 t_max = t + 50
 p = (/p1,p2,0./)
 x = 0.0
 v = 0.0
 
 
 DO
  IF ( x(1) > 0 .and. p(1) > 0) EXIT
  IF (t > t_max) THEN
   write(*,*) 'stuck in vacuum'
   EXIT
  END IF
  CALL velocity(p,v)
  CALL force(x,v,t,f)
  DO i=1,3
   p(i) = p(i) + f(i)*dt
   x(i) = x(i) + v(i)*dt
  END DO
  
  t = t + dt
 END DO
 
 write(*,*) (/0.,0.,0.,p1,p2,0./)
 write(*,*) t_0
 write(*,*) x,p
 write(*,*) t
 
CONTAINS  
  
  SUBROUTINE velocity(p,v)
  IMPLICIT NONE
  REAL, DIMENSION(1:3), INTENT(IN)  :: p
  REAL, DIMENSION(1:3), INTENT(OUT) :: v
  REAL :: gamma_factor
  INTEGER :: i
  gamma_factor = sqrt(1+p(1)**2+p(2)**2+p(3)**2)
  DO i = 1,3
   v(i) = p(i)/gamma_factor
  END DO
 END SUBROUTINE velocity

 SUBROUTINE force(x,v,t,f)
  IMPLICIT NONE
  REAL, DIMENSION(1:3), INTENT(IN) :: x,v
  REAL, INTENT(IN) :: t
  REAL, DIMENSION(1:3), INTENT(OUT) :: f
  REAL, DIMENSION(1:3) :: E,B
  CALL e_field(x,t,E)
  CALL b_field(x,t,B)
  f(1) = E(1) + v(2)*B(3) - v(3)*B(2)
  f(2) = E(2) - v(1)*B(3) + v(3)*B(1)
  f(3) = E(3) + v(1)*B(2) - v(2)*B(1)
 END SUBROUTINE force



 SUBROUTINE b_field(x,t,B)
  IMPLICIT NONE
  REAL, DIMENSION(1:3), INTENT(IN) :: x
  REAL, INTENT(IN) :: t
  REAL, DIMENSION(1:3), INTENT(OUT) :: B
  REAL :: B0 = 6.0, R = 0.0
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
  REAL :: E0 = 6.0, R = 0.0
  !
  E = 0
  IF (x(1) <= 0) THEN
   E(2) = E0 * (sin(x(1)-t) - R *sin(-x(1)-t))
  END IF
 END SUBROUTINE e_field

end program just_one_particle