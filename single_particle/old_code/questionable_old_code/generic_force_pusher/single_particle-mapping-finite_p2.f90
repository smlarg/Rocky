PROGRAM single_particle
 IMPLICIT NONE
 REAL :: p1, t0, pmax
 pmax = -4.5
 p1 = -0.1
 t0 = 0
 DO
  IF (p1 < pmax) THEN
   WRITE (*,*) p1
   EXIT
  END IF
  DO
   IF (t0 > 6.30) EXIT
   CALL orbiter(p1,t0)
   t0 = t0 + .1
  END DO
  t0 = 0
  p1 = p1 - .1
 END DO


CONTAINS

 SUBROUTINE orbiter(p1,t0)
  IMPLICIT NONE
  !Declare variables
  REAL, INTENT(IN) :: p1,t0
  REAL :: t,dt
  REAL, DIMENSION(1:3) :: x,p,v,f
  INTEGER :: i
  
  !Executable statements
  !Initialize variables
  x(1) = 0
  x(2) = 0
  x(3) = 0
  p(1) = p1
  p(2) = 2
  p(3) = 0
  t = t0
  dt = .001
  !Move particle
  DO
   IF (t > 10) THEN
    CALL write_p1_mapping(p1,t0,p(1),x(1))
    EXIT
   END IF
   CALL velocity(p,v)
   CALL force(x,v,t,f)
   DO i=1,3
    p(i) = p(i) + f(i)*dt
    x(i) = x(i) + v(i)*dt
   END DO
   t =  t + dt
  END DO
 END SUBROUTINE orbiter

 SUBROUTINE write_p1_mapping(p1_in,t_in,p1_fin,x_fin)
  REAL, INTENT(IN) :: p1_in,t_in,p1_fin,x_fin
  OPEN(UNIT=9,FILE='p2-2.txt',POSITION='APPEND')
  WRITE(UNIT=9, FMT=*) p1_in,t_in,p1_fin,x_fin
  CLOSE(UNIT=9)
 END SUBROUTINE write_p1_mapping
 
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
  REAL :: B0 = 6.0, R = 1.0
  INTEGER :: i
  DO i = 1,3
   B(i) = 0
  END DO
  IF (x(1) < 0) THEN
   B(3) = B0 * (sin(x(1)-t) + R *sin(-x(1)-t))
  END IF
 END SUBROUTINE b_field

 SUBROUTINE e_field(x,t,E)
  IMPLICIT NONE
  REAL, DIMENSION(1:3), INTENT(IN) :: x
  REAL, INTENT(IN) :: t
  REAL, DIMENSION(1:3), INTENT(OUT) :: E
  REAL :: E0 = 6.0, R = 1.0
  INTEGER :: i
  DO i = 1,3
   E(i) = 0
  END DO
  IF (x(1) < 0) THEN
   E(2) = E0 * (sin(x(1)-t) -  R *sin(-x(1)-t))
  END IF
 END SUBROUTINE e_field

 SUBROUTINE write_phase_space(p,x,t)
  IMPLICIT NONE
  REAL, DIMENSION(1:3), INTENT(IN) :: p,x
  REAL, INTENT(IN) :: t
  OPEN(UNIT=9,FILE="tom.txt",POSITION="APPEND")
  WRITE(UNIT=9,FMT=*) p,x,t
  CLOSE(UNIT=9)
 END SUBROUTINE write_phase_space

END PROGRAM single_particle
