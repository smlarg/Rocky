PROGRAM single_particle_thermal
 IMPLICIT NONE
 REAL :: p1_in, p2_in, x1
 REAL :: p1_out, p2_out, f
 INTEGER :: i,p1_out_int,p2_out_int
 REAL, DIMENSION(200,200) :: p2p1_phasespace = RESHAPE( (/(0,i=1,200*200)/), (/200,200/) ) !I'm not sure I got the syntax right here

 p1_in = -4.0
 p2_in = -4.0
 x1 = 0
 p1_loop : DO
  WRITE (*,*) p1_in
  IF ( p1_in > 4.0 ) EXIT p1_loop
  p2_loop : DO
   IF ( p2_in > 4.0) THEN
    p2_in = -4.0
    EXIT p2_loop
   END IF
   x1_loop : DO
    IF ( x1 > 6.3) THEN
     x1 = 0
     EXIT x1_loop
    END IF
    CALL whatever(p1_in,p2_in,x1,p1_out,p2_out,f)
    p1_out_int = (p1_out + 5 + 0.001) * 10
    p2_out_int = (p2_out + 10 + 0.001) * 10
    IF (p1_out_int > 0 .and. p1_out_int < 201 .and. p2_out_int > 0 .and. p2_out_int < 201) &
     & p2p1_phasespace(p1_out_int,p2_out_int) = p2p1_phasespace(p1_out_int,p2_out_int) + f
    x1 = x1 + 0.1
   END DO x1_loop
   p2_in = p2_in + 0.1
  END DO p2_loop
  p1_in = p1_in + 0.1
 END DO p1_loop

 call write_p2p1(p2p1_phasespace)

 
CONTAINS

 SUBROUTINE whatever(p1_in,p2_in,x1,p1_out,p2_out,f)
  !just to skip where the distribution is extremely small, to save time
  IMPLICIT NONE
  REAL, INTENT(IN) :: p1_in,p2_in,x1
  REAL, INTENT(OUT) :: p1_out,p2_out,f
  
!  IF (p1_in**2 + p2_in**2 > 3.5**2) THEN
!   p1_out = p1_in
!   p2_out = p2_in
!   f = 0.0
!  ELSE
   f = 1
!   p1_out = p1_in
!   p2_out = p2_in
   CALL theotherone(p1_in,p2_in,x1,p1_out,p2_out)
!  END IF
 END SUBROUTINE whatever


 SUBROUTINE theotherone(p1_in,p2_in,x1,p1_out,p2_out)
  IMPLICIT NONE
  !Declare variables
  REAL, INTENT(IN) :: p1_in,p2_in,x1
  REAL, INTENT(OUT) :: p1_out,p2_out
  REAL :: t,dt,t_max
  REAL :: temporary_t
  REAL, DIMENSION(1:3) :: x,p,v,f
  INTEGER :: i
  
  !Executable statements
  !Initialize variables
  x = (/x1,0,0/)
  p = (/p1_in, p2_in,0/)
  t = 0.0
  dt = 0.01
  t_max =  30
  
  !Move particle
  DO
   call velocity(p,v)
   IF (x(1) > 0) THEN
    IF (v(1) < 0 .and. (t + (-x(1)/v(1)) < t_max)) THEN
     temporary_t = (-x(1)/v(1))
     t = t + temporary_t
     x = x + v*temporary_t
    ELSE
     x = x + v*(t_max - t)
     t = t_max
    END IF
   END IF
   IF ( t >= t_max) EXIT
   CALL force(x,v,t,f)
   p = p + f*dt
   x = x + v*dt
   t = t + dt
  END DO
  !set return variables, and then we're done here
  p1_out = p(1)
  p2_out = p(2)
 END SUBROUTINE theotherone


 SUBROUTINE write_p2p1(p2p1)
  REAL,  DIMENSION(:,:), INTENT(IN) :: p2p1
  INTEGER :: i,j
  WRITE(*,FMT='(F10.2)') p2p1(50,100)
  OPEN(UNIT=9,FILE='p2p1-square.txt',POSITION='APPEND')
  DO i=LBOUND(p2p1,2),UBOUND(p2p1,2)
   DO j=LBOUND(p2p1,1),UBOUND(p2p1,1)
    WRITE(UNIT=9, FMT='(F14.8)',ADVANCE='no') p2p1(j,i)
   END DO
   WRITE(UNIT=9, FMT=*) ' '
  END DO
  CLOSE(UNIT=9)
 END SUBROUTINE write_p2p1


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
  INTEGER :: i
  DO i = 1,3
   E(i) = 0
  END DO
  IF (x(1) <= 0) THEN
   E(2) = E0 * (sin(x(1)-t) -  R *sin(-x(1)-t))
  END IF
 END SUBROUTINE e_field


END PROGRAM single_particle_thermal
