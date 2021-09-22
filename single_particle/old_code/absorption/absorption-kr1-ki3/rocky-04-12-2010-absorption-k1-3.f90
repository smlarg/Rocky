
PROGRAM boris_push
 IMPLICIT NONE
 REAL, ALLOCATABLE, DIMENSION(:,:) :: particle_array
 INTEGER :: i
 REAL :: p1, p2, p3, x1
 REAL :: t, dt, t_max
 !
 ALLOCATE(particle_array(1:7,1:300000))
 particle_array = 0
 p1 = -3.45
 p2 = -3.45
 p3 = -3.45
 x1 = 0
 i = 1
 p1_loop : DO
  IF ( p1 > 3.45 ) EXIT p1_loop
  p2_loop : Do
   IF ( p2 > 3.45 ) THEN
    p2 = -3.45
    EXIT p2_loop
   END IF
   p3_loop : DO
    IF ( p3 > 3.45 ) THEN
     p3 = -3.45
     EXIT p3_loop
    END IF
    x1_loop : DO
     IF ( x1 > 6.899) THEN
      x1 = 0.0
      EXIT x1_loop
     END IF
     IF (i > 300000) EXIT p1_loop
     particle_array(1:7,i) = (/x1,0.0,0.0,p1,p2,p3,dist_func((/p1,p2,p3/))/)
     i = i + 1
     x1 = x1 + 0.3
    END DO x1_loop
    p3 = p3 + 0.3
   END DO p3_loop
   p2 = p2 + 0.3
  END DO p2_loop
  p1 = p1 + 0.3
 END DO p1_loop
 !
 t = 0.0
 dt = 0.01
 t_max = 10.0
 DO
  if (t>t_max) EXIT
  DO i = 1, 300000
   call push(particle_array(1:7,i),t,dt)
  END DO
  write (*,*) t
  t =  t + dt
 END DO
 !
 call print_out_p1x1(particle_array)

 CONTAINS

 SUBROUTINE print_out_p1x1(particle_array)
  REAL, DIMENSION(:,:), INTENT(IN) :: particle_array
  REAL, DIMENSION(200,200) :: p1x1
  INTEGER :: p1_int, x1_int
  !
  p1x1 = 0
  do i = 1, UBOUND(particle_array,2)
   p1_int = (particle_array(4,i) + 5  + 0.001) * 10
   x1_int = (particle_array(1,i) + 5  + 0.001) * 10
   IF (p1_int > 0 .and. p1_int < 201  &
    & .and. x1_int > 0 .and. x1_int < 201) &
    & p1x1(x1_int,p1_int) = p1x1(x1_int,p1_int) + particle_array(7,i)
  END do
  !
  OPEN(UNIT=9,FILE='p1x1.dat',STATUS='replace',FORM='unformatted')
  WRITE(UNIT=9) p1x1
  CLOSE(UNIT=9)
 END SUBROUTINE print_out_p1x1

 SUBROUTINE print_out_p2p1x1(particle_array)
  REAL, DIMENSION(:,:), INTENT(IN) :: particle_array
  REAL, DIMENSION(200,200,300) :: p2p1x1
  INTEGER :: p1_int, p2_int, x1_int
  !
  p2p1x1 = 0
  do i = 1, UBOUND(particle_array,2)
   p1_int = (particle_array(4,i) + 5  + 0.001) * 10
   p2_int = (particle_array(5,i) + 10 + 0.001) * 10
   x1_int = (particle_array(1,i) + 5  + 0.001) * 10
   IF (p1_int > 0 .and. p1_int < 201 .and. p2_int > 0 .and. p2_int < 201 &
    & .and. x1_int > 0 .and. x1_int < 301) &
    & p2p1x1(p1_int,p2_int,x1_int) = p2p1x1(p1_int,p2_int,x1_int) + particle_array(7,i)
  END do
  !
  OPEN(UNIT=9,FILE='p2p1x1.dat',STATUS='replace',FORM='unformatted')
  WRITE(UNIT=9) p2p1x1
  CLOSE(UNIT=9)
 END SUBROUTINE print_out_p2p1x1

 SUBROUTINE write_p2p1(p2p1)
  REAL,  DIMENSION(:,:), INTENT(IN) :: p2p1
  INTEGER :: i,j
  OPEN(UNIT=9,FILE='p2p1.dat',STATUS='replace',FORM='unformatted')
  WRITE(UNIT=9) p2p1
  CLOSE(UNIT=9)
 END SUBROUTINE write_p2p1

 REAL FUNCTION dist_func(p)
  IMPLICIT NONE
  REAL, DIMENSION(1:3), INTENT(IN) :: p
  !
  IF (dot_product(p,p) > 3.5**2) THEN
   dist_func = 0
  ELSE
   dist_func = exp(-dot_product(p,p)/1.2)
  END IF
 END FUNCTION dist_func

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
   E(2) =  E0 * ( 6./13.*cos(t+x) + sin(t-x) - 9./13.*sin(t+x))
  ELSE
   E(2) =  E0 * ( 6./13.*exp(-3*x)*cos(t-x) + 4./13.*exp(-3*x)*sin(t-x) )
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
   B(3) =  B0 * ( -6./13.*cos(t+x) + sin(t-x) + 9./13.*sin(t+x) )
  ELSE
   B(3) =  B0 * ( -6./13.*exp(-3*x)*cos(t-x) + 22./13.*exp(-3*x)*sin(t-x) )
  END IF
 END SUBROUTINE b_field

END PROGRAM boris_push
