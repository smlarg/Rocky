PROGRAM rocky

 use data_module_test
 use boris_push
 
 
 IMPLICIT NONE
 REAL, ALLOCATABLE, DIMENSION(:,:) :: particle_array
 INTEGER :: i, particle_count
! REAL, parameter :: p_range = -5.0
 REAL :: p1, p2, x1, pmax , xmin, xmax
 REAL :: t, dt, t_max, dp, dx!, particle_count
 !
 p1 = -p_range
 p2 = p1
 pmax = p_range
 x1 = 6.3
 xmin = x1
 xmax = xmin + 6.999
 dp = 0.15
 dx = 0.1
 i = 1
 particle_count =  int( (2.0*pmax/dp) * (2.0*pmax/dp)* ((xmax - xmin)/dx)*1.05 )
 ALLOCATE(particle_array(1:7,1:particle_count))
 particle_array = 0
 p1_loop : DO
  IF ( p1 > pmax ) EXIT p1_loop
  p2_loop : Do
   IF (p2 > pmax ) THEN
    p2 = -pmax
    EXIT p2_loop
   END IF
   x1_loop : DO
    IF ( x1 > xmax) THEN
     x1 = xmin
     EXIT x1_loop
    END IF
    IF (i > UBOUND(particle_array,2)) EXIT p1_loop
    particle_array(1:7,i) = (/x1,0.0,0.0,p1,p2,0.0,dist_func((/p1,p2,0.0/))/)
    i = i + 1
    x1 = x1 + dx
   END DO x1_loop
   p2 = p2 + dp
  END DO p2_loop
  p1 = p1 + dp
 END DO p1_loop
 !
 t = 0.0
 dt = 0.01
 t_max = 40.0
 DO
  if (t>t_max) EXIT
  DO i = 1, UBOUND(particle_array,2)
   if ( particle_array(7,i) ==  0 ) cycle
   call push(particle_array(1:7,i),t,dt)
  END DO
  write (*,*) t
  t =  t + dt
 END DO
 !
 call print_out_p2p1(particle_array)

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


 SUBROUTINE print_out_p2p1(particle_array)
  REAL, DIMENSION(:,:), INTENT(IN) :: particle_array
  REAL, DIMENSION(200,200) :: p2p1
  INTEGER :: p2_int, p1_int
  !
  p2p1 = 0
  do i = 1, UBOUND(particle_array,2)
   p2_int = (particle_array(5,i) + 10  + 0.001) * 10
   p1_int = (particle_array(4,i) + 5  + 0.001) * 10
   IF (p2_int > 0 .and. p2_int < 201  &
    & .and. p1_int > 0 .and. p1_int < 201) &
    & p2p1(p1_int,p2_int) = p2p1(p1_int,p2_int) + particle_array(7,i)
  END do
  !
  OPEN(UNIT=9,FILE='p2p1.dat',STATUS='replace',FORM='unformatted')
  WRITE(UNIT=9) p2p1
  CLOSE(UNIT=9)
 END SUBROUTINE print_out_p2p1


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


 REAL FUNCTION dist_func(p)
  IMPLICIT NONE
  REAL, DIMENSION(1:3), INTENT(IN) :: p
  !
  IF (dot_product(p,p) > 5.0**2) THEN
   dist_func = 0
  ELSE
   dist_func = exp(-dot_product(p,p)/2.5)
!   dist_func = 1.0
  END IF
 END FUNCTION dist_func

END PROGRAM rocky
