PROGRAM rocky

! use data_module
 use boris_push
! use diagnostics
 
 
 IMPLICIT NONE
 REAL, ALLOCATABLE, DIMENSION(:,:) :: particle_array
 INTEGER :: i, particle_count
 REAL, parameter :: p_range 	= 8.
 real, parameter :: dp 			= 0.1
 REAL, parameter :: t_max 		= 20.0
 real, parameter :: dt			= 0.01
 real, parameter :: x_min 		= 0.0
 real, parameter :: x_max		= x_min + 6.999
 real, parameter :: dx			= 0.1
 REAL :: p1, p2, x1, t
 !
 p1 = -p_range
 p2 = p1
 x1 = x_min
 t = 0.0
 i = 1
 particle_count =  int( (2.0*p_range/dp) * (2.0*p_range/dp)* ((x_max - x_min)/dx)*1.05 )
 ALLOCATE(particle_array(1:13,1:particle_count))
 particle_array = 0
 p1_loop : DO
  IF ( p1 > p_range ) EXIT p1_loop
  p2_loop : Do
   IF (p2 > p_range ) THEN
    p2 = -p_range
    EXIT p2_loop
   END IF
   x1_loop : DO
    IF ( x1 > x_max) THEN
     x1 = x_min
     EXIT x1_loop
    END IF
    IF (i > UBOUND(particle_array,2)) EXIT p1_loop
    particle_array(1:7,i)  = (/x1,0.0,0.0,p1,p2,0.0,dist_func((/p1,p2,0.0/))/)
    particle_array(8:13,i) = (/x1,0.0,0.0,p1,p2,0.0/)
    i = i + 1
    x1 = x1 + dx
   END DO x1_loop
   p2 = p2 + dp
  END DO p2_loop
  p1 = p1 + dp
 END DO p1_loop
 !
 t = 0.0
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
 OPEN(UNIT=9,FILE='particle_array.dat',STATUS='replace',FORM='unformatted')
 WRITE(UNIT=9) particle_array
 CLOSE(UNIT=9)

 open(unit=10,file='array_size.txt',status='replace',form='formatted')
 write(10,*) ubound(particle_array)
 close(10)

 CONTAINS

 REAL FUNCTION dist_func(p)
  IMPLICIT NONE
  REAL, DIMENSION(1:3), INTENT(IN) :: p
  !
  IF (dot_product(p,p) > p_range**2) THEN
   dist_func = 0
  ELSE
   dist_func = exp(-dot_product(p,p)/1.2)
!   dist_func = 1.0
  END IF
 END FUNCTION dist_func

END PROGRAM rocky
