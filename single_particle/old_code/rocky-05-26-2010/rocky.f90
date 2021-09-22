PROGRAM rocky

 use data_module
 use boris_push
 use diagnostics
 
 
 IMPLICIT NONE
 REAL, ALLOCATABLE, DIMENSION(:,:) :: particle_array
 INTEGER :: i, particle_count
! REAL, parameter :: p_range = -5.0
 REAL :: p1, p2, x1, pmax , xmin, xmax
 REAL :: t, dt, dp, dx!, particle_count
 !
 p1 = -p_range
 p2 = p1
 pmax = p_range
 x1 = 6.3
 xmin = x1
 xmax = xmin + 6.999
 dp = 0.1
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
 call write_dist(particle_array)

 CONTAINS

 REAL FUNCTION dist_func(p)
  IMPLICIT NONE
  REAL, DIMENSION(1:3), INTENT(IN) :: p
  !
  IF (dot_product(p,p) > 3.5**2) THEN
   dist_func = 0
  ELSE
   dist_func = exp(-dot_product(p,p)/1.2)
!   dist_func = 1.0
  END IF
 END FUNCTION dist_func

END PROGRAM rocky
