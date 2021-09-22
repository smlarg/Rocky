module diagnostics

implicit none

contains



  SUBROUTINE write_dist(particle_array)
  REAL, DIMENSION(:,:), INTENT(IN) :: particle_array
  REAL, DIMENSION(200) :: dist, dist_fuzz
  INTEGER :: p1_int, i
  real :: smallest_none_zero = 1
  real, parameter :: p_dump_range = 30.
  !
  dist = 0
  dist_fuzz = 0
  do i = 1, UBOUND(particle_array,2)
   p1_int = (particle_array(4,i) + p_dump_range  + 0.001) &
   & / ( p_dump_range*2/ (ubound(dist,1)-1) )
   IF (p1_int > 0 .and. p1_int .le. ubound(dist,1) )  &
    & dist(p1_int) = dist(p1_int) + particle_array(7,i)
  END do
  do i = 2, ubound(dist_fuzz,1) - 1
   dist_fuzz(i) = (dist(i-1) + dist(i) + dist(i+1))/3.0
  end do
  dist =  dist_fuzz
  do i = 1, ubound(dist,1)
   if ( (dist(i) .gt. 0) .and. (dist(i) .lt. smallest_none_zero) ) &
   & smallest_none_zero = dist(i)
  end do
  do i = 1, ubound(dist,1)
   if(dist(i) .le. 0) dist(i) = smallest_none_zero
  end do
  do i = 1, ubound(dist,1)
   dist(i) = log(dist(i))
  end do
  !
  OPEN(UNIT=9,FILE='dist.txt',STATUS='replace',action='write',form='formatted')
  100 format (f10.5)
  WRITE(UNIT=9,fmt=100) dist
  CLOSE(UNIT=9)
  
 END SUBROUTINE write_dist



  SUBROUTINE print_out_p1x1(particle_array)
  REAL, DIMENSION(:,:), INTENT(IN) :: particle_array
  REAL, DIMENSION(200,200) :: p1x1
  INTEGER :: p1_int, x1_int, i
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
  INTEGER :: p2_int, p1_int, i
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
  INTEGER :: p1_int, p2_int, x1_int, i
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
 
 
end module diagnostics