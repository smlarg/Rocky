program raw_reader

 use diagnostics
 
 implicit none
 
 real, allocatable, dimension(:,:) :: particle_array
 integer :: i
 
 allocate(particle_array(1:13,1:1881331))
 
 open(unit=9,file='particle_array.dat',status='old',form='unformatted',action='read')
 read(unit=9) particle_array
 close(unit=9)
 
 call write_dist_energy(particle_array)
 
! open(unit=10,file='scatter_plot.txt',status='replace',form='formatted',action='write')
! 
! do i =1, ubound(particle_array,2)
!  !if(particle_array(4,i) .gt. 20. .and. particle_array(4,i) .lt. 25.) then
!  if(particle_array(4,i) .gt. 15.) then
!   write(10,*) (/energy(particle_array(4:6,i)),particle_array(11,i)/)
!  end if
! end do
! 
! close(10)
 
  
 
! contains
! 
! real function energy(p)
!  implicit none
!  REAL, DIMENSION(1:3), INTENT(IN)  :: p
!  REAL :: gamma_factor
!  !
!  gamma_factor = sqrt(1+p(1)**2+p(2)**2+p(3)**2)
!  energy = gamma_factor - 1
! END FUNCTION energy
 
 
 
 
end program raw_reader
