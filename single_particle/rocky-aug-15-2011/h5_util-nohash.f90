
module hdf5_util

use hdf5 
   
implicit none

!include 'mpif.h'

private

integer, parameter :: p_double = kind(1.0d0)
integer, parameter :: p_single = kind(1.0e0)

integer, parameter :: p_lower = 1
integer, parameter :: p_upper = 2

type :: t_h5_tune
  
  ! Use gpfs hints / optimizations
  logical :: gpfs = .false.

  ! use MPI+POSIX instead of MPI-IO for parallel I/O
  logical :: posix = .false.                
  
  ! metadata block size
  integer(hsize_t) :: meta_block_size = -1
  
  ! cache size (I think this is only for read operations...)
  integer(size_t)  :: cache = -1
  
  ! data alignment in file
  integer(hsize_t), dimension(2) :: alignment = -1
  
  ! sieve buffer size
  integer(size_t) :: sieve_buf_size = -1 
  
  ! MPI-IO hints
  character( len = 1024 ) :: access_style = '-'
  logical :: collective_buffering = .false.
  integer :: cb_block_size = -1
  integer :: cb_buffer_size = -1



  ! dataset tune
  
  ! Internal buffer size
  integer(hsize_t) :: buffer_size = -1
  
  ! Use chunked dataset
  logical :: chunked = .true.
  logical :: no_chunk_last_dim = .false.
  
  ! use independent or collective transfers
  integer :: transferMode 

end type t_h5_tune

type(t_h5_tune), save :: h5_tune

interface add_h5_atribute
  module procedure add_h5_atribute_str
  module procedure add_h5_atribute_str_v1
  module procedure add_h5_atribute_logical
  module procedure add_h5_atribute_logical_v1

  module procedure add_h5_atribute_single
  module procedure add_h5_atribute_v1_single

  module procedure add_h5_atribute_double
  module procedure add_h5_atribute_v1_double

  module procedure add_h5_atribute_int
  module procedure add_h5_atribute_v1_int

end interface

interface add_h5_dataset
  module procedure add_h5_dataset_1d_single
  module procedure add_h5_dataset_2d_single
  module procedure add_h5_dataset_3d_single
  
  module procedure add_h5_dataset_1d_parallel_single
  module procedure add_h5_dataset_2d_parallel_single
  module procedure add_h5_dataset_3d_parallel_single

  module procedure add_h5_dataset_1d_double
  module procedure add_h5_dataset_2d_double
  module procedure add_h5_dataset_3d_double
  
  module procedure add_h5_dataset_1d_parallel_double
  module procedure add_h5_dataset_2d_parallel_double
  module procedure add_h5_dataset_3d_parallel_double

  module procedure add_h5_dataset_1d_int
  module procedure add_h5_dataset_2d_int
  module procedure add_h5_dataset_3d_int
  
  module procedure add_h5_dataset_1d_parallel_int
  module procedure add_h5_dataset_2d_parallel_int
  module procedure add_h5_dataset_3d_parallel_int
end interface

interface open_hdf5
  module procedure open_hdf5
end interface

interface close_hdf5
  module procedure close_hdf5
end interface

interface create_hdf5_file
  module procedure create_hdf5_file
end interface

interface close_hdf5_file
  module procedure close_hdf5_file
end interface

public :: add_h5_atribute, add_h5_dataset
public :: open_hdf5, close_hdf5
public :: create_hdf5_file, close_hdf5_file
public :: t_h5_tune, h5_tune

contains

!---------------------------------------------------
subroutine open_hdf5( ierr )
!---------------------------------------------------
!
!---------------------------------------------------
   implicit none
   
   integer, intent(out) :: ierr
   
   call h5open_f(ierr) 

   ! this needs to be done after h5open_f
   h5_tune%transferMode = H5FD_MPIO_COLLECTIVE_F

end subroutine open_hdf5
!---------------------------------------------------

!---------------------------------------------------
subroutine close_hdf5( ierr )
!---------------------------------------------------
!
!---------------------------------------------------
   implicit none
   
   integer, intent(out) :: ierr
   
   call h5close_f(ierr) 

end subroutine close_hdf5
!---------------------------------------------------

!---------------------------------------------------
subroutine create_hdf5_file( filename, file_id, comm )
!---------------------------------------------------
   
   implicit none
   
   ! dummy variables
   
   character(len=*),  intent(in)  :: filename
   integer(hid_t),    intent(out) :: file_id
   integer,           intent(in) :: comm
   
   ! local variables
   
   integer :: ierr
   
   integer(hid_t) :: plist_id
   integer(size_t) :: rdcc_nelmts
   integer :: fileInfo


   
   ! executable statements

DEBUGMSG( "Creating file ", trim(filename))

   fileInfo = MPI_INFO_NULL
   call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, ierr)
    
    
   ! default for sieve size is 64 kB
   if ( h5_tune%sieve_buf_size /= -1 ) then
     call h5pset_sieve_buf_size_f( plist_id, h5_tune%sieve_buf_size, ierr )
   endif
   
   ! Metadata block size
   if ( h5_tune%meta_block_size /= -1 ) then
      call h5pset_meta_block_size_f( plist_id, h5_tune%meta_block_size, ierr )
   endif
   
   ! Cache
   if ( h5_tune%cache /= -1 ) then
      rdcc_nelmts = 1
      call h5pset_cache_f( plist_id, 0, rdcc_nelmts, h5_tune%cache, 0.75, ierr )
   endif
   
   ! Create the file
   call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, ierr, access_prp = plist_id)

   ! Close the property list
   call h5pclose_f(plist_id, ierr)
   
   ! free fileInfo if necessary
   if ( fileInfo /= MPI_INFO_NULL ) then
      call MPI_INFO_FREE( fileInfo, ierr )
   endif

end subroutine create_hdf5_file
!---------------------------------------------------

!---------------------------------------------------
subroutine close_hdf5_file( file_id )
!---------------------------------------------------

   implicit none
   
   ! dummy variables
   
   integer(HID_T),    intent(in) :: file_id
   
   ! local variables
   
   integer :: ierr
	   
   ! executable statements
   call h5fclose_f(file_id, ierr)

end subroutine close_hdf5_file
!---------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Logical and string attributes need special treatment. Remaining interfaces are generated through
! template functions
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine add_h5_atribute_str( objID, name, attribute )
!---------------------------------------------------------------------------------------------------

  implicit none
  
  integer(hid_t), intent(in) :: objID
  character( len = * ), intent(in) :: name
  character( len = * ), intent(in) :: attribute
  
  integer(hid_t) :: dataspaceID, typeID, attrID
  integer(hsize_t), dimension(1) :: dims
  integer(size_t) :: size
  
  integer :: ierr
  
  dims(1) = 1
  call h5screate_simple_f(1, dims, dataspaceID, ierr ) 
  call h5tcopy_f(H5T_NATIVE_CHARACTER, typeID, ierr)
  
  size = len(attribute)
  call h5tset_size_f(typeID, size, ierr)
  
  call h5acreate_f( objID, name, typeID, dataspaceID, attrID, ierr )
  call h5awrite_f( attrID, typeID, attribute, dims, ierr)
  call h5aclose_f( attrID, ierr )
  call h5tclose_f( typeID, ierr )
  call h5sclose_f( dataspaceID, ierr )

end subroutine add_h5_atribute_str
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine add_h5_atribute_str_v1( objID, name, attribute )
!---------------------------------------------------------------------------------------------------

  implicit none
  
  integer(hid_t), intent(in) :: objID
  character( len = * ), intent(in) :: name
  character( len = * ), dimension(:), intent(in) :: attribute
  
  integer(hid_t) :: dataspaceID, typeID, attrID
  integer(hsize_t), dimension(1) :: dims
  integer(size_t) :: maxlen
  integer :: i, ierr
  
  dims(1) = size(attribute)
  call h5screate_simple_f(1, dims, dataspaceID, ierr ) 
  
  maxlen = 0
  do i = 1, size(attribute)-1
    if (len(attribute(i)) > maxlen) maxlen = len(attribute(i))
  enddo
  
  call h5tcopy_f(H5T_NATIVE_CHARACTER, typeID, ierr)
  call h5tset_size_f(typeID, maxlen, ierr)
 
  call h5acreate_f( objID, name, typeID, dataspaceID, attrID, ierr )
  call h5awrite_f( attrID, typeID, attribute, dims, ierr)
  call h5aclose_f( attrID, ierr )
  call h5tclose_f( typeID, ierr )
  call h5sclose_f( dataspaceID, ierr )

end subroutine add_h5_atribute_str_v1
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine add_h5_atribute_logical( objID, name, attribute )
!---------------------------------------------------------------------------------------------------

  implicit none
  
  integer(hid_t), intent(in) :: objID
  character( len = * ), intent(in) :: name
  logical, intent(in) :: attribute
  
  integer(hid_t) :: dataspaceID, attrID
  integer(hsize_t), dimension(1) :: dims
  integer :: bool, ierr
  
  dims(1) = 1
  if ( attribute ) then
    bool = 1
  else
    bool = 0
  endif
  call h5screate_simple_f(1, dims, dataspaceID, ierr ) 
  call h5acreate_f( objID, name, H5T_NATIVE_INTEGER, dataspaceID, attrID, ierr )
  call h5awrite_f( attrID, H5T_NATIVE_INTEGER, bool, dims, ierr)
  call h5aclose_f( attrID, ierr )
  call h5sclose_f( dataspaceID, ierr )

end subroutine add_h5_atribute_logical
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine add_h5_atribute_logical_v1( objID, name, attribute )
!---------------------------------------------------------------------------------------------------

  implicit none
  
  integer(hid_t), intent(in) :: objID
  character( len = * ), intent(in) :: name
  logical, dimension(:), intent(in) :: attribute
  
  integer(hid_t) :: dataspaceID, attrID
  integer :: i, ierr
  integer(hsize_t), dimension(1) :: dims
  integer, dimension(:), allocatable :: bool
  
  dims(1) = size(attribute)
  allocate( bool(dims(1)) )
  do i = 1, dims(1)
	if ( attribute(i) ) then
	  bool(i) = 1
	else
	  bool(i) = 0
	endif
  enddo
  
  call h5screate_simple_f(1, dims, dataspaceID, ierr ) 
  call h5acreate_f( objID, name, H5T_NATIVE_INTEGER, dataspaceID, attrID, ierr )
  call h5awrite_f( attrID, H5T_NATIVE_INTEGER, bool, dims, ierr)
  call h5aclose_f( attrID, ierr )
  call h5sclose_f( dataspaceID, ierr )

  deallocate(bool)
  
end subroutine add_h5_atribute_logical_v1
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
!  Generate specific template functions for single, double and integer datatypes.
!  Note that the module interfaces are not generated automatically and must be explicity written
!  in the module header above.
!
!  - add_h5_atribute
!  - add_h5_atribute_v1
!  - add_h5_dataset_1d
!  - add_h5_dataset_2d
!  - add_h5_dataset_3d
!  - add_h5_dataset_1d_parallel  
!  - add_h5_dataset_2d_parallel
!  - add_h5_dataset_3d_parallel
!---------------------------------------------------------------------------------------------------


end module hdf5_util


!---------------------------------------------------------------------------------------------------
! end of hdf5_util module
!---------------------------------------------------------------------------------------------------
