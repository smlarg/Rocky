PROGRAM file_name_test
 IMPLICIT NONE
 REAL :: real_number
 INTEGER :: i
 CHARACTER(len=5) :: nuts
 CHARACTER(LEN=20) :: name
 
 real_number = 10.01
 i = 325
 WRITE(nuts,'(F5.2)') real_number
 
! WRITE(*,'(F5.2)') real_number
! WRITE(*,*) nuts



 name = 'number' // nuts // '.dat'
 
 OPEN(UNIT=9,FILE=name,STATUS='replace',FORM='unformatted')
 WRITE(UNIT=9) real_number
 CLOSE(UNIT=9)


END PROGRAM file_name_test