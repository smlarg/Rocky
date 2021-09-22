program write_test
 INTEGER, DIMENSION(3,4)::c
 INTEGER :: i,j
 c = reshape((/(i,i=1,12)/),(/3,4/))
 OPEN(UNIT=9,FILE='harry.txt',POSITION='APPEND')
 DO i=LBOUND(c,2),UBOUND(c,2)
  DO j=LBOUND(c,1),UBOUND(c,1)
   WRITE(UNIT=9,FMT='(I3)', ADVANCE='no') c(j,i)
  END DO
  WRITE(UNIT=9,FMT=*) ' '
 END DO
! PRINT '(8F8.2)',((c(i,j),i=1,3),j=1,3)
end program write_test
