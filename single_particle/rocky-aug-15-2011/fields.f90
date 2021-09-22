module fields

implicit none

contains

 SUBROUTINE e_field(Y,t,E)
  IMPLICIT NONE
  REAL, DIMENSION(1:3), INTENT(IN) :: Y
  REAL, INTENT(IN) :: t
  REAL, DIMENSION(1:3), INTENT(OUT) :: E
  REAL :: x
  REAL :: E0 = 6.0 
  !
  E = 0.
  x = Y(1)
  IF (x <= 0) THEN
   E(2) =  E0 * (sin(x-t) -  sin(-x-t))
  ELSE
   E(2) =  0.
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
  B = 0.
  x = Y(1)
  IF (x <= 0) THEN
   B(3) =  B0 * (sin(x-t) + sin(-x-t))
  ELSE
   B(3) =  0.
  END IF
 END SUBROUTINE b_field
 
end module fields