
subroutine twod_ps( ps_n1, ps_n2, x1min, x1max, x2min, x2max, np, data1, data2, w, q, ps )
  implicit none
  integer, intent(in) :: ps_n1, ps_n2
  real(8), intent(in) :: x1min, x1max, x2min, x2max
  integer, intent(in) :: np
  real(8), intent(in), dimension(np) :: data1, data2, w, q
  real(8), intent(out), dimension(ps_n1,ps_n2) :: ps
  !
  integer :: i, i1, i2
  real(8) :: rdx1, rdx2, lq, nx1, nx2, eps1, eps2
  !
  rdx1 = ps_n1/(x1max - x1min)
  rdx2 = ps_n2/(x2max - x2min)
  
  ps = 0.
  
  do i = 1, np
    nx1 = ( data1(i) - x1min )*rdx1 + 1
    nx2 = ( data2(i) - x2min )*rdx2 + 1
    
    i1 = nint(nx1)
    i2 = nint(nx2)
    
    eps1 = nx1 - i1 + 0.5
    eps2 = nx2 - i2 + 0.5
    
    lq = w(i)*q(i)
    
    if ( i1 > 0 .and. i1 .le. ps_n1 .and. i2 > 0 .and. i2 .le. ps_n2 ) &
      ps(i1, i2) = ps(i1, i2) + (1-eps1)*(1-eps2)*lq
    i1 = i1 + 1
    if ( i1 > 0 .and. i1 .le. ps_n1 .and. i2 > 0 .and. i2 .le. ps_n2 ) &
      ps(i1, i2) = ps(i1, i2) + (eps1)*(1-eps2)*lq
    i1 = i1 - 1
    i2 = i2 + 1
    if ( i1 > 0 .and. i1 .le. ps_n1 .and. i2 > 0 .and. i2 .le. ps_n2 ) &
      ps(i1, i2) = ps(i1, i2) + (1-eps1)*(eps2)*lq
    i1 = i1 + 1
    if ( i1 > 0 .and. i1 .le. ps_n1 .and. i2 > 0 .and. i2 .le. ps_n2 ) &
      ps(i1, i2) = ps(i1, i2) + (eps1)*(eps2)*lq
    
  enddo
  
end subroutine twod_ps

 subroutine get_block_of_fields( block_size, X, t, E, B )
  implicit none
  integer, intent(in) :: block_size
  real(8), intent(in), dimension(block_size, 3) :: X
  real(8), intent(in) :: t
  real(8), intent(out), dimension(block_size, 3) :: E, B
  !
  integer :: i
  !
  do i = 1, block_size
    call e_field( X(i,:), t, E(i,:) )
    call b_field( X(i,:), t, B(i,:) )
  end do
  !
 end subroutine get_block_of_fields
 
 subroutine push_block( block_size, x, p, dt, E, B, xo, po )
  implicit none
  integer, intent(in) :: block_size
  real(8), intent(in), dimension(block_size, 3) :: x, p
  real(8), intent(in) :: dt
  real(8), intent(in), dimension(block_size, 3) :: E, B
  real(8), intent(out), dimension(block_size, 3) :: xo, po
  !
  integer :: i
  !
  xo = x
  po = p
  do i = 1, block_size
    call push( xo(i,:), po(i,:), dt, E(i,:), B(i,:) )
  enddo
  !
 end subroutine push_block
 
 subroutine pushall( particle_count, x, p, t, dt, xo, po )
  implicit none
  integer, intent(in) :: particle_count
  real(8), intent(in), dimension(particle_count, 3) :: x, p
  real(8), intent(in) :: t, dt
  real(8), intent(out), dimension(particle_count, 3) :: xo, po
  !
  integer :: i
  !
  xo = x
  po = p
  do i = 1, particle_count
    call push2( xo(i,:), po(i,:), t, dt)
  enddo
  !
 end subroutine pushall
 
 SUBROUTINE e_field(Y,t,E)
  IMPLICIT NONE
  REAL(8), DIMENSION(1:3), INTENT(IN) :: Y
  REAL(8), INTENT(IN) :: t
  REAL(8), DIMENSION(1:3), INTENT(OUT) :: E
  REAL(8) :: x
  !
  E = 0
  x = Y(1)
  if (x <= 0) then
    E(2) = 6.0*sin(1.0*t - 1.0*x) - 4.0327868852459*sin(1.0*t + 1.0*x) + 2.36065573770492*cos(1.0*t + 1.0*x)
  else
    E(2) = (1.9672131147541*sin(1.0*t - 1.5*x) + 2.36065573770492*cos(1.0*t - 1.5*x))*exp(-3.0*x)
  endif
 END SUBROUTINE e_field
 
 SUBROUTINE b_field(Y,t,B)
  IMPLICIT NONE
  REAL(8), DIMENSION(1:3), INTENT(IN) :: Y
  REAL(8), INTENT(IN) :: t
  REAL(8), DIMENSION(1:3), INTENT(OUT) :: B
  REAL(8) :: x
  !
  B = 0
  x = Y(1)
  if (x <= 0) then
    B(3) = 6.0*sin(1.0*t - 1.0*x) + 4.0327868852459*sin(1.0*t + 1.0*x) - 2.36065573770492*cos(1.0*t + 1.0*x)
  else
    B(3) = (10.0327868852459*sin(1.0*t - 1.5*x) - 2.36065573770492*cos(1.0*t - 1.5*x))*exp(-3.0*x)
  endif
 END SUBROUTINE b_field
 
 SUBROUTINE push(x, p, dt, E, B)
  implicit none
  real(8), intent(inout), dimension(3) :: x, p
  REAL(8), INTENT(IN) :: dt
  real(8), intent(in), dimension(3) :: E, B
  !
  REAL(8), DIMENSION(3) :: u_minus, t, s, u_plus
  real(8) :: rg
  !
  u_minus = p + E*dt*0.5
  rg = sqrt(1.0 + u_minus(1)**2 + u_minus(2)**2 + u_minus(3)**2)
  rg = 1.0/rg
  t = B*dt*0.5*rg
  s = 2*t/(1 + dot_product(t,t))
  call cross_product(u_minus, s, u_plus)
  u_plus = u_plus + u_minus + t*dot_product(u_minus,s)
  p = u_plus + E*dt*0.5
  rg = sqrt(1.0 + p(1)**2 + p(2)**2 + p(3)**2)
  rg = 1.0/rg
  x = x + p*dt*rg
 END SUBROUTINE push
 
 SUBROUTINE push2(x, p, time, dt)
  implicit none
  REAL(8), DIMENSION(3), INTENT(INOUT) :: x, p
  REAL(8), INTENT(IN) :: time, dt
  !
  real(8) :: inv_gamma
  !
  REAL(8), DIMENSION(3) :: u_minus, E, t, B, s, u_plus
  INTEGER :: i
  ! 
  CALL e_field(x,time,E)
  u_minus = p + E*dt*0.5
  CALL b_field(x,time,B)
  t = B*dt*0.5*inv_gamma(u_minus)
  s = 2*t/(1 + dot_product(t,t))
  call cross_product(u_minus, s, u_plus)
  u_plus = u_plus + u_minus + t*dot_product(u_minus,s)
  p = u_plus + E*dt*0.5
  x = x + p*dt*inv_gamma(p)
 END SUBROUTINE push2


 SUBROUTINE cross_product(A,B,C)
  REAL(8), DIMENSION(1:3), INTENT(IN) :: A,B
  REAL(8), DIMENSION(1:3), INTENT(OUT) :: C
  !AxB = C
  C(1) =  A(2)*B(3) - A(3)*B(2)
  C(2) = -A(1)*B(3) + A(3)*B(1)
  C(3) =  A(1)*B(2) - A(2)*B(1)
 END SUBROUTINE cross_product
 
 FUNCTION inv_gamma(p)
  IMPLICIT NONE
  real(8) inv_gamma
  !
  REAL(8), DIMENSION(1:3), INTENT(IN)  :: p
  REAL(8) :: gamma_factor
  !
  gamma_factor = sqrt(1+p(1)**2+p(2)**2+p(3)**2)
  inv_gamma = 1/gamma_factor
 END FUNCTION inv_gamma
