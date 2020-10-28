module type_defs
  integer, parameter:: sp = kind(1.0),&
    dp = selected_real_kind(2*precision(1.0_sp)),&
    qp = selected_real_kind(2*precision(1.0_dp))
end module type_defs

module problem_setup
  use type_defs
  implicit none
  integer,  parameter :: Nx = 20
end module problem_setup

module arrs
  use type_defs
  implicit none
  real(dp), allocatable, dimension(:) :: u,b,x
end module arrs

module afuns
  
contains
  subroutine apply_1D_laplacian(au,u,n,hx)
    use type_defs
    implicit none
    integer, intent(in) :: n
    real(dp), intent(out) :: au(n)
    real(dp), intent(in)  ::  u(n),hx
    real(dp) :: hxi2
    integer :: i
    hxi2 = 1.0_dp/hx**2
    Au(1) = hxi2*(u(2) - 2.0_dp*u(1)         )
    Au(n) = hxi2*(     - 2.0_dp*u(n) + u(n-1))
    do i = 2,n-1
     Au(i)= hxi2*(u(i+1) - 2.0_dp*u(i) + u(i-1))
    end do
    
  end subroutine apply_1D_laplacian

end module afuns

module iterative_solvers
  use type_defs
  implicit none
  real(dp), parameter :: TOL = 1.0e-3_dp
  
contains
  
  subroutine richardson(x,b,nx,hx,l)
    use type_defs
    use afuns
    implicit none
    integer, intent(in) :: nx
    real(dp), intent(inout)  :: hx
    real(dp), intent(inout)  :: x(nx)
    real(dp), intent(inout)  :: b(nx)
    integer, intent(out) :: l
    real(dp) :: ax(nx),residual(nx)
    real(dp) :: res_norm

    ! loop
    res_norm = 2*TOL
    l = 0
    x = 0.0d0
    do while (res_norm .gt. TOL)
       ! One iteration
       call apply_1D_laplacian(ax,x,nx+1,hx)
       residual = (ax-b)
       x  = x + 0.1*hx**2*residual
       res_norm = sqrt(sum(residual**2))
       l = l+1
    end do
  end subroutine richardson
  
end module iterative_solvers

program ins
  use type_defs
  use problem_setup
  use arrs
  use afuns
  use iterative_solvers
  implicit none
  ! This program solves u_xx = b 
  ! on the domain [x] \in [0,1] with zero boundary conditions 
  ! hx = 1/Nx
  real(dp) :: hx
  logical, parameter :: iterative_method = .true.
  integer :: i,n_iter,N_sys,info  
  real(dp), allocatable, dimension(:,:) :: A
  real(dp), allocatable, dimension(:) :: action_of_A,u_inner
  integer, allocatable, dimension(:) ::  ipiv
  ! Set up the grid
  hx = 1.0_dp/real(Nx,dp)
  allocate(x(0:nx))
  do i = 0,nx
   x(i) = real(i,dp)*hx
  end do
  allocate(u(0:nx),b(1:nx-1))
  allocate(A(nx-1,nx-1),&
       action_of_A(nx-1),&
       u_inner(nx-1),&
       ipiv(nx-1))
  
  n_sys = nx-1
  do i = 1,nx-1
     u_inner = 0.0
     u_inner(i) = 1.0d0
     call apply_1D_laplacian(action_of_A,u_inner,nx-1,hx)
     A(:,i) = action_of_A
  end do
  b = -2.0_dp
  ! u^{n+1} = u^{n} + (A*u^{n}-b), n = 0, 1, 2... 
  
!!$  u_inner = 0.0
!!$  write(*,*) u_inner
!!$  do i = 1,10000
!!$     ! One iteration
!!$     call apply_1D_laplacian(action_of_A,u_inner,nx-1,hx)
!!$     ! Fixed point iteration
!!$     u_inner  = u_inner + 0.1*hx**2*(action_of_A-b)
!!$     write(*,*) u_inner-x(1:nx-1)*(1.0_dp-x(1:nx-1))
!!$  end do
!!$  stop 123

  if (iterative_method) then
     u = 0.0_dp
     call richardson(u_inner,b,nx-1,hx,n_iter)
     u(1:nx-1) = u_inner
  else
     CALL DGETRF(N_sys,N_sys,A,N_sys,ipiv,INFO)
     CALL DGETRS('N',N_sys,1,A,N_sys,IPIV,b,N_sys,INFO)
     u = 0.d0
     u(1:nx-1) = b 
     n_iter = -1
  end if
 write(*,*) maxval(abs(u - x*(1.0_dp-x))), n_iter
end program ins
