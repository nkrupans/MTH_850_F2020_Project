module type_defs
  integer, parameter:: sp = kind(1.0),&
    dp = selected_real_kind(2*precision(1.0_sp)),&
    qp = selected_real_kind(2*precision(1.0_dp))
end module type_defs

module problem_setup
  use type_defs
  implicit none
  real(dp), parameter :: Lx = 1.0_dp
  real(dp), parameter :: Ly = 1.0_dp
  integer,  parameter :: Nx = 50          ! Number of gridpoints in x
  integer,  parameter :: Ny = 50          ! Number of gridpoints in y
  integer,  parameter :: Nsteps = 5000    ! Number of timesteps  
  logical, parameter :: do_plot = .true.  ! Plot?
  integer,  parameter :: Nplot = 100      ! If so plot every Nplot steps  
  real(dp), parameter :: k = 0.01_dp      ! Timestep 
  real(dp), parameter :: alpha = 0.1_dp/k 
  real(dp), parameter :: nu = 0.01_dp     ! Viscosity 
  real(dp), parameter :: pi = acos(-1.d0)
  
end module problem_setup

module arrs
  use type_defs
  implicit none
  real(dp), allocatable, dimension(:) :: x,y
  real(dp), allocatable, dimension(:,:) :: u,v,p
  real(dp), allocatable, dimension(:,:) :: uold,vold,pold
  real(dp), allocatable, dimension(:,:) :: leu,lev,liu,liv
  real(dp), allocatable, dimension(:,:) :: leuold,levold,liuold,livold
  real(dp), allocatable, dimension(:,:) :: gux,guy,gvx,gvy,pbx,pby
end module arrs

module matrices
  use type_defs
  implicit none
  real(dp), allocatable, dimension(:,:) :: LapUV,LapP,LapPBig
  integer,  allocatable, dimension(:) :: ipiv_pbig,ipiv_uv
  real(dp), allocatable, dimension(:) :: pvec,lpvec,uvec,vvec
  real(dp), allocatable, dimension(:) :: pbvecbig,pvecbig
end module matrices

module afuns
contains

  subroutine apply_pressure_laplacian(lp,p,nx,ny,hx,hy)
    use type_defs
    implicit none
    integer, intent(in) :: nx,ny
    real(dp), intent(out) :: lp((nx+1)*(ny+1))
    real(dp), intent(in)  ::  p((nx+1)*(ny+1)),hx,hy
    real(dp) :: hxi2,hyi2
    integer :: mx,my,i,j,k,kp,km,kt,kb
    ! Suboptimal quick and dirty implementation ...
    ! It would be better to avoid all the conditional statements
    !
    hxi2 = 1.0_dp/hx**2
    hyi2 = 1.0_dp/hy**2

    mx = (nx+1)
    my = (ny+1)
    do j=1,my
     do i=1,mx
      k = i+(j-1)*mx
      kp = k+1
      km = k-1
      kb = k-mx
      kt = k+mx
      lp(k) = -2.0_dp*(hxi2+hyi2)*p(k)
      if (i.lt.mx) then
       if (i.eq.1) then
        lp(k) = lp(k) + 2.0_dp*hxi2*p(kp)
       else
        lp(k) = lp(k) + hxi2*p(kp)
       end if
      end if
      if (i.gt.1) then
       if (i.eq.mx) then
        lp(k) = lp(k) + 2.0_dp*hxi2*p(km)
       else
        lp(k) = lp(k) + hxi2*p(km)
       end if
      end if
      if (j.lt.my) then
       if (j.eq.1) then
        lp(k) = lp(k) + 2.0_dp*hyi2*p(kt)
       else
        lp(k) = lp(k) + hyi2*p(kt)
       end if
      end if
      if (j.gt.1) then
       if (j.eq.my) then
        lp(k) = lp(k) + 2.0_dp*hyi2*p(kb)
       else
        lp(k) = lp(k) + hyi2*p(kb)
       end if
      end if
     end do
    end do
  end subroutine apply_pressure_laplacian

  subroutine apply_velocity_laplacian(lu,u,nx,ny,hx,hy)
    use type_defs
    implicit none
    integer, intent(in) :: nx,ny
    real(dp), intent(out) :: lu((nx-1)*(ny-1))
    real(dp), intent(in)  ::  u((nx-1)*(ny-1)),hx,hy
    real(dp) :: hxi2,hyi2
    integer :: mx,my,i,j,k,kp,km,kt,kb

    ! Suboptimal quick and dirty implementation ...
    ! It would be better to avoid all the conditional statements

    hxi2 = 1.0_dp/hx**2
    hyi2 = 1.0_dp/hy**2
    mx = (nx-1)
    my = (ny-1)

    do j=1,my
     do i=1,mx
      k = i+(j-1)*mx
      kp = k+1
      km = k-1
      kb = k-mx
      kt = k+mx
      lu(k) = -2.0_dp*(hxi2+hyi2)*u(k)
      if (i.lt.mx) lu(k) = lu(k) + hxi2*u(kp)
      if (i.gt.1)  lu(k) = lu(k) + hxi2*u(km)
      if (j.lt.my) lu(k) = lu(k) + hyi2*u(kt)
      if (j.gt.1)  lu(k) = lu(k) + hyi2*u(kb)
     end do
    end do
  end subroutine apply_velocity_laplacian

end module afuns


program ins
  use type_defs
  use problem_setup
  use arrs
  use matrices
  use afuns
  implicit none
  ! This program solves the incompressible Navier-Stokes equations on
  ! the domain [x,y] \in [0,Lx] \times [0,Ly] using a grid spacing
  ! hx=Lx/Nx, hy = Ly/Ny.
  real(dp) :: hx,hy,time1,time2
  integer :: i,j,sys_size_p,sys_size_pbig,info,nt,sys_size_uv
  character(100) :: str

  ! Set up the grid, we include the ghost points.
  hx = Lx/real(Nx,dp)
  hy = Ly/real(Ny,dp)
  allocate(x(-1:nx+1),y(-1:ny+1))
  do i = -1,nx+1
   x(i) = real(i,dp)*hx
  end do
  do j = -1,ny+1
   y(j) = real(j,dp)*hy
  end do

  ! We store u,v,p as two dimensional arrays.
  allocate(u(-1:nx+1,-1:ny+1),v(-1:nx+1,-1:ny+1),p(-1:nx+1,-1:ny+1))
  allocate(uold(-1:nx+1,-1:ny+1),vold(-1:nx+1,-1:ny+1),pold(-1:nx+1,-1:ny+1))
  ! Arrays for forcing and other stuff
  allocate(leu(0:nx,0:ny),lev(0:nx,0:ny),liu(0:nx,0:ny),liv(0:nx,0:ny))
  allocate(leuold(0:nx,0:ny),levold(0:nx,0:ny),liuold(0:nx,0:ny),livold(0:nx,0:ny))

  write(*,*) 'Setting up Laplacians'
  call cpu_time(time1)
  ! setup pressure laplacian
  allocate(LapP((Nx+1)*(Ny+1),(Nx+1)*(Ny+1)))
  allocate(pvec((Nx+1)*(Ny+1)),lpvec((Nx+1)*(Ny+1)))
  allocate(LapUV((Nx-1)*(Ny-1),(Nx-1)*(Ny-1)))
  allocate(uvec((Nx-1)*(Ny-1)),vvec((Nx-1)*(Ny-1)))
  allocate(pbvecbig((Nx+1)*(Ny+1)+1),pvecbig((Nx+1)*(Ny+1)+1))
  !
  pvec = 0.0_dp
  do i = 1,(nx+1)*(ny+1)
   pvec(i) = 1.0_dp
   call apply_pressure_laplacian(lpvec,pvec,nx,ny,hx,hy)
   pvec(i) = 0.0_dp
   LapP(:,i) = lpvec
  end do
  ! call printdble2d(LapP,1,(nx+1)*(ny+1),1,(nx+1)*(ny+1),'LapP.txt')
  ! setup uv laplacian
  uvec = 0.0_dp
  do i = 1,(nx-1)*(ny-1)
   uvec(i) = 1.0_dp
   call apply_velocity_laplacian(vvec,uvec,nx,ny,hx,hy)
   uvec(i) = 0.0_dp
   LapUV(:,i) = -0.5_dp*nu*vvec
   LapUV(i,i) = LapUV(i,i) + 1.0_dp/k
  end do
  ! call printdble2d(LapUV,1,(nx-1)*(ny-1),1,(nx-1)*(ny-1),'LapUV.txt')
  allocate(LapPBig((Nx+1)*(Ny+1)+1,(Nx+1)*(Ny+1)+1))
  LapPBig = 1.0_dp
  LapPBig((Nx+1)*(Ny+1)+1,(Nx+1)*(Ny+1)+1) = 0.0_dp
  LapPBig(1:(Nx+1)*(Ny+1),1:(Nx+1)*(Ny+1)) = LapP
  call cpu_time(time2)
  write(*,*) '... done setting up Laplacians it took ',time2-time1 ,' seconds'

  ! Set up dense linear algebra stuff.
  sys_size_p = (nx+1)*(ny+1)
  sys_size_pbig = sys_size_p + 1
  allocate(ipiv_pbig(sys_size_pbig))
  sys_size_uv = (nx-1)*(ny-1)
  allocate(ipiv_uv(sys_size_uv))
  write(*,*) 'Factoring matrices'
  call cpu_time(time1)
  ! Factor
  CALL DGETRF(sys_size_pbig,sys_size_pbig,LapPbig,sys_size_pbig,ipiv_pbig,INFO)
  CALL DGETRF(sys_size_uv,sys_size_uv,LapUV,sys_size_uv,ipiv_uv,INFO)
  call cpu_time(time2)
  write(*,*) '... factorization took ',time2-time1 ,' seconds'

  ! These are the arrays for the boundary forcing u=gu, v=gv on the
  ! boundary. We order them
  ! gux(1,:) = "left"
  ! gux(2,:) = "right"
  ! guy(1,:) = "bottom"
  ! guy(2,:) = "top"
  allocate(gux(2,0:ny),guy(2,0:nx),gvx(2,0:ny),gvy(2,0:nx))
  ! These hold boundary forcings
  allocate(pbx(2,0:ny),pby(2,0:nx))
  ! Boundary conditions must be constant in time in current implementation.
  gux = 0.0_dp
  guy = 0.0_dp
  gvx = 0.0_dp
  gvy = 0.0_dp

  ! Various initial data
  u = 0.d0
  v = 0.d0

  call random_number(u(0:nx,0:ny))
  call random_number(v(0:nx,0:ny))
  u(0:nx,0:ny) = u(0:nx,0:ny) - 0.5_dp
  v(0:nx,0:ny) = v(0:nx,0:ny) - 0.5_dp

  do j = 0,ny
   do i = 0,nx
    ! u(i,j) = sin(pi*x(i))*sin(pi*y(j))
    ! v(i,j) = sin(pi*x(i))*sin(pi*y(j))
    call taylor(u(i,j),v(i,j),X(i),Y(j),0.5d0,0.5d0,0.1d0,1.0_dp)
   end do
  end do

  ! Initial data for a lid-driven cavity flow
  guy(2,:) = 1.0_dp
  u = 0.d0
  v = 0.d0

  ! start up the computation with a single Euler step (backwards in time).
  ! We need to find the current pressure to compute the advection term.
  uold = 0.0_dp
  vold = 0.0_dp
  pold = 0.0_dp
  call updateBCforU(u,v,gux,guy,gvx,gvy,nx,ny)
  call computeAndUpdateGPforU(u,v,hx,hy,nx,ny)
  call computeGPforP(pbx,pby,u,v,gux,guy,gvx,gvy,hx,hy,nu,nx,ny)
  call setupRhsideP(pbvecbig,nx,ny,hx,hy,u,v,alpha,pbx,pby)

  ! !!! YOUR CODE REPLACES THIS !!!!
  ! Solve for p
  CALL DGETRS('N',sys_size_pbig,1,LapPbig,sys_size_pbig,IPIV_pbig,&
    pbvecbig,sys_size_pbig,INFO)
  ! !!! END YOUR CODE REPLACES THIS !!!!

  ! Swap long vector into 2D array
  p = 0.0d0
  do j = 0,ny
   do i = 0,nx
    p(i,j) = pbvecbig(1+i+j*(nx+1))
   end do
  end do


  call computeLE(Leu,Lev,u,v,p,hx,hy,nx,ny)
  call computeLI(Liu,Liv,u,v,nu,hx,hy,nx,ny)
  uold(0:nx,0:ny) = u(0:nx,0:ny) - k*(Leu + Liu)
  vold(0:nx,0:ny) = v(0:nx,0:ny) - k*(Lev + Liv)
  call updateBCforU(uold,vold,gux,guy,gvx,gvy,nx,ny)
  call computeAndUpdateGPforU(uold,vold,hx,hy,nx,ny)
  call computeGPforP(pbx,pby,uold,vold,gux,guy,gvx,gvy,hx,hy,nu,nx,ny)
  call setupRhsideP(pbvecbig,nx,ny,hx,hy,uold,vold,alpha,pbx,pby)

  ! !!! YOUR CODE REPLACES THIS !!!!
  ! Solve for p
  CALL DGETRS('N',sys_size_pbig,1,LapPbig,sys_size_pbig,IPIV_pbig,pbvecbig,&
    sys_size_pbig,INFO)
  ! !!! END YOUR CODE REPLACES THIS !!!!

  ! Swap long vector into 2D array
  do j = 0,ny
   do i = 0,nx
    pold(i,j) = pbvecbig(1+i+j*(nx+1))
   end do
  end do
  call computeLE(Leuold,Levold,uold,vold,pold,hx,hy,nx,ny)
  call computeLI(Liuold,Livold,uold,vold,nu,hx,hy,nx,ny)

  write(*,*) 'Starting Time Loop....'
  call cpu_time(time1)
  do nt = 1,Nsteps
   ! semi-implicit method
   ! Get boundary conditions at time n+1 ;
   ! Boundary conditions are constant in time in current implementation so no update
   call setupRhsideUV(Uvec,Vvec,nx,ny,hx,hy,gux,guy,gvx,gvy,&
     Leu,Lev,Leuold,Levold,Liu,Liv,u,v,k,nu)
   ! Swap
   uold = u
   vold = v
   ! solve for new u and v

   ! !!! YOUR CODE REPLACES THIS !!!!
   CALL DGETRS('N',sys_size_uv,1,Lapuv,sys_size_uv,IPIV_uv,uvec,sys_size_uv,INFO)
   CALL DGETRS('N',sys_size_uv,1,Lapuv,sys_size_uv,IPIV_uv,vvec,sys_size_uv,INFO)
   ! !!! END YOUR CODE REPLACES THIS !!!!

   ! Swap long vector into 2D array
   do j = 1,ny-1
    do i = 1,nx-1
     u(i,j) = uvec(i+(j-1)*(nx-1))
     v(i,j) = vvec(i+(j-1)*(nx-1))
    end do
   end do
   !
   call updateBCforU(u,v,gux,guy,gvx,gvy,nx,ny)
   call computeAndUpdateGPforU(u,v,hx,hy,nx,ny)
   call computeGPforP(pbx,pby,u,v,gux,guy,gvx,gvy,hx,hy,nu,nx,ny)
   call setupRhsideP(pbvecbig,nx,ny,hx,hy,u,v,alpha,pbx,pby)

   ! !!! YOUR CODE REPLACES THIS !!!!
   ! Solve for p
   CALL DGETRS('N',sys_size_pbig,1,LapPbig,sys_size_pbig,IPIV_pbig,&
     pbvecbig,sys_size_pbig,INFO)
   ! !!! END YOUR CODE REPLACES THIS !!!!

   ! Swap long vector into 2D array
   do j = 0,ny
    do i = 0,nx
     p(i,j) = pbvecbig(1+i+j*(nx+1))
    end do
   end do
   ! Swap
   Liuold=Liu
   Livold=Liv
   Leuold=Leu
   Levold=Lev
   ! Compute new "operators"
   call computeLE(Leu,Lev,u,v,p,hx,hy,nx,ny)
   call computeLI(Liu,Liv,u,v,nu,hx,hy,nx,ny)
   ! Plot stuff
   if (do_plot) then
    if (mod(nt,nplot) .eq. 0 ) then
     WRITE(str,'("u",I8.8,".txt")') nt
     call printdble2d(u(0:nx,0:ny),0,nx,0,ny,trim(str))
     WRITE(str,'("v",I8.8,".txt")') nt
     call printdble2d(v(0:nx,0:ny),0,nx,0,ny,trim(str))
    end if
   end if

   if (maxval(abs(u(0:nx,0:ny))) .gt. 1.0d6) then
    write(*,*) "Solution blew up?"
    stop 123
   end if
   if (mod(nt,100).eq.0) then
    write(*,*) "Timestep ", nt, " out of ", nsteps
   end if

  end do
  call cpu_time(time2)
  write(*,*) 'Time loop took ',time2-time1 ,' seconds'

end program ins


subroutine computeLE(Leu,Lev,u,v,p,hx,hy,nx,ny)
  ! This routine computes the explicit part of the right hand side
  use type_defs
  implicit none
  integer, intent(in)   :: nx,ny
  real(dp), intent(in)  :: hx,hy
  real(dp), intent(out) :: leu(0:nx,0:ny),lev(0:nx,0:ny)
  real(dp), intent(in)  :: u(-1:nx+1,-1:ny+1),v(-1:nx+1,-1:ny+1),p(-1:nx+1,-1:ny+1)
  real(dp) :: hx2i,hy2i
  integer :: i,j

  hx2i = 0.5_dp/hx
  hy2i = 0.5_dp/hy
  do j = 0,ny
   do i = 0,nx
    Leu(i,j) = -(&
      + u(i,j)*(u(i+1,j)-u(i-1,j))*hx2i &
      + v(i,j)*(u(i,j+1)-u(i,j-1))*hy2i &
      + (p(i+1,j)-p(i-1,j))*hx2i)
    Lev(i,j) = -(&
      + u(i,j)*(v(i+1,j)-v(i-1,j))*hx2i &
      + v(i,j)*(v(i,j+1)-v(i,j-1))*hy2i &
      + (p(i,j+1)-p(i,j-1))*hy2i)
   end do
  end do

end subroutine computeLE

subroutine computeLI(Liu,Liv,u,v,nu,hx,hy,nx,ny)
  ! This routine computes the implicit part of the right hand side
  use type_defs
  implicit none
  integer, intent(in)   :: nx,ny
  real(dp), intent(in)  :: hx,hy,nu
  real(dp), intent(out) :: liu(0:nx,0:ny),liv(0:nx,0:ny)
  real(dp), intent(in)  :: u(-1:nx+1,-1:ny+1),v(-1:nx+1,-1:ny+1)
  real(dp) :: hxi2,hyi2
  integer :: i,j

  hxi2 = nu/hx**2
  hyi2 = nu/hy**2
  do j = 0,ny
   do i = 0,nx
    Liu(i,j) = (u(i+1,j)-2.0_dp*u(i,j)+u(i-1,j))*hxi2 &
      +(u(i,j+1)-2.0_dp*u(i,j)+u(i,j-1))*hyi2
    Liv(i,j) = (v(i+1,j)-2.0_dp*v(i,j)+v(i-1,j))*hxi2 &
      +(v(i,j+1)-2.0_dp*v(i,j)+v(i,j-1))*hyi2
   end do
  end do
end subroutine computeLI

subroutine updateBCforU(u,v,gux,guy,gvx,gvy,nx,ny)
  ! We update u and v to the left and right.
  use type_defs
  implicit none
  integer, intent(in)     :: nx,ny
  real(dp), intent(in)    :: gux(2,0:ny),guy(2,0:nx),gvx(2,0:ny),gvy(2,0:nx)
  real(dp), intent(inout) :: u(-1:nx+1,-1:ny+1),v(-1:nx+1,-1:ny+1)

  ! Left and right.
  u(0,0:ny)  = gux(1,0:ny)
  u(nx,0:ny) = gux(2,0:ny)
  v(0,0:ny)  = gvx(1,0:ny)
  v(nx,0:ny) = gvx(2,0:ny)
  ! Then we do the top and bottom.
  u(0:nx,0)  = guy(1,0:nx)
  u(0:nx,ny) = guy(2,0:nx)
  v(0:nx,0)  = gvy(1,0:nx)
  v(0:nx,ny) = gvy(2,0:nx)

end subroutine updateBCforU

subroutine computeAndUpdateGPforU(u,v,hx,hy,nx,ny)
  ! We update u and v to the left and right.
  use type_defs
  implicit none
  integer, intent(in)     :: nx,ny
  real(dp), intent(in)    :: hx,hy
  real(dp), intent(inout) :: u(-1:nx+1,-1:ny+1),v(-1:nx+1,-1:ny+1)
  real(dp) :: vbx(2,0:ny), uby(2,0:nx), ubx(2,0:ny), vby(2,0:nx)
  ! This routine assumes the interior point and the physical boundary
  ! conditions have been updated.
  ! We extrapolate v to the left and right
  vbx(1,0:ny) = 3.0_dp*v(0,0:ny)  - 3.0_dp*v(1,0:ny)    + v(2,0:ny)
  vbx(2,0:ny) = 3.0_dp*v(nx,0:ny) - 3.0_dp*v(nx-1,0:ny) + v(nx-2,0:ny)
  ! We extrapolate u to the bottom and top
  uby(1,0:nx) = 3.0_dp*u(0:nx,0)  - 3.0_dp*u(0:nx,1)    + u(0:nx,2)
  uby(2,0:nx) = 3.0_dp*u(0:nx,ny) - 3.0_dp*u(0:nx,ny-1) + u(0:nx,ny-2)
  ! update ghost point values
  v(0,0:ny)  = vbx(1,0:ny)
  v(nx,0:ny) = vbx(2,0:ny)
  u(0:nx,0)  = uby(1,0:nx)
  u(0:nx,ny) = uby(2,0:nx)
  ! we also need to extrapolate data for the stencil
  u(-1,0)    = 3.0_dp*u(0,0)   - 3.0_dp*u(1,0)     + u(2,0)
  u(nx+1,0)  = 3.0_dp*u(nx,0)  - 3.0_dp*u(nx-1,0)  + u(nx-2,0)
  u(-1,ny)   = 3.0_dp*u(0,ny)  - 3.0_dp*u(1,ny)    + u(2,ny)
  u(nx+1,ny) = 3.0_dp*u(nx,ny) - 3.0_dp*u(nx-1,ny) + u(nx-2,ny)
  !
  v(0,-1)    = 3.0_dp*v(0,0)   - 3.0_dp*v(0,1)     + v(0,2)
  v(0,ny+1)  = 3.0_dp*v(0,ny)  - 3.0_dp*v(0,ny-1)  + v(0,ny-2)
  v(nx,-1)   = 3.0_dp*v(nx,0)  - 3.0_dp*v(nx,1)    + v(nx,2)
  v(nx,ny+1) = 3.0_dp*v(nx,ny) - 3.0_dp*v(nx,ny-1) + v(nx,ny-2)
  ! To the left and right we get u from the zero divergence condition
  ! to the left
  ubx(1,0:ny) = u(1,0:ny)    + (hx/hy)*(v(0,1:ny+1)-v(0,-1:ny-1))
  ! to the right
  ubx(2,0:ny) = u(nx-1,0:ny) - (hx/hy)*(v(nx,1:ny+1)-v(nx,-1:ny-1))
  ! At the bottom and top we get v from the zero divergence condition
  ! at the bottom
  vby(1,0:nx) = v(0:nx,1)    + (hy/hx)*(u(1:nx+1,0)-u(-1:nx-1,0))
  vby(2,0:nx) = v(0:nx,ny-1) - (hy/hx)*(u(1:nx+1,ny)-u(-1:nx-1,ny))
  ! update ghost point values
  u(-1,0:ny)   = ubx(1,0:ny)
  u(nx+1,0:ny) = ubx(2,0:ny)
  v(0:nx,-1)   = vby(1,0:nx)
  v(0:nx,ny+1) = vby(2,0:nx)
  ! Finally we extrapolate to the corners
  u(-1,-1)     = 3.0_dp*u(0,-1)    - 3.0_dp*u(1,-1)      + u(2,-1)
  u(nx+1,-1)   = 3.0_dp*u(nx,-1)   - 3.0_dp*u(nx-1,-1)   + u(nx-2,-1)
  u(-1,ny+1)   = 3.0_dp*u(0,ny+1)  - 3.0_dp*u(1,ny+1)    + u(2,ny+1)
  u(nx+1,ny+1) = 3.0_dp*u(nx,ny+1) - 3.0_dp*u(nx-1,ny+1) + u(nx-2,ny+1)
  v(-1,-1)     = 3.0_dp*v(-1,0)    - 3.0_dp*v(-1,1)      + v(-1,2)
  v(-1,ny+1)   = 3.0_dp*v(-1,ny)   - 3.0_dp*v(-1,ny-1)   + v(-1,ny-2)
  v(nx+1,-1)   = 3.0_dp*v(nx+1,0)  - 3.0_dp*v(nx+1,1)    + v(nx+1,2)
  v(nx+1,ny+1) = 3.0_dp*v(nx+1,ny) - 3.0_dp*v(nx+1,ny-1) + v(nx+1,ny-2)

end subroutine computeAndUpdateGPforU

subroutine computeGPforP(pbx,pby,u,v,gux,guy,gvx,gvy,hx,hy,nu,nx,ny)
  use type_defs
  implicit none
  integer, intent(in)     :: nx,ny
  real(dp), intent(in)    :: hx,hy,nu
  real(dp), intent(inout) :: u(-1:nx+1,-1:ny+1),v(-1:nx+1,-1:ny+1)
  real(dp), intent(in)    :: gux(2,0:ny),guy(2,0:nx),gvx(2,0:ny),gvy(2,0:nx)
  real(dp), intent(out)   :: pbx(2,0:ny), pby(2,0:nx)


  ! We compute "ghost values" for p
  ! pbx will be used to the right hand side in the Poisson eq.
  ! This approximation uses the curl-curl condition
  ! left
  pbx(1,0:ny) = -2.0_dp*hx*(&
    -gux(1,0:ny)*(u(1,0:ny)-u(-1,0:ny))/(2.0_dp**hx) &
    -gvx(1,0:ny)*(u(0,1:ny+1)-u(0,-1:ny-1))/(2.0_dp*hy) &
    +nu*(&
    -(v(1,1:ny+1)-v(-1,1:ny+1)-v(1,-1:ny-1)+v(-1,-1:ny-1))/(4.0_dp*hx*hy) &
    +(u(0,1:ny+1)-2.0_dp*u(0,0:ny)+u(0,-1:ny-1))/(hy**2)))
  ! right
  pbx(2,0:ny) = +2.0_dp*hx*(&
    -gux(2,0:ny)*(u(nx+1,0:ny)-u(nx-1,0:ny))/(2.0_dp*hx) &
    -gvx(2,0:ny)*(u(nx,1:ny+1)-u(nx,-1:ny-1))/(2.0_dp*hy) &
    +nu*(&
    -(v(nx+1,1:ny+1)-v(nx-1,1:ny+1)-v(nx+1,-1:ny-1)+v(nx-1,-1:ny-1))/(4.0_dp*hx*hy) &
    +(u(nx,1:ny+1)-2.0_dp*u(nx,0:ny)+u(nx,-1:ny-1))/(hy**2)))

  ! bottom
  pby(1,0:nx) = -2.0_dp*hy*(&
    -gvy(1,0:nx)*(v(0:nx,1)   - v(0:nx,-1))/(2.0_dp*hy) &
    -guy(1,0:nx)*(v(1:nx+1,0) - v(-1:nx-1,0))/(2.0_dp*hx) &
    +nu*(&
    +(v(1:nx+1,1)-v(-1:nx-1,1)-v(1:nx+1,-1)+v(-1:nx-1,-1))/(4.0_dp*hy*hx) &
    +(v(1:nx+1,0)-2.0_dp*v(0:nx,0)+v(-1:nx-1,0))/(hx**2)))
  ! top
  pby(2,0:nx) = +2.0_dp*hy*(&
    -gvy(2,0:nx)*(v(0:nx,ny+1)-v(0:nx,ny-1))/(2.0_dp*hy) &
    -guy(2,0:nx)*(v(1:nx+1,ny)-v(-1:nx-1,ny))/(2.0_dp*hx) &
    +nu*(&
    + (v(1:nx+1,ny+1)-v(1:nx+1,ny-1)-v(-1:nx-1,ny+1)+v(-1:nx-1,ny-1))/(4.0_dp*hy*hx) &
    + (v(1:nx+1,ny)-2.0_dp*v(0:nx,ny)+v(-1:nx-1,ny))/(hx**2)))

end subroutine computeGPforP

subroutine setupRhsideP(bP,nx,ny,hx,hy,u,v,alpha,pbx,pby)
  use type_defs
  implicit none
  integer, intent(in)     :: nx,ny
  real(dp), intent(in)    :: hx,hy,alpha
  real(dp), intent(in)    :: u(-1:nx+1,-1:ny+1),v(-1:nx+1,-1:ny+1)
  real(dp), intent(in)    :: pbx(2,0:ny), pby(2,0:nx)
  real(dp), intent(out)    :: bP((Nx+1)*(Ny+1))
  real(dp) :: D0xu(0:nx,0:ny),D0xv(0:nx,0:ny),D0yu(0:nx,0:ny),D0yv(0:nx,0:ny),F(0:nx,0:ny)
  real(dp) :: hxi2,hyi2
  integer :: i,j,k

  hxi2 = 1.0_dp/hx**2
  hyi2 = 1.0_dp/hy**2

  D0xu = (u(1:nx+1,0:ny)-u(-1:nx-1,0:ny))*(1.0_dp/(2.0_dp*hx))
  D0xv = (v(1:nx+1,0:ny)-v(-1:nx-1,0:ny))*(1.0_dp/(2.0_dp*hx))
  D0yu = (u(0:nx,1:ny+1)-u(0:nx,-1:ny-1))*(1.0_dp/(2.0_dp*hy))
  D0yv = (v(0:nx,1:ny+1)-v(0:nx,-1:ny-1))*(1.0_dp/(2.0_dp*hy))

  F = alpha*(D0xu+D0yv)-(D0xu**2+D0yv**2+2*D0xv*D0yu)

  ! inner points
  do j = 1,ny-1
   do i = 1,nx-1
    k = 1+i + j*(nx+1)
    bP(k) = F(i,j)
   end do
  end do

  ! Bottom and top sides
  do j = 0,ny,ny
   do i=1,nx-1
    k = 1+i + j*(nx+1)
    bP(k)=F(i,j)
    if (j.eq.0)  bP(k) = bP(k) - pby(1,i)*hyi2
    if (j.eq.ny) bP(k) = bP(k) - pby(2,i)*hyi2
   end do
  end do
  ! Left and right sides
  do j=1,ny-1
   do i=0,nx,nx
    k = 1 + i + j*(nx+1)
    bP(k)=F(i,j)
    if (i.eq.0)  bP(k) = bP(k) - pbx(1,j)*hxi2
    if (i.eq.nx) bP(k) = bP(k) - pbx(2,j)*hxi2
   end do
  end do

  ! corners
  do j = 0,ny,ny
   do i = 0,nx,nx
    k =  1 + i + j*(nx+1)
    bP(k) = F(i,j)
    if (i.eq.0)  bP(k)=bP(k)-pbx(1,j)*hxi2
    if (i.eq.nx) bP(k)=bP(k)-pbx(2,j)*hxi2
    if (j.eq.0)  bP(k)=bP(k)-pby(1,i)*hyi2
    if (j.eq.ny) bP(k)=bP(k)-pby(2,i)*hyi2
   end do
  end do

end subroutine setupRhsideP

subroutine setupRhsideUV(bU,bV,nx,ny,hx,hy,gux,guy,gvx,gvy,Leu,Lev,Leuold,Levold,&
  Liu,Liv,u,v,k,nu)
  use type_defs
  implicit none
  integer, intent(in)   :: nx,ny
  real(dp), intent(in)  :: hx,hy,nu,k
  real(dp), intent(in)  :: leu(0:nx,0:ny),lev(0:nx,0:ny),leuold(0:nx,0:ny),levold(0:nx,0:ny)
  real(dp), intent(in)  :: liu(0:nx,0:ny),liv(0:nx,0:ny)
  real(dp), intent(in)  :: u(-1:nx+1,-1:ny+1),v(-1:nx+1,-1:ny+1)
  real(dp), intent(in)  :: gux(2,0:ny),guy(2,0:nx),gvx(2,0:ny),gvy(2,0:nx)
  real(dp), intent(out) :: bu((Nx-1)*(Ny-1)),bv((Nx-1)*(Ny-1))
  real(dp)  :: ffu(1:nx-1,1:ny-1),ffv(1:nx-1,1:ny-1)
  real(dp) :: hxi2,hyi2
  integer :: i,j,l

  bu = 0.0_dp
  bv = 0.0_dp
  hxi2 = 1.0_dp/hx**2
  hyi2 = 1.0_dp/hy**2
  ! loop over side
  do j = 2,ny-1
   do i = 1,nx-1,nx-1
    l = i + (j-1)*(nx-1)
    if (i.eq.1) then
     bU(l) = bU(l) - gux(1,j)*hxi2
     bV(l) = bV(l) - gvx(1,j)*hxi2
    end if
    if (i.eq.nx-1) then
     bU(l) = bU(l) - gux(2,j)*hxi2
     bV(l) = bV(l) - gvx(2,j)*hxi2
    end if
   end do
  end do

  ! loop over side
  do j = 1,ny-1,ny-1
   do i = 1,nx-1
    l = i + (j-1)*(nx-1)
    if (j.eq.1) then
     bU(l) = bU(l) - guy(1,i)*hyi2
     bV(l) = bV(l) - gvy(1,i)*hyi2
    end if
    if (j.eq.ny-1) then
     bU(l) = bU(l) - guy(2,i)*hyi2
     bV(l) = bV(l) - gvy(2,i)*hyi2
    end if
   end do
  end do
  ! loop over corners
  do j = 1,ny-1,ny-1
   do i = 1,nx-1,nx-1
    l = i + (j-1)*(nx-1)
    if (i.eq.1) then
     bU(l) = bU(l) - gux(1,j)*hxi2
     bV(l) = bV(l) - gvx(1,j)*hxi2
    end if
    if (i.eq.nx-1) then
     bU(l) = bU(l) - gux(2,j)*hxi2
     bV(l) = bV(l) - gvx(2,j)*hxi2
    end if
    if (j.eq.1) then
     bU(l) = bU(l) - guy(1,i)*hyi2
     bV(l) = bV(l) - gvy(1,i)*hyi2
    end if
    if (j.eq.ny-1) then
     bU(l) = bU(l) - guy(2,i)*hyi2
     bV(l) = bV(l) - gvy(2,i)*hyi2
    end if
   end do
  end do

  bU = -0.5*bU*nu
  bV = -0.5*bV*nu
  !
  FFu = u(1:nx-1,1:ny-1)*(1.0_dp/k) &
    + 1.5_dp*Leu(1:nx-1,1:ny-1) &
    - 0.5_dp*Leuold(1:nx-1,1:ny-1) &
    + 0.5_dp*Liu(1:nx-1,1:ny-1)
  FFv = v(1:nx-1,1:ny-1)*(1.0_dp/k) &
    + 1.5_dp*Lev(1:nx-1,1:ny-1) &
    - 0.5_dp*Levold(1:nx-1,1:ny-1)  &
    + 0.5_dp*Liv(1:nx-1,1:ny-1)
  !
  do j = 1,ny-1
   do i = 1,nx-1
    l = i + (j-1)*(nx-1)
    bU(l) = bU(l) + FFu(i,j)
    bV(l) = bV(l) + FFv(i,j)
   end do
  end do
end subroutine setupRhsideUV

subroutine taylor(u,v,X,Y,Lx,Ly,r0,gamma)
  use type_defs
  implicit none
  real(dp) :: u,v,x,y,lx,ly,r0,gamma,r

  R = sqrt((X-Lx)**2+(Y-Ly)**2)
  R = gamma*exp(-(R/r0)**2)
  u = -(Y-Ly)*R
  v =  (X-Lx)*R
end subroutine taylor


subroutine printdble2d(u,nx1,nx2,ny1,ny2,str)
  use type_defs
  implicit none
  integer, intent(in) :: nx1,nx2,ny1,ny2
  real(dp), intent(in) :: u(nx1:nx2,ny1:ny2)
  character(len=*), intent(in) :: str
  integer :: i,j
  open(2,file=trim(str),status='unknown')
  do j=ny1,ny2,1
   do i=nx1,nx2,1
    if(abs(u(i,j)) .lt. 1e-40) then
     write(2,fmt='(E24.16)',advance='no') 0.d0
    else
     write(2,fmt='(E24.16)',advance='no') u(i,j)
    end if
   end do
   write(2,'()')
  end do
  close(2)
end subroutine printdble2d
