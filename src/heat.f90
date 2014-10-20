! Heat equation solver (experimental)
subroutine heat(time,temp)
  use p3dfft_wrapper
  use basic_operators
  use vars
  implicit none
  ! Input/output
  type(timetype),intent(inout) :: time
  real(kind=pr),intent(inout) :: temp(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  ! Local variables
  integer :: ix,iy,iz,mpicommdir,mpiszdir,radir,rbdir,gadir,gbdir,mpirankdir,&
             mpicode
  real(kind=pr) :: xx,yy,zz,t1,norminf1,norminfloc1,norminf2,&
                   norminfloc2,norminf,norminfloc,&
                   tempdxdx,tempdydy,tempdzdz,dx2inv,dy2inv,dz2inv,h2inv,dt
  real(kind=pr) :: cnmaty(ra(2):rb(2),ra(2):rb(2)),vly(ra(2):rb(2)),vry(ra(2):rb(2)),&
                   utmpy(ga(2):gb(2)),bcmaty(2*mpidims(2),2*mpidims(2))
  real(kind=pr) :: cnmatz(ra(3):rb(3),ra(3):rb(3)),vlz(ra(3):rb(3)),vrz(ra(3):rb(3)),&
                   utmpz(ga(3):gb(3)),bcmatz(2*mpidims(1),2*mpidims(1))
  real(kind=pr) :: utmpx(ga(1):gb(1))

  ! This subroutine assumes that the domain decomposition is 2D
  ! The domain is NOT split in x direction
  ! In y direction, it is split in mpidims(2) parts
  ! In z direction, it is split in mpidims(1) parts

  ! Only the latest value of time step size is used
  dt = time%dt_new

  ! Y DIRECTION
  ! Set up parameters in y direction
  mpicommdir = mpicommy
  mpiszdir = mpidims(2)
  h2inv =  1.d0/(dy**2)
  radir = ra(2)
  rbdir = rb(2)
  gadir = ga(2)
  gbdir = gb(2)
  ! Cases if # subdomains = 1 or >=2
  if (mpiszdir>1) then 
    ! Parallel 1d solver init
    call heat_cn_1d_mpi_init (mpicommdir,mpiszdir,mpirankdir,h2inv,radir,rbdir,dt,&
                              bcmaty,cnmaty,vly,vry)
    ! Loop for all lines y=const
    do iz = ga(3),gb(3)
      !zz = dble(iz)*dz
      do ix = ga(1),gb(1)
        !xx = dble(ix)*dx 
        ! Vector to be processed
        utmpy(:) = temp(ix,gadir:gbdir,iz,1)
        ! Solve linear system
        call heat_cn_1d_mpi_solver (mpicommdir,mpiszdir,mpirankdir,h2inv,radir,rbdir,gadir,gbdir,dt,&
                                    bcmaty,cnmaty,vly,vry,utmpy)
        ! Vector returned
        temp(ix,radir:rbdir,iz,1) = utmpy(radir:rbdir)
      enddo 
    enddo
    ! Synchronize ghost points
    call synchronize_ghosts_FD_y_mpi_heat (temp(:,:,:,1))
  else
    do iz=ga(3),gb(3)
      !zz = dble(iz)*dz
      do ix=ga(1),gb(1)
        !xx = dble(ix)*dx 
        ! Vector to be processed
        utmpy(:) = temp(ix,gadir:gbdir,iz,1)
        ! Solve linear system
        call heat_cn_1d_serial_solver (h2inv,radir,rbdir,gadir,gbdir,dt,utmpy)
        ! Vector returned
        temp(ix,radir:rbdir,iz,1) = utmpy(radir:rbdir)
      enddo
    enddo    
    ! Synchronize ghost points
    call synchronize_ghosts_FD_y_serial_heat (temp(:,:,:,1))
  endif

  ! Z DIRECTION
  ! Set up parameters in z direction
  mpicommdir = mpicommz
  mpiszdir = mpidims(1)
  h2inv =  1.d0/(dz**2)
  radir = ra(3)
  rbdir = rb(3)
  gadir = ga(3)
  gbdir = gb(3)
  ! Cases if # subdomains = 1 or >=2
  if (mpiszdir>1) then 
    ! Parallel 1d solver init
    call heat_cn_1d_mpi_init (mpicommdir,mpiszdir,mpirankdir,h2inv,radir,rbdir,dt,&
                              bcmatz,cnmatz,vlz,vrz)
    ! Loop for all lines y=const
    do iy=ga(2),gb(2)
      !yy = dble(iy)*dy
      do ix=ga(1),gb(1)
        !xx = dble(ix)*dx 
        ! Vector to be processed
        utmpz(:) = temp(ix,iy,gadir:gbdir,1)
        call heat_cn_1d_mpi_solver (mpicommdir,mpiszdir,mpirankdir,h2inv,radir,rbdir,gadir,gbdir,dt,&
                                    bcmatz,cnmatz,vlz,vrz,utmpz)
        ! Vector returned
        temp(ix,iy,radir:rbdir,1) = utmpz(radir:rbdir)
      enddo 
    enddo
    ! Synchronize ghost points
    call synchronize_ghosts_FD_z_mpi_heat (temp(:,:,:,1))
  else
    do iy=ga(2),gb(2)
      !yy = dble(iy)*dy
      do ix=ga(1),gb(1)
        !xx = dble(ix)*dx 
        ! Vector to be processed
        utmpz(:) = temp(ix,iy,gadir:gbdir,1)
        ! Solve linear system
        call heat_cn_1d_serial_solver (h2inv,radir,rbdir,gadir,gbdir,dt,utmpz)
        ! Vector returned
        temp(ix,iy,radir:rbdir,1) = utmpz(radir:rbdir)
      enddo
    enddo    
    ! Synchronize ghost points
    call synchronize_ghosts_FD_z_serial_heat (temp(:,:,:,1))
  endif

  ! X DIRECTION
  ! Set up parameters in x direction
  h2inv =  1.d0/(dx**2)
  radir = ra(1)
  rbdir = rb(1)
  gadir = ga(1)
  gbdir = gb(1)
  ! Loop for all lines x=const. This is local.
  do iz=ga(3),gb(3)
    !zz = dble(iz)*dz
    do iy=ga(2),gb(2)
      !xx = dble(ix)*dx 
      ! Vector to be processed
      utmpx(:) = temp(gadir:gbdir,iy,iz,1)
      ! Solve linear system
      call heat_cn_1d_serial_solver (h2inv,radir,rbdir,gadir,gbdir,dt,utmpx)
      ! Vector returned
      temp(radir:rbdir,iy,iz,1) = utmpx(radir:rbdir)
    enddo
  enddo    
  ! Synchronize ghost points
  call synchronize_ghosts_FD_x_serial_heat (temp(:,:,:,1))

!!!!!!!!!!!!!!!!!!!!!
!  call synchronize_ghosts_FD_x_heat (temp(:,:,:,1))
!  if (mpidims(2)>1) then
!    call synchronize_ghosts_FD_y_heat (temp(:,:,:,1))
!  else
!    call synchronize_ghosts_FD_y_serial_heat (temp(:,:,:,1))
!  endif
!  if (mpidims(1)>1) then
!    call synchronize_ghosts_FD_z_heat (temp(:,:,:,1))
!  else
!    call synchronize_ghosts_FD_z_serial_heat (temp(:,:,:,1))
!  endif
!!!!!!!!!!!!!!!!!!!!!

  norminfloc = 0.0d0
  norminfloc1 = 0.0d0
  norminfloc2 = 0.0d0
  do iz=ra(3),rb(3)
    zz = dble(iz)*dz
    do iy=ra(2),rb(2)
      yy = dble(iy)*dy
      do ix=ra(1),rb(1)
        xx = dble(ix)*dx 

        temp(ix,iy,iz,2) = dexp(-2*pi**2*nu*time%time)*(dexp(-15*pi**2*nu*time%time)*dcos(4*pi*yy)+dcos(pi*yy))*dcos(pi*zz)

        norminfloc = max(norminfloc,abs(temp(ix,iy,iz,1)-temp(ix,iy,iz,2)))
        norminfloc1 = max(norminfloc1,abs(temp(ix,iy,iz,1)))
        norminfloc2 = max(norminfloc2,abs(temp(ix,iy,iz,2)))
      enddo
    enddo
  enddo    
  call MPI_ALLREDUCE (norminfloc,norminf,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,mpicode) 
  call MPI_ALLREDUCE (norminfloc1,norminf1,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,mpicode) 
  call MPI_ALLREDUCE (norminfloc2,norminf2,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,mpicode) 

  if (mpirank==0) then 
    open(14,file='heatnorm.t',status='unknown',position='append')
    write(14,'(i,5(es15.8,1x))') time%it,time%time,nu,norminf1,norminf2,norminf
    print *, 't=',time%time,'norminf=',norminf,'norminf1=',norminf1,'norminf2=',norminf2
    close(14)
  endif
end subroutine heat

subroutine heat_init(temp)
  use p3dfft_wrapper
  use basic_operators
  use vars
  implicit none

  real(kind=pr),intent(inout)::temp(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:2)
  real(kind=pr)::xx,yy,zz
  integer::ix,iy,iz

  ! for linear solver
!  integer :: nn
!  real(kind=pr),allocatable::cnmat(:,:)
!  real(kind=pr),allocatable::rhs(:),vf(:)

  do iz=ra(3),rb(3)
    zz = dble(iz)*dz
    do iy=ra(2),rb(2)
      yy = dble(iy)*dy
      do ix=ra(1),rb(1)
        xx = dble(ix)*dx 
        temp(ix,iy,iz,1) = (dcos(4*pi*yy)+dcos(pi*yy))*dcos(pi*zz)
        temp(ix,iy,iz,2) = 0.0
      enddo
    enddo
  enddo     
  
  ! Test linear solver
!  if (mpirank==0) then
!    allocate(cnmat(ra(1):rb(1),ra(1):rb(1)))
!    allocate(rhs(ra(1):rb(1)))
!    allocate(vf(ra(1):rb(1)))
!    nn = rb(1)-ra(1)+1
!    cnmat = 0.0d0
!    do ix=ra(1),rb(1)
!      cnmat(ix,ix) = 1.0d0
!      rhs(ix) = dble(ix)
!    enddo
!    call solve_loc1d ( cnmat, rhs, vf, nn )
!    print *, vf
!    deallocate(cnmat,rhs,vf)
!  endif

  ! Setup line communicators
  call setup_cart_groups

  ! Synchronize ghost points
  call synchronize_ghosts_FD_heat (temp(:,:,:,1))
end subroutine heat_init

subroutine solve_loc1d ( mat, rhs, x, nn )
  !--------------------------------------------
  ! solves the linear system J*x = F
  !--------------------------------------------
  use vars
  implicit none
  integer, intent(in) :: nn
  real(kind=pr),dimension(1:nn,1:nn), intent(in) :: mat
  real(kind=pr),dimension(1:nn), intent(out) :: x
  real(kind=pr),dimension(1:nn), intent(in) :: rhs
  real(kind=pr),dimension(1:nn,1:nn) :: mat2
  real(kind=pr) :: t0
  integer :: error,ipiv(1:nn)  
  t0 = MPI_wtime()
  
  !mat2 = transpose(mat)
  mat2 = mat
  call dgetrf (nn,nn,mat2,nn,ipiv,error)
  if (error .ne. 0) then
    write(*,*) "!!! Crutial: dgetrf error.", error
    call abort()
  endif
  
  x = rhs
  call dgetrs ('N',nn,1,mat2,nn,ipiv,x,nn,error)
  if (error .ne. 0) then 
    write(*,*) "!!! Crutial: dgetrs error.", error
    call abort()
  endif
  
  time_LAPACK = time_LAPACK + MPI_wtime() - t0
end subroutine solve_loc1d


! LOD splitting. Initialization of the 1d implicit MPI solver
subroutine heat_cn_1d_mpi_init(mpicommdir,mpiszdir,mpirankdir,h2inv,radir,rbdir,dt,&
                           bcmat,cnmat,vl,vr)
  use p3dfft_wrapper
  use basic_operators
  use vars
  implicit none
  ! Input/output
  integer,intent(inout) :: mpicommdir,mpiszdir,mpirankdir,radir,rbdir
  real(kind=pr),intent(inout) :: h2inv,dt
  real(kind=pr),intent(inout) :: cnmat(radir:rbdir,radir:rbdir),bcmat(2*mpiszdir,2*mpiszdir),&
                                 vl(radir:rbdir),vr(radir:rbdir)
  ! Local variables
  integer :: nn,j,mpicode
  real(kind=pr) :: sendfoo(4),recvfoo(4*mpiszdir)
  real(kind=pr) :: vl1(mpiszdir),vlN(mpiszdir),vr1(mpiszdir),vrN(mpiszdir)
  real(kind=pr) :: rhs(radir:rbdir)
  ! Get local ranks in the line
  call MPI_COMM_RANK(mpicommdir,mpirankdir,mpicode)
  ! Crank-Nicolson matrix in x direction
  cnmat(:,:) = 0.d0
  do j = radir,rbdir
    cnmat(j,j) = 1.d0 + 1.d0*dt*nu*h2inv
  enddo
  do j = radir,rbdir-1
    cnmat(j,j+1) = - 0.5d0*dt*nu*h2inv
    cnmat(j+1,j) = - 0.5d0*dt*nu*h2inv
  enddo
  ! Boundary conditions for domain decomposition
  ! BC influence basis
  nn = rbdir-radir+1
  rhs(:) = 0.d0
  rhs(radir) = 1.d0
  call solve_loc1d (cnmat,rhs,vl,nn)
  vl(:) = (-0.5d0*dt*nu*h2inv)*vl(:)
  rhs(rbdir) = 1.d0
  rhs(radir) = 0.d0
  call solve_loc1d (cnmat,rhs,vr,nn)
  vr(:) = (-0.5d0*dt*nu*h2inv)*vr(:)
  ! BC influence matrix
  ! It is only stored by one process
  ! Communicate values at the interface to rank 0
  sendfoo(1) = vl(radir)
  sendfoo(2) = vl(rbdir)
  sendfoo(3) = vr(radir)
  sendfoo(4) = vr(rbdir)
  call MPI_GATHER (sendfoo,4,MPI_DOUBLE_PRECISION,recvfoo,4,MPI_DOUBLE_PRECISION,0,mpicommdir,mpicode) 
  do j = 1,mpiszdir
    vl1(j) = recvfoo(4*j-3)
    vlN(j) = recvfoo(4*j-2)
    vr1(j) = recvfoo(4*j-1)
    vrN(j) = recvfoo(4*j)
  enddo
  ! BC influence matrix is only stored by one process
  if (mpirankdir == 0) then
    bcmat(:,:) = 0.d0
    do j = 1,mpiszdir
        bcmat(2*j-1,2*j-1) = 1.d0
        bcmat(2*j-1,modulo(2*j-3,2*mpiszdir)+1) = vl1(j)
        bcmat(2*j-1,modulo(2*j,2*mpiszdir)+1) = vr1(j)
        bcmat(2*j,2*j) = 1.d0
        bcmat(2*j,modulo(2*j-3,2*mpiszdir)+1) = vlN(j)
        bcmat(2*j,modulo(2*j,2*mpiszdir)+1) = vrN(j)
    enddo   
  endif
end subroutine heat_cn_1d_mpi_init


! LOD splitting. 1d MPI solver
subroutine heat_cn_1d_mpi_solver(mpicommdir,mpiszdir,mpirankdir,h2inv,radir,rbdir,gadir,gbdir,dt,&
                             bcmat,cnmat,vl,vr,utmp)
  use p3dfft_wrapper
  use basic_operators
  use vars
  implicit none
  ! Input/output
  integer,intent(inout) :: mpicommdir,mpiszdir,mpirankdir,radir,rbdir,gadir,gbdir
  real(kind=pr),intent(inout) :: h2inv,dt
  real(kind=pr),intent(inout) :: cnmat(radir:rbdir,radir:rbdir),bcmat(2*mpiszdir,2*mpiszdir),&
                                 vl(radir:rbdir),vr(radir:rbdir),utmp(gadir:gbdir)
  ! local variables
  integer :: j,mpicode
  real(kind=pr) :: bcxl,bcxr,shortfoo(2*mpiszdir),longfoo(4*mpiszdir)
  real(kind=pr) :: bcrhs(2*mpiszdir),bcx(2*mpiszdir),&
                   rhs(radir:rbdir),vf(radir:rbdir),bcxls(mpiszdir),&
                   bcxrs(mpiszdir),vf1(mpiszdir),vfN(mpiszdir)
  
  ! Crank-Nicolson explicit part
  rhs(:) = utmp(radir:rbdir)+0.5d0*dt*nu*(utmp((radir-1):(rbdir-1))-2.d0*utmp(radir:rbdir)+utmp((radir+1):(rbdir+1)))*h2inv
  ! Solve local system
  call solve_loc1d (cnmat,rhs,vf,rbdir-radir+1)
  ! Communicate rhs to rank 0 in the line
  shortfoo(1) = vf(radir)
  shortfoo(2) = vf(rbdir)
  call MPI_GATHER (shortfoo,2,MPI_DOUBLE_PRECISION,longfoo,2,MPI_DOUBLE_PRECISION,0,mpicommdir,mpicode) 
  do j = 1,mpiszdir
    vf1(j) = longfoo(2*j-1)
    vfN(j) = longfoo(2*j)
  enddo
  ! BC influence RHS
  if (mpirankdir == 0) then
    do j = 1,mpiszdir
      bcrhs(2*j-1) = vf1(j)
      bcrhs(2*j) = vfN(j)
    enddo
    ! Solve BC influence system
    call solve_loc1d (bcmat,bcrhs,bcx,2*mpiszdir)
    ! Rearrange for mpi scatter
    do j = 1,mpiszdir
      bcxls(j) = bcx(modulo(2*j-3,2*mpiszdir)+1)
      bcxrs(j) = bcx(modulo(2*j,2*mpiszdir)+1)
    enddo
  endif
  ! Scatter from rank 0 in the line to all ranks
  do j = 1,mpiszdir
    longfoo(2*j-1) = bcxls(j)
    longfoo(2*j) = bcxrs(j)
  enddo
  call MPI_SCATTER (longfoo,2,MPI_DOUBLE_PRECISION,shortfoo,2,MPI_DOUBLE_PRECISION,0,mpicommdir,mpicode) 
  bcxl = shortfoo(1)
  bcxr = shortfoo(2)
  ! Superpose local solution and BC influence
  utmp(radir:rbdir) = vf(:)-bcxl*vl(:)-bcxr*vr(:)
end subroutine heat_cn_1d_mpi_solver


! LOD splitting. 1d serial solver
subroutine heat_cn_1d_serial_solver(h2inv,radir,rbdir,gadir,gbdir,dt,utmp)
  use p3dfft_wrapper
  use basic_operators
  use vars
  implicit none
  ! Input/output
  integer,intent(inout) :: radir,rbdir,gadir,gbdir
  real(kind=pr),intent(inout) :: h2inv,dt
  real(kind=pr),intent(inout) :: utmp(gadir:gbdir)
  ! local variables
  integer :: j
  real(kind=pr) :: cnmat(radir:rbdir,radir:rbdir),rhs(radir:rbdir)

  ! Crank-Nicolson matrix in x direction
  cnmat(:,:) = 0.d0
  do j = radir,rbdir
    cnmat(j,j) = 1.d0 + 1.d0*dt*nu*h2inv
  enddo
  do j = radir,rbdir-1
    cnmat(j,j+1) = - 0.5d0*dt*nu*h2inv
    cnmat(j+1,j) = - 0.5d0*dt*nu*h2inv
  enddo
  ! This matrix is circulant
  cnmat(radir,rbdir) = - 0.5d0*dt*nu*h2inv
  cnmat(rbdir,radir) = - 0.5d0*dt*nu*h2inv
  ! Crank-Nicolson explicit part
  rhs(:) = utmp(radir:rbdir)+0.5d0*dt*nu*(utmp((radir-1):(rbdir-1))-2.d0*utmp(radir:rbdir)+utmp((radir+1):(rbdir+1)))*h2inv
  ! Solve local system
  call solve_loc1d (cnmat,rhs,utmp(radir:rbdir),rbdir-radir+1)
end subroutine heat_cn_1d_serial_solver

