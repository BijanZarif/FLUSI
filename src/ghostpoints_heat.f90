! Setup 1d communicators
subroutine setup_cart_groups
  use p3dfft_wrapper
  use basic_operators
  use vars
  implicit none
  integer :: mpicolor,mpikey,mpicode
  integer :: mpicommtmp1,mpicommtmp2
  logical :: mpiperiods(2),period,reorder
  ! Set parameters
  period=.true.
  reorder=.false.
  ! Get Cartesian topology information
  call MPI_CART_GET(mpicommcart,2,mpidims,mpiperiods,mpicoords,mpicode)
  ! Communicator for line in y direction
  mpicolor = mpicoords(2) 
  mpikey = mpicoords(1)
  call MPI_COMM_SPLIT (mpicommcart,mpicolor,mpikey,mpicommtmp1,mpicode)
  call MPI_CART_CREATE(mpicommtmp1,1,mpidims(1),period,reorder,mpicommz,mpicode)
  ! Communicator for line in z direction
  mpicolor = mpicoords(1) 
  mpikey = mpicoords(2)
  call MPI_COMM_SPLIT (mpicommcart,mpicolor,mpikey,mpicommtmp2,mpicode)
  call MPI_CART_CREATE(mpicommtmp2,1,mpidims(2),period,reorder,mpicommy,mpicode)
end subroutine setup_cart_groups

! Ghost point synchronization in z
subroutine synchronize_ghosts_FD_z_mpi_heat(fld)
  use p3dfft_wrapper
  use basic_operators
  use vars
  implicit none
  ! Input/output
  real(kind=pr),intent(inout)::fld(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  ! Local variables
  integer :: statu(MPI_STATUS_SIZE),nnl,nnr,disp,dir,source,dest,mpicode
  real(kind=pr) :: fld_send_l(ga(1):gb(1),ga(2):gb(2),ra(3):(2*ra(3)-ga(3)-1)),&
                   fld_send_r(ga(1):gb(1),ga(2):gb(2),(2*rb(3)-gb(3)+1):rb(3))
  real(kind=pr) :: fld_recv_l(ga(1):gb(1),ga(2):gb(2),ga(3):(ra(3)-1)),&
                   fld_recv_r(ga(1):gb(1),ga(2):gb(2),(rb(3)+1):gb(3))
  ! Size of buffer arrays
  nnl = (ra(3)-ga(3))*(gb(2)-ga(2)+1)*(gb(1)-ga(1)+1)
  nnr = (gb(3)-rb(3))*(gb(2)-ga(2)+1)*(gb(1)-ga(1)+1)
  ! Copy data to buffer arrays
  fld_send_l(:,:,:) = fld(ga(1):gb(1),ga(2):gb(2),ra(3):(2*ra(3)-ga(3)-1))
  fld_send_r(:,:,:) = fld(ga(1):gb(1),ga(2):gb(2),(2*rb(3)-gb(3)+1):rb(3))
  ! Find rank of neighbors on the right
  disp=1                 !immediate neighbors
  dir=0
  call MPI_CART_SHIFT(mpicommz,dir,disp,source,dest,mpicode)
  ! Send data to the neighbors on the right
  call MPI_SENDRECV(fld_send_r,nnr,MPI_DOUBLE_PRECISION,dest,0,&
                    fld_recv_l,nnr,MPI_DOUBLE_PRECISION,source,0,&
                    mpicommz,statu,mpicode)
  ! Find rank of neighbors on the left
  disp=-1                !immediate neighbors
  dir=0
  call MPI_CART_SHIFT(mpicommz,dir,disp,source,dest,mpicode)
  ! Send data to the neighbors on the left
  call MPI_SENDRECV(fld_send_l,nnl,MPI_DOUBLE_PRECISION,dest,0,&
                    fld_recv_r,nnl,MPI_DOUBLE_PRECISION,source,0,&
                    mpicommz,statu,mpicode)
  ! Copy data to output array
  fld(ga(1):gb(1),ga(2):gb(2),ga(3):(ra(3)-1)) = fld_recv_l(:,:,:)
  fld(ga(1):gb(1),ga(2):gb(2),(rb(3)+1):gb(3)) = fld_recv_r(:,:,:)
end subroutine synchronize_ghosts_FD_z_mpi_heat

! Ghost point synchronization in y
subroutine synchronize_ghosts_FD_y_mpi_heat(fld)
  use p3dfft_wrapper
  use basic_operators
  use vars
  implicit none
  ! Input/output
  real(kind=pr),intent(inout)::fld(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  ! Local variables
  integer :: statu(MPI_STATUS_SIZE),nnl,nnr,disp,dir,source,dest,mpicode
  real(kind=pr) :: fld_send_l(ga(1):gb(1),ra(2):(2*ra(2)-ga(2)-1),ga(3):gb(3)),&
                   fld_send_r(ga(1):gb(1),(2*rb(2)-gb(2)+1):rb(2),ga(3):gb(3))
  real(kind=pr) :: fld_recv_l(ga(1):gb(1),ga(2):(ra(2)-1),ga(3):gb(3)),&
                   fld_recv_r(ga(1):gb(1),(rb(2)+1):gb(2),ga(3):gb(3))
  ! Size of buffer arrays
  nnl = (ra(2)-ga(2))*(gb(3)-ga(3)+1)*(gb(1)-ga(1)+1)
  nnr = (gb(2)-rb(2))*(gb(3)-ga(3)+1)*(gb(1)-ga(1)+1)
!print *,mpirank,nnl,nnr,size(fld_send_l),size(fld_send_r),size(fld_recv_l),size(fld_recv_r)
  ! Copy data to buffer arrays
  fld_send_l(:,:,:) = fld(ga(1):gb(1),ra(2):(2*ra(2)-ga(2)-1),ga(3):gb(3))
  fld_send_r(:,:,:) = fld(ga(1):gb(1),(2*rb(2)-gb(2)+1):rb(2),ga(3):gb(3))

  ! Find rank of neighbors on the right
  disp=1                 !immediate neighbors
  dir=0
  call MPI_CART_SHIFT(mpicommy,dir,disp,source,dest,mpicode)
  ! Send data to the neighbors on the right
  call MPI_SENDRECV(fld_send_r,nnr,MPI_DOUBLE_PRECISION,dest,0,&
                    fld_recv_l,nnr,MPI_DOUBLE_PRECISION,source,0,&
                    mpicommy,statu,mpicode)
  ! Find rank of neighbors on the left
  disp=-1                !immediate neighbors
  dir=0
  call MPI_CART_SHIFT(mpicommy,dir,disp,source,dest,mpicode)
  ! Send data to the neighbors on the left
  call MPI_SENDRECV(fld_send_l,nnl,MPI_DOUBLE_PRECISION,dest,0,&
                    fld_recv_r,nnl,MPI_DOUBLE_PRECISION,source,0,&
                    mpicommy,statu,mpicode)
  ! Copy data to output array
  fld(ga(1):gb(1),ga(2):(ra(2)-1),ga(3):gb(3)) = fld_recv_l(:,:,:)
  fld(ga(1):gb(1),(rb(2)+1):gb(2),ga(3):gb(3)) = fld_recv_r(:,:,:)
end subroutine synchronize_ghosts_FD_y_mpi_heat

! Ghost point synchronization in x. This is always local
subroutine synchronize_ghosts_FD_x_serial_heat(fld)
  use p3dfft_wrapper
  use basic_operators
  use vars
  implicit none
  ! Input/output
  real(kind=pr),intent(inout)::fld(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  ! Copy data to ghost points. X direction is local
  fld((rb(1)+1):gb(1),ga(2):gb(2),ga(3):gb(3)) = fld(ra(1):(2*ra(1)-ga(1)-1),ga(2):gb(2),ga(3):gb(3))
  fld(ga(1):(ra(1)-1),ga(2):gb(2),ga(3):gb(3)) = fld((2*rb(1)-gb(1)+1):rb(1),ga(2):gb(2),ga(3):gb(3))
end subroutine synchronize_ghosts_FD_x_serial_heat

! Ghost point synchronization in y on one process
subroutine synchronize_ghosts_FD_y_serial_heat(fld)
  use p3dfft_wrapper
  use basic_operators
  use vars
  implicit none
  ! Input/output
  real(kind=pr),intent(inout)::fld(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  ! Copy data to ghost points. Y direction is local in this case
  fld(ga(1):gb(1),(rb(2)+1):gb(2),ga(3):gb(3)) = fld(ga(1):gb(1),ra(2):(2*ra(2)-ga(2)-1),ga(3):gb(3))
  fld(ga(1):gb(1),ga(2):(ra(2)-1),ga(3):gb(3)) = fld(ga(1):gb(1),(2*rb(2)-gb(2)+1):rb(2),ga(3):gb(3))
end subroutine synchronize_ghosts_FD_y_serial_heat

! Ghost point synchronization in z on one process
subroutine synchronize_ghosts_FD_z_serial_heat(fld)
  use p3dfft_wrapper
  use basic_operators
  use vars
  implicit none
  ! Input/output
  real(kind=pr),intent(inout)::fld(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  ! Copy data to ghost points. Z direction is local in this case
  fld(ga(1):gb(1),ga(2):gb(2),(rb(3)+1):gb(3)) = fld(ga(1):gb(1),ga(2):gb(2),ra(3):(2*ra(3)-ga(3)-1))
  fld(ga(1):gb(1),ga(2):gb(2),ga(3):(ra(3)-1)) = fld(ga(1):gb(1),ga(2):gb(2),(2*rb(3)-gb(3)+1):rb(3))
end subroutine synchronize_ghosts_FD_z_serial_heat

! Ghost point synchronization in all directions
subroutine synchronize_ghosts_FD_heat(fld)
  use p3dfft_wrapper
  use basic_operators
  use vars
  implicit none
  ! Input/output
  real(kind=pr),intent(inout)::fld(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))

  call synchronize_ghosts_FD_x_serial_heat (fld(:,:,:))
  if (mpidims(2)>1) then
    call synchronize_ghosts_FD_y_mpi_heat (fld(:,:,:))
  else
    call synchronize_ghosts_FD_y_serial_heat (fld(:,:,:))
  endif
  if (mpidims(1)>1) then
    call synchronize_ghosts_FD_z_mpi_heat (fld(:,:,:))
  else
    call synchronize_ghosts_FD_z_serial_heat (fld(:,:,:))
  endif
end subroutine synchronize_ghosts_FD_heat









