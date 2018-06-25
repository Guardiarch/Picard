

subroutine CalcEfield(E,U,Nx_local,Ny_local,Nz_local,dxyz,xflow,yflow,zflow,iprocs,jprocs,kprocs,myid)
  use mpi
  implicit none


  logical, intent(in) :: xflow,yflow,zflow
  integer, intent(in) :: Nx_local,Ny_local,Nz_local,iprocs,jprocs,kprocs,myid
  real*8, intent(in) :: U(Nx_local+2,Ny_local+2,Nz_local+2), dxyz(3)
  real*8, intent(out) :: E(3,Nx_local+2,Ny_local+2,Nz_local+2)

! local variables
  integer :: ii, idpx, idnx, idpy, idny, idpz, idnz, tag, ierr, req(4), size_sendrecv
  real*8, allocatable :: real_neg_send(:,:,:), real_neg_recv(:,:,:), real_pos_send(:,:,:), real_pos_recv(:,:,:)

  
  idpx = myid + (mod(myid+iprocs+1,                           iprocs) - mod(myid+iprocs, iprocs))
  idnx = myid + (mod(myid+iprocs-1,                           iprocs) - mod(myid+iprocs, iprocs))
  idpy = myid + (mod(myid+iprocs*jprocs+iprocs,               iprocs*jprocs) - mod(myid+iprocs*jprocs, iprocs*jprocs))
  idny = myid + (mod(myid+iprocs*jprocs-iprocs,               iprocs*jprocs) - mod(myid+iprocs*jprocs, iprocs*jprocs))
  idpz = myid + (mod(myid+iprocs*jprocs*kprocs+iprocs*jprocs, iprocs*jprocs*kprocs) &
               - mod(myid+iprocs*jprocs*kprocs,               iprocs*jprocs*kprocs))
  idnz = myid + (mod(myid+iprocs*jprocs*kprocs-iprocs*jprocs, iprocs*jprocs*kprocs) &
               - mod(myid+iprocs*jprocs*kprocs,               iprocs*jprocs*kprocs))

  E = 0.0d0

  E(3,2:Nx_local+1,2:Ny_local+1,2:Nz_local+1) = -( U(2:Nx_local+1,2:Ny_local+1,3:Nz_local+2) &
                                                  -U(2:Nx_local+1,2:Ny_local+1,1:Nz_local+0) )*0.5d0/dxyz(3)
  E(2,2:Nx_local+1,2:Ny_local+1,2:Nz_local+1) = -( U(2:Nx_local+1,3:Ny_local+2,2:Nz_local+1) &
                                                  -U(2:Nx_local+1,1:Ny_local+0,2:Nz_local+1) )*0.5d0/dxyz(2)
  E(1,2:Nx_local+1,2:Ny_local+1,2:Nz_local+1) = -( U(3:Nx_local+2,2:Ny_local+1,2:Nz_local+1) &
                                                  -U(1:Nx_local+0,2:Ny_local+1,2:Nz_local+1) )*0.5d0/dxyz(1)

!  E(3,2:Nx_local+1,2:Ny_local+1,2:Nz_local+1) = -( 4.0d0*U(2:Nx_local+1,2:Ny_local+1,3:Nz_local+2) &
!                                                 + 2.0d0*U(2:Nx_local+1,3:Ny_local+2,3:Nz_local+2) &
!                                                 + 2.0d0*U(3:Nx_local+2,2:Ny_local+1,3:Nz_local+2) &
!                                                 + 2.0d0*U(2:Nx_local+1,1:Ny_local+0,3:Nz_local+2) &
!                                                 + 2.0d0*U(1:Nx_local+0,2:Ny_local+1,3:Nz_local+2) &
!                                                 + 1.0d0*U(3:Nx_local+2,3:Ny_local+2,3:Nz_local+2) &
!                                                 + 1.0d0*U(3:Nx_local+2,1:Ny_local+0,3:Nz_local+2) &
!                                                 + 1.0d0*U(1:Nx_local+0,3:Ny_local+2,3:Nz_local+2) &
!                                                 + 1.0d0*U(1:Nx_local+0,1:Ny_local+0,3:Nz_local+2) &
!                                                 - 4.0d0*U(2:Nx_local+1,2:Ny_local+1,1:Nz_local+0) &
!                                                 - 2.0d0*U(2:Nx_local+1,3:Ny_local+2,1:Nz_local+0) &
!                                                 - 2.0d0*U(3:Nx_local+2,2:Ny_local+1,1:Nz_local+0) &
!                                                 - 2.0d0*U(2:Nx_local+1,1:Ny_local+0,1:Nz_local+0) &
!                                                 - 2.0d0*U(1:Nx_local+0,2:Ny_local+1,1:Nz_local+0) &
!                                                 - 1.0d0*U(3:Nx_local+2,3:Ny_local+2,1:Nz_local+0) &
!                                                 - 1.0d0*U(3:Nx_local+2,1:Ny_local+0,1:Nz_local+0) &
!                                                 - 1.0d0*U(1:Nx_local+0,3:Ny_local+2,1:Nz_local+0) &
!                                                 - 1.0d0*U(1:Nx_local+0,1:Ny_local+0,1:Nz_local+0) )*0.03125d0*dxyz(3)**(-1)
!  E(2,2:Nx_local+1,2:Ny_local+1,2:Nz_local+1) = -( 4.0d0*U(2:Nx_local+1,3:Ny_local+2,2:Nz_local+1) &
!                                                 + 2.0d0*U(2:Nx_local+1,3:Ny_local+2,3:Nz_local+2) &
!                                                 + 2.0d0*U(2:Nx_local+1,3:Ny_local+2,1:Nz_local+0) &
!                                                 + 2.0d0*U(3:Nx_local+2,3:Ny_local+2,2:Nz_local+1) &
!                                                 + 2.0d0*U(1:Nx_local+0,3:Ny_local+2,2:Nz_local+1) &
!                                                 + 1.0d0*U(3:Nx_local+2,3:Ny_local+2,3:Nz_local+2) &
!                                                 + 1.0d0*U(3:Nx_local+2,3:Ny_local+2,1:Nz_local+0) &
!                                                 + 1.0d0*U(1:Nx_local+0,3:Ny_local+2,3:Nz_local+2) &
!                                                 + 1.0d0*U(1:Nx_local+0,3:Ny_local+2,1:Nz_local+0) &
!                                                 - 4.0d0*U(2:Nx_local+1,1:Ny_local+0,2:Nz_local+1) &
!                                                 - 2.0d0*U(2:Nx_local+1,1:Ny_local+0,3:Nz_local+2) &
!                                                 - 2.0d0*U(2:Nx_local+1,1:Ny_local+0,1:Nz_local+0) &
!                                                 - 2.0d0*U(3:Nx_local+2,1:Ny_local+0,2:Nz_local+1) &
!                                                 - 2.0d0*U(1:Nx_local+0,1:Ny_local+0,2:Nz_local+1) &
!                                                 - 1.0d0*U(3:Nx_local+2,1:Ny_local+0,3:Nz_local+2) &
!                                                 - 1.0d0*U(3:Nx_local+2,1:Ny_local+0,1:Nz_local+0) &
!                                                 - 1.0d0*U(1:Nx_local+0,1:Ny_local+0,3:Nz_local+2) &
!                                                 - 1.0d0*U(1:Nx_local+0,1:Ny_local+0,1:Nz_local+0) )*0.03125d0*dxyz(2)**(-1)
!  E(1,2:Nx_local+1,2:Ny_local+1,2:Nz_local+1) = -( 4.0d0*U(3:Nx_local+2,2:Ny_local+1,2:Nz_local+1) &
!                                                 + 2.0d0*U(3:Nx_local+2,3:Ny_local+2,2:Nz_local+1) &
!                                                 + 2.0d0*U(3:Nx_local+2,2:Ny_local+1,3:Nz_local+2) &
!                                                 + 2.0d0*U(3:Nx_local+2,1:Ny_local+0,2:Nz_local+1) &
!                                                 + 2.0d0*U(3:Nx_local+2,2:Ny_local+1,1:Nz_local+0) &
!                                                 + 1.0d0*U(3:Nx_local+2,3:Ny_local+2,3:Nz_local+2) &
!                                                 + 1.0d0*U(3:Nx_local+2,1:Ny_local+0,3:Nz_local+2) &
!                                                 + 1.0d0*U(3:Nx_local+2,3:Ny_local+2,1:Nz_local+0) &
!                                                 + 1.0d0*U(3:Nx_local+2,1:Ny_local+0,1:Nz_local+0) &
!                                                 - 4.0d0*U(1:Nx_local+0,2:Ny_local+1,2:Nz_local+1) &
!                                                 - 2.0d0*U(1:Nx_local+0,3:Ny_local+2,2:Nz_local+1) &
!                                                 - 2.0d0*U(1:Nx_local+0,2:Ny_local+1,3:Nz_local+2) &
!                                                 - 2.0d0*U(1:Nx_local+0,1:Ny_local+0,2:Nz_local+1) &
!                                                 - 2.0d0*U(1:Nx_local+0,2:Ny_local+1,1:Nz_local+0) &
!                                                 - 1.0d0*U(1:Nx_local+0,3:Ny_local+2,3:Nz_local+2) &
!                                                 - 1.0d0*U(1:Nx_local+0,1:Ny_local+0,3:Nz_local+2) &
!                                                 - 1.0d0*U(1:Nx_local+0,3:Ny_local+2,1:Nz_local+0) &
!                                                 - 1.0d0*U(1:Nx_local+0,1:Ny_local+0,1:Nz_local+0) )*0.03125d0*dxyz(1)**(-1)


! Z DIRECTION

  allocate( real_neg_send(3,  Nx_local+2,  Ny_local+2) )
  allocate( real_pos_send(3,  Nx_local+2,  Ny_local+2) )
  allocate( real_neg_recv(3,  Nx_local+2,  Ny_local+2) )
  allocate( real_pos_recv(3,  Nx_local+2,  Ny_local+2) )
  size_sendrecv        = (3)*(Nx_local+2)*(Ny_local+2)
  real_neg_send(:, :, :) = E(:, :, :, 2)
  real_pos_send(:, :, :) = E(:, :, :, Nz_local+1)

  
  tag = 1
  call MPI_IRECV(real_pos_recv, size_sendrecv, MPI_DOUBLE_PRECISION, idpz, tag, MPI_COMM_WORLD, req(1), ierr)
  call MPI_ISEND(real_neg_send, size_sendrecv, MPI_DOUBLE_PRECISION, idnz, tag, MPI_COMM_WORLD, req(3), ierr)
  tag = 2
  call MPI_IRECV(real_neg_recv, size_sendrecv, MPI_DOUBLE_PRECISION, idnz, tag, MPI_COMM_WORLD, req(2), ierr)
  call MPI_ISEND(real_pos_send, size_sendrecv, MPI_DOUBLE_PRECISION, idpz, tag, MPI_COMM_WORLD, req(4), ierr)

! Wait to receive
  call MPI_WAITALL(2, req(1:2), MPI_STATUSES_IGNORE, ierr)

  
  if (zflow) then

    if (myid - idnz .gt. 0) then
      E(:,:,:,1)          = real_neg_recv(:,:,:)
    end if

    if (idpz - myid .gt. 0) then
      E(:,:,:,Nz_local+2) = real_pos_recv(:,:,:)
    end if
    
  else

  ! periodic
    E(:,:,:,1)            = real_neg_recv(:,:,:)
    E(:,:,:,Nz_local+2)   = real_pos_recv(:,:,:)
    
  end if
 
  call MPI_WAITALL(2, req(3:4), MPI_STATUSES_IGNORE, ierr)

  deallocate( real_neg_send )
  deallocate( real_pos_send )
  deallocate( real_neg_recv )
  deallocate( real_pos_recv )


! Y DIRECTION

  allocate( real_neg_send(3,  Nx_local+2,  Nz_local+2) )
  allocate( real_pos_send(3,  Nx_local+2,  Nz_local+2) )
  allocate( real_neg_recv(3,  Nx_local+2,  Nz_local+2) )
  allocate( real_pos_recv(3,  Nx_local+2,  Nz_local+2) )
  size_sendrecv        = (3)*(Nx_local+2)*(Nz_local+2)
  real_neg_send(:, :, :) = E(:, :, 2, :)
  real_pos_send(:, :, :) = E(:, :, Ny_local+1, :)

  
  tag = 1
  call MPI_IRECV(real_pos_recv, size_sendrecv, MPI_DOUBLE_PRECISION, idpy, tag, MPI_COMM_WORLD, req(1), ierr)
  call MPI_ISEND(real_neg_send, size_sendrecv, MPI_DOUBLE_PRECISION, idny, tag, MPI_COMM_WORLD, req(3), ierr)
  tag = 2
  call MPI_IRECV(real_neg_recv, size_sendrecv, MPI_DOUBLE_PRECISION, idny, tag, MPI_COMM_WORLD, req(2), ierr)
  call MPI_ISEND(real_pos_send, size_sendrecv, MPI_DOUBLE_PRECISION, idpy, tag, MPI_COMM_WORLD, req(4), ierr)

! Wait to receive
  call MPI_WAITALL(2, req(1:2), MPI_STATUSES_IGNORE, ierr)

  
  if (yflow) then

    if (myid - idny .gt. 0) then
      E(:,:,1,:)          = real_neg_recv(:,:,:)
    end if

    if (idpy - myid .gt. 0) then
      E(:,:,Ny_local+2,:) = real_pos_recv(:,:,:)
    end if
    
  else

  ! periodic
    E(:,:,1,:)            = real_neg_recv(:,:,:)
    E(:,:,Ny_local+2,:)   = real_pos_recv(:,:,:)
    
  end if
 
  call MPI_WAITALL(2, req(3:4), MPI_STATUSES_IGNORE, ierr)

  deallocate( real_neg_send )
  deallocate( real_pos_send )
  deallocate( real_neg_recv )
  deallocate( real_pos_recv )

  
! X DIRECTION

  allocate( real_neg_send(3,  Ny_local+2,  Nz_local+2) )
  allocate( real_pos_send(3,  Ny_local+2,  Nz_local+2) )
  allocate( real_neg_recv(3,  Ny_local+2,  Nz_local+2) )
  allocate( real_pos_recv(3,  Ny_local+2,  Nz_local+2) )
  size_sendrecv        = (3)*(Ny_local+2)*(Nz_local+2)
  real_neg_send(:, :, :) = E(:, 2, :, :)
  real_pos_send(:, :, :) = E(:, Nx_local+1, :, :)

  
  tag = 1
  call MPI_IRECV(real_pos_recv, size_sendrecv, MPI_DOUBLE_PRECISION, idpx, tag, MPI_COMM_WORLD, req(1), ierr)
  call MPI_ISEND(real_neg_send, size_sendrecv, MPI_DOUBLE_PRECISION, idnx, tag, MPI_COMM_WORLD, req(3), ierr)
  tag = 2
  call MPI_IRECV(real_neg_recv, size_sendrecv, MPI_DOUBLE_PRECISION, idnx, tag, MPI_COMM_WORLD, req(2), ierr)
  call MPI_ISEND(real_pos_send, size_sendrecv, MPI_DOUBLE_PRECISION, idpx, tag, MPI_COMM_WORLD, req(4), ierr)

! Wait to receive
  call MPI_WAITALL(2, req(1:2), MPI_STATUSES_IGNORE, ierr)

  
  if (xflow) then

    if (myid - idnx .gt. 0) then
      E(:,1,:,:)          = real_neg_recv(:,:,:)
    end if

    if (idpx - myid .gt. 0) then
      E(:,Nx_local+2,:,:) = real_pos_recv(:,:,:)
    end if
    
  else

  ! periodic
    E(:,1,:,:)            = real_neg_recv(:,:,:)
    E(:,Nx_local+2,:,:)   = real_pos_recv(:,:,:)
    
  end if
 
  call MPI_WAITALL(2, req(3:4), MPI_STATUSES_IGNORE, ierr)


  return
end subroutine CalcEfield
