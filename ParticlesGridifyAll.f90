subroutine ParticlesGridifyAll(F,N,C,real_particles,int_particles,Lx_min,Lx_max,Ly_min,Ly_max,Lz_min,Lz_max,dxyz,&
                   xmin,xmax,ymin,ymax,zmin,zmax,Nspecies,Nx_local,Ny_local,Nz_local,max_per_proc,num_local,&
                   xflow,yflow,zflow,iprocs,jprocs,kprocs,charge,n2p,myid) 
  use mpi
  implicit none


! Parameters
  integer, intent(in) :: Nspecies, Nx_local, Ny_local, Nz_local, max_per_proc, iprocs, jprocs, kprocs, myid, num_local
  logical, intent(in) :: xflow, yflow, zflow
  real*8, intent(in) :: Lx_min,Lx_max,Ly_min,Ly_max,Lz_min,Lz_max,dxyz(3),xmin,xmax,ymin,ymax,zmin,zmax,charge(8),n2p(8)
  integer, intent(in) :: int_particles(max_per_proc)
  real*8, intent(in) :: real_particles(8,max_per_proc)
  real*8, intent(out) :: F( Nspecies, 3, Nx_local+2, Ny_local+2, Nz_local+2 )
  real*8, intent(out) :: N( Nspecies, Nx_local+2, Ny_local+2, Nz_local+2 )
  real*8, intent(out) :: C( Nx_local+2, Ny_local+2, Nz_local+2 )
!  real*8, intent(out) :: J( 3, Nx_local+2, Ny_local+2, Nz_local+2 )  ! Current density
      

! Local variables
  real*8 rest(3), drest(3), F1, F2, F3, F4, &
       F5, F6, F7, F8, dV
  integer :: ii, jj, kk, pp, index_local(3), size_sendrecv, size_sendrecv2
  real*8, allocatable :: real_neg_send(:,:,:,:), real_neg_recv(:,:,:,:)
  real*8, allocatable :: real_pos_send(:,:,:,:), real_pos_recv(:,:,:,:)
  real*8, allocatable :: real_neg_send2(:,:,:,:,:), real_neg_recv2(:,:,:,:,:)
  real*8, allocatable :: real_pos_send2(:,:,:,:,:), real_pos_recv2(:,:,:,:,:)
  integer :: idpx, idnx, idpy, idny, idpz, idnz

! For the implementation of non-blocking communication
  integer :: tag, ierr, req(8)

  
  idpx = myid + (mod(myid+iprocs+1,                           iprocs) - mod(myid+iprocs, iprocs))
  idnx = myid + (mod(myid+iprocs-1,                           iprocs) - mod(myid+iprocs, iprocs))
  idpy = myid + (mod(myid+iprocs*jprocs+iprocs,               iprocs*jprocs) - mod(myid+iprocs*jprocs, iprocs*jprocs))
  idny = myid + (mod(myid+iprocs*jprocs-iprocs,               iprocs*jprocs) - mod(myid+iprocs*jprocs, iprocs*jprocs))
  idpz = myid + (mod(myid+iprocs*jprocs*kprocs+iprocs*jprocs, iprocs*jprocs*kprocs) &
               - mod(myid+iprocs*jprocs*kprocs,               iprocs*jprocs*kprocs))
  idnz = myid + (mod(myid+iprocs*jprocs*kprocs-iprocs*jprocs, iprocs*jprocs*kprocs) &
               - mod(myid+iprocs*jprocs*kprocs,               iprocs*jprocs*kprocs))


!
! This is the 1d interpretation of 'rest' and 'drest'
!
!        cell i            cell i+1
!  |                 |                 |
!  |--------|-----o--|--------|--------|
!  |                 |                 |   
!
!           <-----|----------->
!             rest   drest
!
!           <----------------->
!                 dx
!

! Set a few local variables
  dV = dxyz(1)*dxyz(2)*dxyz(3)
      
! First we reinitiate the matrices to zero
  F = 0.0d0
  N = 0.0d0
  C = 0.0d0
      

! Loop over all species and particles
      
  do pp = num_local, 1, -1 
            

! Find the index corresponding to the lowest x, y, z, which the particle contributes charge to.

    rest(1) = dmod( real_particles(1,pp)-Lx_min+0.5d0*dxyz(1), dxyz(1) )
    rest(2) = dmod( real_particles(2,pp)-Ly_min+0.5d0*dxyz(2), dxyz(2) )
    rest(3) = dmod( real_particles(3,pp)-Lz_min+0.5d0*dxyz(3), dxyz(3) )
    drest   = dxyz-rest
    
    index_local(1) = 1+idint( (real_particles(1,pp)-Lx_min)/dxyz(1)+0.5d0 )
    index_local(2) = 1+idint( (real_particles(2,pp)-Ly_min)/dxyz(2)+0.5d0 )
    index_local(3) = 1+idint( (real_particles(3,pp)-Lz_min)/dxyz(3)+0.5d0 )
    

! Start interpolation         
    F1 = drest(1)*drest(2)*drest(3)/dV
    F2 =  rest(1)*drest(2)*drest(3)/dV
    F3 = drest(1)* rest(2)*drest(3)/dV
    F4 = drest(1)*drest(2)* rest(3)/dV
    F5 =  rest(1)* rest(2)*drest(3)/dV
    F6 = drest(1)* rest(2)* rest(3)/dV
    F7 =  rest(1)*drest(2)* rest(3)/dV
    F8 =  rest(1)* rest(2)* rest(3)/dV

 !   if (int_particles(pp) .lt. 1 .or. index_local(1) .lt. 1 .or. index_local(2) .lt. 1 .or. index_local(3) .lt. 1) then
 !     write(*,*) pp, num_local
 !     write(*,*) int_particles(pp)
 !     write(*,*) real_particles(:,pp)
 !     write(*,*) index_local(:)
 !     write(*,*) rest(:)
 !   end if

    N( int_particles(pp), index_local(1)+0, index_local(2)+0, index_local(3)+0 ) = &
    N( int_particles(pp), index_local(1)+0, index_local(2)+0, index_local(3)+0 ) + &
      F1*n2p( int_particles(pp) )
    N( int_particles(pp), index_local(1)+1, index_local(2)+0, index_local(3)+0 ) = &
    N( int_particles(pp), index_local(1)+1, index_local(2)+0, index_local(3)+0 ) + &
      F2*n2p( int_particles(pp) )
    N( int_particles(pp), index_local(1)+0, index_local(2)+1, index_local(3)+0 ) = &
    N( int_particles(pp), index_local(1)+0, index_local(2)+1, index_local(3)+0 ) + &
      F3*n2p( int_particles(pp) )
    N( int_particles(pp), index_local(1)+0, index_local(2)+0, index_local(3)+1 ) = &
    N( int_particles(pp), index_local(1)+0, index_local(2)+0, index_local(3)+1 ) + &
      F4*n2p( int_particles(pp) )
    N( int_particles(pp), index_local(1)+1, index_local(2)+1, index_local(3)+0 ) = &
    N( int_particles(pp), index_local(1)+1, index_local(2)+1, index_local(3)+0 ) + &
      F5*n2p( int_particles(pp) )
    N( int_particles(pp), index_local(1)+0, index_local(2)+1, index_local(3)+1 ) = &
    N( int_particles(pp), index_local(1)+0, index_local(2)+1, index_local(3)+1 ) + &
      F6*n2p( int_particles(pp) )
    N( int_particles(pp), index_local(1)+1, index_local(2)+0, index_local(3)+1 ) = &
    N( int_particles(pp), index_local(1)+1, index_local(2)+0, index_local(3)+1 ) + &
      F7*n2p( int_particles(pp) )
    N( int_particles(pp), index_local(1)+1, index_local(2)+1, index_local(3)+1 ) = &
    N( int_particles(pp), index_local(1)+1, index_local(2)+1, index_local(3)+1 ) + &
      F8*n2p( int_particles(pp) )

    F( int_particles(pp), 1, index_local(1)+0, index_local(2)+0, index_local(3)+0 ) = &
    F( int_particles(pp), 1, index_local(1)+0, index_local(2)+0, index_local(3)+0 ) + &
      F1*n2p( int_particles(pp) )*real_particles(4,pp)
    F( int_particles(pp), 1, index_local(1)+1, index_local(2)+0, index_local(3)+0 ) = &
    F( int_particles(pp), 1, index_local(1)+1, index_local(2)+0, index_local(3)+0 ) + &
      F2*n2p( int_particles(pp) )*real_particles(4,pp)
    F( int_particles(pp), 1, index_local(1)+0, index_local(2)+1, index_local(3)+0 ) = &
    F( int_particles(pp), 1, index_local(1)+0, index_local(2)+1, index_local(3)+0 ) + &
      F3*n2p( int_particles(pp) )*real_particles(4,pp)
    F( int_particles(pp), 1, index_local(1)+0, index_local(2)+0, index_local(3)+1 ) = &
    F( int_particles(pp), 1, index_local(1)+0, index_local(2)+0, index_local(3)+1 ) + &
      F4*n2p( int_particles(pp) )*real_particles(4,pp)
    F( int_particles(pp), 1, index_local(1)+1, index_local(2)+1, index_local(3)+0 ) = &
    F( int_particles(pp), 1, index_local(1)+1, index_local(2)+1, index_local(3)+0 ) + &
      F5*n2p( int_particles(pp) )*real_particles(4,pp)
    F( int_particles(pp), 1, index_local(1)+0, index_local(2)+1, index_local(3)+1 ) = &
    F( int_particles(pp), 1, index_local(1)+0, index_local(2)+1, index_local(3)+1 ) + &
      F6*n2p( int_particles(pp) )*real_particles(4,pp)
    F( int_particles(pp), 1, index_local(1)+1, index_local(2)+0, index_local(3)+1 ) = &
    F( int_particles(pp), 1, index_local(1)+1, index_local(2)+0, index_local(3)+1 ) + &
      F7*n2p( int_particles(pp) )*real_particles(4,pp)
    F( int_particles(pp), 1, index_local(1)+1, index_local(2)+1, index_local(3)+1 ) = &
    F( int_particles(pp), 1, index_local(1)+1, index_local(2)+1, index_local(3)+1 ) + &
      F8*n2p( int_particles(pp) )*real_particles(4,pp)

    F( int_particles(pp), 2, index_local(1)+0, index_local(2)+0, index_local(3)+0 ) = &
    F( int_particles(pp), 2, index_local(1)+0, index_local(2)+0, index_local(3)+0 ) + &
      F1*n2p( int_particles(pp) )*real_particles(5,pp)
    F( int_particles(pp), 2, index_local(1)+1, index_local(2)+0, index_local(3)+0 ) = &
    F( int_particles(pp), 2, index_local(1)+1, index_local(2)+0, index_local(3)+0 ) + &
      F2*n2p( int_particles(pp) )*real_particles(5,pp)
    F( int_particles(pp), 2, index_local(1)+0, index_local(2)+1, index_local(3)+0 ) = &
    F( int_particles(pp), 2, index_local(1)+0, index_local(2)+1, index_local(3)+0 ) + &
      F3*n2p( int_particles(pp) )*real_particles(5,pp)
    F( int_particles(pp), 2, index_local(1)+0, index_local(2)+0, index_local(3)+1 ) = &
    F( int_particles(pp), 2, index_local(1)+0, index_local(2)+0, index_local(3)+1 ) + &
      F4*n2p( int_particles(pp) )*real_particles(5,pp)
    F( int_particles(pp), 2, index_local(1)+1, index_local(2)+1, index_local(3)+0 ) = &
    F( int_particles(pp), 2, index_local(1)+1, index_local(2)+1, index_local(3)+0 ) + &
      F5*n2p( int_particles(pp) )*real_particles(5,pp)
    F( int_particles(pp), 2, index_local(1)+0, index_local(2)+1, index_local(3)+1 ) = &
    F( int_particles(pp), 2, index_local(1)+0, index_local(2)+1, index_local(3)+1 ) + &
      F6*n2p( int_particles(pp) )*real_particles(5,pp)
    F( int_particles(pp), 2, index_local(1)+1, index_local(2)+0, index_local(3)+1 ) = &
    F( int_particles(pp), 2, index_local(1)+1, index_local(2)+0, index_local(3)+1 ) + &
      F7*n2p( int_particles(pp) )*real_particles(5,pp)
    F( int_particles(pp), 2, index_local(1)+1, index_local(2)+1, index_local(3)+1 ) = &
    F( int_particles(pp), 2, index_local(1)+1, index_local(2)+1, index_local(3)+1 ) + &
      F8*n2p( int_particles(pp) )*real_particles(5,pp)

    F( int_particles(pp), 3, index_local(1)+0, index_local(2)+0, index_local(3)+0 ) = &
    F( int_particles(pp), 3, index_local(1)+0, index_local(2)+0, index_local(3)+0 ) + &
      F1*n2p( int_particles(pp) )*real_particles(6,pp)
    F( int_particles(pp), 3, index_local(1)+1, index_local(2)+0, index_local(3)+0 ) = &
    F( int_particles(pp), 3, index_local(1)+1, index_local(2)+0, index_local(3)+0 ) + &
      F2*n2p( int_particles(pp) )*real_particles(6,pp)
    F( int_particles(pp), 3, index_local(1)+0, index_local(2)+1, index_local(3)+0 ) = &
    F( int_particles(pp), 3, index_local(1)+0, index_local(2)+1, index_local(3)+0 ) + &
      F3*n2p( int_particles(pp) )*real_particles(6,pp)
    F( int_particles(pp), 3, index_local(1)+0, index_local(2)+0, index_local(3)+1 ) = &
    F( int_particles(pp), 3, index_local(1)+0, index_local(2)+0, index_local(3)+1 ) + &
      F4*n2p( int_particles(pp) )*real_particles(6,pp)
    F( int_particles(pp), 3, index_local(1)+1, index_local(2)+1, index_local(3)+0 ) = &
    F( int_particles(pp), 3, index_local(1)+1, index_local(2)+1, index_local(3)+0 ) + &
      F5*n2p( int_particles(pp) )*real_particles(6,pp)
    F( int_particles(pp), 3, index_local(1)+0, index_local(2)+1, index_local(3)+1 ) = &
    F( int_particles(pp), 3, index_local(1)+0, index_local(2)+1, index_local(3)+1 ) + &
      F6*n2p( int_particles(pp) )*real_particles(6,pp)
    F( int_particles(pp), 3, index_local(1)+1, index_local(2)+0, index_local(3)+1 ) = &
    F( int_particles(pp), 3, index_local(1)+1, index_local(2)+0, index_local(3)+1 ) + &
      F7*n2p( int_particles(pp) )*real_particles(6,pp)
    F( int_particles(pp), 3, index_local(1)+1, index_local(2)+1, index_local(3)+1 ) = &
    F( int_particles(pp), 3, index_local(1)+1, index_local(2)+1, index_local(3)+1 ) + &
      F8*n2p( int_particles(pp) )*real_particles(6,pp)

  end do

!!!! MPI part for N, F

! Z BOUNDARY

  allocate( real_neg_send(Nspecies,  Nx_local+2,  Ny_local+2,  2) )
  allocate( real_pos_send(Nspecies,  Nx_local+2,  Ny_local+2,  2) )
  allocate( real_neg_recv(Nspecies,  Nx_local+2,  Ny_local+2,  2) )
  allocate( real_pos_recv(Nspecies,  Nx_local+2,  Ny_local+2,  2) )
  size_sendrecv        = (Nspecies)*(Nx_local+2)*(Ny_local+2)*(2)

  allocate( real_neg_send2(Nspecies,  3,  Nx_local+2,  Ny_local+2,  2) )
  allocate( real_pos_send2(Nspecies,  3,  Nx_local+2,  Ny_local+2,  2) )
  allocate( real_neg_recv2(Nspecies,  3,  Nx_local+2,  Ny_local+2,  2) )
  allocate( real_pos_recv2(Nspecies,  3,  Nx_local+2,  Ny_local+2,  2) )
  size_sendrecv2        = (Nspecies)*(3)*(Nx_local+2)*(Ny_local+2)*(2)

  real_neg_send(:,:,:,:)    = N(:,:,:,1:2)
  real_pos_send(:,:,:,:)    = N(:,:,:,Nz_local+1:Nz_local+2)
  real_neg_send2(:,:,:,:,:) = F(:,:,:,:,1:2)
  real_pos_send2(:,:,:,:,:) = F(:,:,:,:,Nz_local+1:Nz_local+2)

  tag = 1
  call MPI_IRECV(real_pos_recv, size_sendrecv, MPI_DOUBLE_PRECISION, idpz, tag, MPI_COMM_WORLD, req(1), ierr)
  call MPI_ISEND(real_neg_send, size_sendrecv, MPI_DOUBLE_PRECISION, idnz, tag, MPI_COMM_WORLD, req(5), ierr)
  tag = 2
  call MPI_IRECV(real_neg_recv, size_sendrecv, MPI_DOUBLE_PRECISION, idnz, tag, MPI_COMM_WORLD, req(2), ierr)
  call MPI_ISEND(real_pos_send, size_sendrecv, MPI_DOUBLE_PRECISION, idpz, tag, MPI_COMM_WORLD, req(6), ierr)
  tag = 3
  call MPI_IRECV(real_pos_recv2, size_sendrecv2, MPI_DOUBLE_PRECISION, idpz, tag, MPI_COMM_WORLD, req(3), ierr)
  call MPI_ISEND(real_neg_send2, size_sendrecv2, MPI_DOUBLE_PRECISION, idnz, tag, MPI_COMM_WORLD, req(7), ierr)
  tag = 4
  call MPI_IRECV(real_neg_recv2, size_sendrecv2, MPI_DOUBLE_PRECISION, idnz, tag, MPI_COMM_WORLD, req(4), ierr)
  call MPI_ISEND(real_pos_send2, size_sendrecv2, MPI_DOUBLE_PRECISION, idpz, tag, MPI_COMM_WORLD, req(8), ierr)

! Wait to receive
  call MPI_WAITALL(4, req(1:4), MPI_STATUSES_IGNORE, ierr)


  if (zflow) then

    if (myid - idnz .gt. 0) then
      N(:,:,:,1:2)                   = N(:,:,:,1:2) + real_neg_recv(:,:,:,:)
      F(:,:,:,:,1:2)                 = F(:,:,:,:,1:2) + real_neg_recv2(:,:,:,:,:)
    else
      N(:,:,:,2)                     = N(:,:,:,2) + N(:,:,:,1)
      N(:,:,:,1)                     = 0.0d0
      F(:,:,:,:,2)                   = F(:,:,:,:,2) + F(:,:,:,:,1)
      F(:,:,:,:,1)                   = 0.0d0
    end if

    if (idpz - myid .gt. 0) then
      N(:,:,:,Nz_local+1:Nz_local+2)    = N(:,:,:,Nz_local+1:Nz_local+2) + real_pos_recv(:,:,:,:)
      F(:,:,:,:,Nz_local+1:Nz_local+2)  = F(:,:,:,:,Nz_local+1:Nz_local+2) + real_pos_recv2(:,:,:,:,:)
    else
      N(:,:,:,Nz_local+1)            = N(:,:,:,Nz_local+1) + N(:,:,:,Nz_local+2)
      N(:,:,:,Nz_local+2)            = 0.0d0
      F(:,:,:,:,Nz_local+1)          = F(:,:,:,:,Nz_local+1) + F(:,:,:,:,Nz_local+2)
      F(:,:,:,:,Nz_local+2)          = 0.0d0
    end if
    
  else

  ! periodic
    N(:,:,:,1:2)                     = N(:,:,:,1:2) + real_neg_recv(:,:,:,:)
    N(:,:,:,Nz_local+1:Nz_local+2)   = N(:,:,:,Nz_local+1:Nz_local+2) + real_pos_recv(:,:,:,:)
    F(:,:,:,:,1:2)                   = F(:,:,:,:,1:2) + real_neg_recv2(:,:,:,:,:)
    F(:,:,:,:,Nz_local+1:Nz_local+2) = F(:,:,:,:,Nz_local+1:Nz_local+2) + real_pos_recv2(:,:,:,:,:)
    
  end if
 
  call MPI_WAITALL(4, req(5:8), MPI_STATUSES_IGNORE, ierr)

  deallocate( real_neg_send )
  deallocate( real_pos_send )
  deallocate( real_neg_recv )
  deallocate( real_pos_recv )
  deallocate( real_neg_send2 )
  deallocate( real_pos_send2 )
  deallocate( real_neg_recv2 )
  deallocate( real_pos_recv2 )


! Y BOUNDARY

  allocate( real_neg_send(Nspecies,  Nx_local+2,  2,  Nz_local+2) )
  allocate( real_pos_send(Nspecies,  Nx_local+2,  2,  Nz_local+2) )
  allocate( real_neg_recv(Nspecies,  Nx_local+2,  2,  Nz_local+2) )
  allocate( real_pos_recv(Nspecies,  Nx_local+2,  2,  Nz_local+2) )
  size_sendrecv        = (Nspecies)*(Nx_local+2)*(2)*(Nz_local+2)

  allocate( real_neg_send2(Nspecies,  3,  Nx_local+2,  2,  Nz_local+2) )
  allocate( real_pos_send2(Nspecies,  3,  Nx_local+2,  2,  Nz_local+2) )
  allocate( real_neg_recv2(Nspecies,  3,  Nx_local+2,  2,  Nz_local+2) )
  allocate( real_pos_recv2(Nspecies,  3,  Nx_local+2,  2,  Nz_local+2) )
  size_sendrecv2        = (Nspecies)*(3)*(Nx_local+2)*(2)*(Nz_local+2)

  real_neg_send(:,:,:,:)    = N(:,:,1:2,:)
  real_pos_send(:,:,:,:)    = N(:,:,Ny_local+1:Ny_local+2,:)
  real_neg_send2(:,:,:,:,:) = F(:,:,:,1:2,:)
  real_pos_send2(:,:,:,:,:) = F(:,:,:,Ny_local+1:Ny_local+2,:)

  tag = 1
  call MPI_IRECV(real_pos_recv, size_sendrecv, MPI_DOUBLE_PRECISION, idpy, tag, MPI_COMM_WORLD, req(1), ierr)
  call MPI_ISEND(real_neg_send, size_sendrecv, MPI_DOUBLE_PRECISION, idny, tag, MPI_COMM_WORLD, req(5), ierr)
  tag = 2
  call MPI_IRECV(real_neg_recv, size_sendrecv, MPI_DOUBLE_PRECISION, idny, tag, MPI_COMM_WORLD, req(2), ierr)
  call MPI_ISEND(real_pos_send, size_sendrecv, MPI_DOUBLE_PRECISION, idpy, tag, MPI_COMM_WORLD, req(6), ierr)
  tag = 3
  call MPI_IRECV(real_pos_recv2, size_sendrecv2, MPI_DOUBLE_PRECISION, idpy, tag, MPI_COMM_WORLD, req(3), ierr)
  call MPI_ISEND(real_neg_send2, size_sendrecv2, MPI_DOUBLE_PRECISION, idny, tag, MPI_COMM_WORLD, req(7), ierr)
  tag = 4
  call MPI_IRECV(real_neg_recv2, size_sendrecv2, MPI_DOUBLE_PRECISION, idny, tag, MPI_COMM_WORLD, req(4), ierr)
  call MPI_ISEND(real_pos_send2, size_sendrecv2, MPI_DOUBLE_PRECISION, idpy, tag, MPI_COMM_WORLD, req(8), ierr)

! Wait to receive
  call MPI_WAITALL(4, req(1:4), MPI_STATUSES_IGNORE, ierr)


  if (yflow) then

    if (myid - idny .gt. 0) then
      N(:,:,1:2,:)                   = N(:,:,1:2,:) + real_neg_recv(:,:,:,:)
      F(:,:,:,1:2,:)                 = F(:,:,:,1:2,:) + real_neg_recv2(:,:,:,:,:)
    else
      N(:,:,2,:)                     = N(:,:,2,:) + N(:,:,1,:)
      N(:,:,1,:)                     = 0.0d0
      F(:,:,:,2,:)                   = F(:,:,:,2,:) + F(:,:,:,1,:)
      F(:,:,:,1,:)                   = 0.0d0
    end if

    if (idpy - myid .gt. 0) then
      N(:,:,Ny_local+1:Ny_local+2,:)   = N(:,:,Ny_local+1:Ny_local+2,:) + real_pos_recv(:,:,:,:)
      F(:,:,:,Ny_local+1:Ny_local+2,:) = F(:,:,:,Ny_local+1:Ny_local+2,:) + real_pos_recv2(:,:,:,:,:)
    else
      N(:,:,Ny_local+1,:)            = N(:,:,Ny_local+1,:) + N(:,:,Ny_local+2,:)
      N(:,:,Ny_local+2,:)            = 0.0d0
      F(:,:,:,Ny_local+1,:)          = F(:,:,:,Ny_local+1,:) + F(:,:,:,Ny_local+2,:)
      F(:,:,:,Ny_local+2,:)          = 0.0d0
    end if
    
  else

  ! periodic
    N(:,:,1:2,:)                     = N(:,:,1:2,:) + real_neg_recv(:,:,:,:)
    N(:,:,Ny_local+1:Ny_local+2,:)   = N(:,:,Ny_local+1:Ny_local+2,:) + real_pos_recv(:,:,:,:)
    F(:,:,:,1:2,:)                   = F(:,:,:,1:2,:) + real_neg_recv2(:,:,:,:,:)
    F(:,:,:,Ny_local+1:Ny_local+2,:) = F(:,:,:,Ny_local+1:Ny_local+2,:) + real_pos_recv2(:,:,:,:,:)
    
  end if
 
  call MPI_WAITALL(4, req(5:8), MPI_STATUSES_IGNORE, ierr)

  deallocate( real_neg_send )
  deallocate( real_pos_send )
  deallocate( real_neg_recv )
  deallocate( real_pos_recv )
  deallocate( real_neg_send2 )
  deallocate( real_pos_send2 )
  deallocate( real_neg_recv2 )
  deallocate( real_pos_recv2 )

! X BOUNDARY

  allocate( real_neg_send(Nspecies,  2,  Ny_local+2,  Nz_local+2) )
  allocate( real_pos_send(Nspecies,  2,  Ny_local+2,  Nz_local+2) )
  allocate( real_neg_recv(Nspecies,  2,  Ny_local+2,  Nz_local+2) )
  allocate( real_pos_recv(Nspecies,  2,  Ny_local+2,  Nz_local+2) )
  size_sendrecv        = (Nspecies)*(2)*(Ny_local+2)*(Nz_local+2)

  allocate( real_neg_send2(Nspecies,  3,  2,  Ny_local+2,  Nz_local+2) )
  allocate( real_pos_send2(Nspecies,  3,  2,  Ny_local+2,  Nz_local+2) )
  allocate( real_neg_recv2(Nspecies,  3,  2,  Ny_local+2,  Nz_local+2) )
  allocate( real_pos_recv2(Nspecies,  3,  2,  Ny_local+2,  Nz_local+2) )
  size_sendrecv2        = (Nspecies)*(3)*(2)*(Ny_local+2)*(Nz_local+2)

  real_neg_send(:,:,:,:)    = N(:,1:2,:,:)
  real_pos_send(:,:,:,:)    = N(:,Nx_local+1:Nx_local+2,:,:)
  real_neg_send2(:,:,:,:,:) = F(:,:,1:2,:,:)
  real_pos_send2(:,:,:,:,:) = F(:,:,Nx_local+1:Nx_local+2,:,:)

  tag = 1
  call MPI_IRECV(real_pos_recv, size_sendrecv, MPI_DOUBLE_PRECISION, idpx, tag, MPI_COMM_WORLD, req(1), ierr)
  call MPI_ISEND(real_neg_send, size_sendrecv, MPI_DOUBLE_PRECISION, idnx, tag, MPI_COMM_WORLD, req(5), ierr)
  tag = 2
  call MPI_IRECV(real_neg_recv, size_sendrecv, MPI_DOUBLE_PRECISION, idnx, tag, MPI_COMM_WORLD, req(2), ierr)
  call MPI_ISEND(real_pos_send, size_sendrecv, MPI_DOUBLE_PRECISION, idpx, tag, MPI_COMM_WORLD, req(6), ierr)
  tag = 3
  call MPI_IRECV(real_pos_recv2, size_sendrecv2, MPI_DOUBLE_PRECISION, idpx, tag, MPI_COMM_WORLD, req(3), ierr)
  call MPI_ISEND(real_neg_send2, size_sendrecv2, MPI_DOUBLE_PRECISION, idnx, tag, MPI_COMM_WORLD, req(7), ierr)
  tag = 4
  call MPI_IRECV(real_neg_recv2, size_sendrecv2, MPI_DOUBLE_PRECISION, idnx, tag, MPI_COMM_WORLD, req(4), ierr)
  call MPI_ISEND(real_pos_send2, size_sendrecv2, MPI_DOUBLE_PRECISION, idpx, tag, MPI_COMM_WORLD, req(8), ierr)

! Wait to receive
  call MPI_WAITALL(4, req(1:4), MPI_STATUSES_IGNORE, ierr)


  if (xflow) then

    if (myid - idnx .gt. 0) then
      N(:,1:2,:,:)                   = N(:,1:2,:,:) + real_neg_recv(:,:,:,:)
      F(:,:,1:2,:,:)                 = F(:,:,1:2,:,:) + real_neg_recv2(:,:,:,:,:)
    else
      N(:,2,:,:)                     = N(:,2,:,:) + N(:,1,:,:)
      N(:,1,:,:)                     = 0.0d0
      F(:,:,2,:,:)                   = F(:,:,2,:,:) + F(:,:,1,:,:)
      F(:,:,1,:,:)                   = 0.0d0
    end if

    if (idpx - myid .gt. 0) then
      N(:,Nx_local+1:Nx_local+2,:,:)   = N(:,Nx_local+1:Nx_local+2,:,:) + real_pos_recv(:,:,:,:)
      F(:,:,Nx_local+1:Nx_local+2,:,:) = F(:,:,Nx_local+1:Nx_local+2,:,:) + real_pos_recv2(:,:,:,:,:)
    else
      N(:,Nx_local+1,:,:)            = N(:,Nx_local+1,:,:) + N(:,Nx_local+2,:,:)
      N(:,Nx_local+2,:,:)            = 0.0d0
      F(:,:,Nx_local+1,:,:)          = F(:,:,Nx_local+1,:,:) + F(:,:,Nx_local+2,:,:)
      F(:,:,Nx_local+2,:,:)          = 0.0d0
    end if
    
  else

  ! periodic
    N(:,1:2,:,:)                     = N(:,1:2,:,:) + real_neg_recv(:,:,:,:)
    N(:,Nx_local+1:Nx_local+2,:,:)   = N(:,Nx_local+1:Nx_local+2,:,:) + real_pos_recv(:,:,:,:)
    F(:,:,1:2,:,:)                   = F(:,:,1:2,:,:) + real_neg_recv2(:,:,:,:,:)
    F(:,:,Nx_local+1:Nx_local+2,:,:) = F(:,:,Nx_local+1:Nx_local+2,:,:) + real_pos_recv2(:,:,:,:,:)
    
  end if
 
  call MPI_WAITALL(4, req(5:8), MPI_STATUSES_IGNORE, ierr)

  deallocate( real_neg_send )
  deallocate( real_pos_send )
  deallocate( real_neg_recv )
  deallocate( real_pos_recv )
  deallocate( real_neg_send2 )
  deallocate( real_pos_send2 )
  deallocate( real_neg_recv2 )
  deallocate( real_pos_recv2 )


  
  do pp = 1, Nspecies
    C(:,:,:) = C(:,:,:) + N(pp,:,:,:)*charge(pp)
!    J(1,:,:,:) = J(1,:,:,:) + F(pp,1,:,:,:)*charge(pp)  ! Current in x
!    J(2,:,:,:) = J(2,:,:,:) + F(pp,2,:,:,:)*charge(pp)  ! Current in y
!    J(3,:,:,:) = J(3,:,:,:) + F(pp,3,:,:,:)*charge(pp)  ! Current in z
  end do

  return
end subroutine ParticlesGridifyAll
