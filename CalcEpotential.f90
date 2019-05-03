! Conjugate gradient method for solving Poisson's equation by Yousef Saad,
!Iterative Methods for Sparse Linear Systems, p.200 a.6.18.

subroutine CalcEpotential(U,C,Nx_local,Ny_local,Nz_local, &
     Lx_min,Lx_max,Ly_min,Ly_max,Lz_min,Lz_max,xflow,yflow,zflow,&
     xmin,xmax,ymin,ymax,zmin,zmax,dxyz,maxerr,maxit, &
     iprocs,jprocs,kprocs,Nspecies,myid)
  use mpi
  implicit none


  logical, intent(in) :: xflow, yflow, zflow
  integer, intent(in) :: Nx_local,Ny_local,Nz_local,maxit, &
       iprocs,jprocs,kprocs,Nspecies,myid
  real*8, intent(in) :: Lx_min,Lx_max,Ly_min,Ly_max,Lz_min,Lz_max, &
       xmin,xmax,ymin,ymax,zmin,zmax,maxerr,dxyz(3)
  real*8, intent(in) :: C(Nx_local+2,Ny_local+2,Nz_local+2)
  real*8, intent(inout) :: U(Nx_local+2,Ny_local+2,Nz_local+2)

  real*8, parameter :: eps0 = 8.854187817d-12  ! Electric constant [ SI ]

! local variables
  integer :: ss, idpx, idnx, idpy, idny, idpz, idnz, tag, ierr, req(8), &
       size_sendrecv
  real*8, allocatable :: real_neg_send(:,:), real_neg_recv(:,:), &
       real_pos_send(:,:), real_pos_recv(:,:)
  real*8 :: dV, b_local, b_global, b_global_sqrt, r_local, r_global, &
       r_old, maxerr_sqrt

! For the algorithm
  real*8 :: b(Nx_local+2,Ny_local+2,Nz_local+2)
  real*8 :: res(Nx_local+2,Ny_local+2,Nz_local+2)
  real*8 :: ax, ay, az, af

  b = 0.0d0
  res = 0.0d0
  r_local = 0.0d0
  dV = dxyz(1)*dxyz(2)*dxyz(3)
  maxerr_sqrt = dsqrt(maxerr)

  idpx = myid + (mod(myid+iprocs+1,                           iprocs) - mod(myid+iprocs, iprocs))
  idnx = myid + (mod(myid+iprocs-1,                           iprocs) - mod(myid+iprocs, iprocs))
  idpy = myid + (mod(myid+iprocs*jprocs+iprocs,               iprocs*jprocs) - mod(myid+iprocs*jprocs, iprocs*jprocs))
  idny = myid + (mod(myid+iprocs*jprocs-iprocs,               iprocs*jprocs) - mod(myid+iprocs*jprocs, iprocs*jprocs))
  idpz = myid + (mod(myid+iprocs*jprocs*kprocs+iprocs*jprocs, iprocs*jprocs*kprocs) &
               - mod(myid+iprocs*jprocs*kprocs,               iprocs*jprocs*kprocs))
  idnz = myid + (mod(myid+iprocs*jprocs*kprocs-iprocs*jprocs, iprocs*jprocs*kprocs) &
               - mod(myid+iprocs*jprocs*kprocs,               iprocs*jprocs*kprocs))


  ax = 0.5d0*( (dxyz(1)*dxyz(2))**2.0d0 + (dxyz(1)*dxyz(3))**2.0d0 + (dxyz(2)*dxyz(3))**2.0d0 )**(-1.0d0) * (dxyz(2)*dxyz(3))**2.0d0
  ay = 0.5d0*( (dxyz(1)*dxyz(2))**2.0d0 + (dxyz(1)*dxyz(3))**2.0d0 + (dxyz(2)*dxyz(3))**2.0d0 )**(-1.0d0) * (dxyz(1)*dxyz(3))**2.0d0
  az = 0.5d0*( (dxyz(1)*dxyz(2))**2.0d0 + (dxyz(1)*dxyz(3))**2.0d0 + (dxyz(2)*dxyz(3))**2.0d0 )**(-1.0d0) * (dxyz(1)*dxyz(2))**2.0d0
  af = 0.5d0*( (dxyz(1)*dxyz(2))**2.0d0 + (dxyz(1)*dxyz(3))**2.0d0 + (dxyz(2)*dxyz(3))**2.0d0 )**(-1.0d0) * (dxyz(1)*dxyz(2)*dxyz(3))**2.0d0

  
  b(:,:,:) = -C(:,:,:)*(dV*eps0)**(-1.0d0)

  b_local = sum( b(2:Nx_local+1,2:Ny_local+1,2:Nz_local+1)**2.0d0 )
  call MPI_allreduce( b_local, b_global, 1, MPI_DOUBLE_PRECISION, &
       MPI_SUM, MPI_comm_world, ierr )
  if (ierr>0) then
     write (*,*) 'CalcEpotential, 1, ierr=',ierr,' myid=',myid
  end if
  b_global_sqrt = dsqrt(b_global)
  r_global = 1.0d0*b_global
  r_old = 0.0d0

  do ss = 1, 2*maxit+1
    
    if (ss .le. maxit) then
      if (dsqrt(r_global) .le. maxerr*b_global_sqrt ) then
        exit
!!!! REMOVE THE FOLLOWING COMMENTS IF YOU WANT TO ALLOW LOCAL ERRORS
!      elif (dsqrt(abs(r_global-r_old)) .le. maxerr*b_global_sqrt ) then
!        exit
      end if
    else
      write (*,*) 'Maximum iteration reached!'
      ! call MPI_abort( MPI_comm_world, ierr )
    end if

    r_old = 1.0d0*r_global

    U(2:Nx_local+1,2:Ny_local+1,2:Nz_local+1) = ax*( U(3:Nx_local+2,2:Ny_local+1,2:Nz_local+1) + &
                                                     U(1:Nx_local+0,2:Ny_local+1,2:Nz_local+1) ) &
                                              + ay*( U(2:Nx_local+1,3:Ny_local+2,2:Nz_local+1) + &
                                                     U(2:Nx_local+1,1:Ny_local+0,2:Nz_local+1) ) &
                                              + az*( U(2:Nx_local+1,2:Ny_local+1,3:Nz_local+2) + &
                                                     U(2:Nx_local+1,2:Ny_local+1,1:Nz_local+0) ) &
                                              - af*( b(2:Nx_local+1,2:Ny_local+1,2:Nz_local+1) )

 
    res(2:Nx_local+1,2:Ny_local+1,2:Nz_local+1) = b(2:Nx_local+1,2:Ny_local+1,2:Nz_local+1) - &
                                              ( ( U(3:Nx_local+2,2:Ny_local+1,2:Nz_local+1) &
                                                + U(1:Nx_local+0,2:Ny_local+1,2:Nz_local+1) &
                                          - 2.0d0*U(2:Nx_local+1,2:Ny_local+1,2:Nz_local+1) )/dxyz(1)**2.0d0 &
                                              + ( U(2:Nx_local+1,3:Ny_local+2,2:Nz_local+1) &
                                                + U(2:Nx_local+1,1:Ny_local+0,2:Nz_local+1) &
                                          - 2.0d0*U(2:Nx_local+1,2:Ny_local+1,2:Nz_local+1) )/dxyz(2)**2.0d0 &
                                              + ( U(2:Nx_local+1,2:Ny_local+1,3:Nz_local+2) &
                                                + U(2:Nx_local+1,2:Ny_local+1,1:Nz_local+0) &
                                          - 2.0d0*U(2:Nx_local+1,2:Ny_local+1,2:Nz_local+1) )/dxyz(3)**2.0d0 )

! Residual
    r_local = sum( res(2:Nx_local+1,2:Ny_local+1,2:Nz_local+1)**2.0d0 )
    call MPI_allreduce( r_local, r_global, 1, MPI_DOUBLE_PRECISION, &
         MPI_SUM, MPI_comm_world, ierr )
    if (ierr>0) then
       write (*,*) 'CalcEpotential, 2, ierr=',ierr,' myid=',myid
    end if

! Send potential around
    
! Z DIRECTION

    allocate( real_neg_send(Nx_local+2,  Ny_local+2) )
    allocate( real_pos_send(Nx_local+2,  Ny_local+2) )
    allocate( real_neg_recv(Nx_local+2,  Ny_local+2) )
    allocate( real_pos_recv(Nx_local+2,  Ny_local+2) )
    size_sendrecv        = (Nx_local+2)*(Ny_local+2)
    real_neg_send(:, :)  = U(:, :, 2)
    real_pos_send(:, :)  = U(:, :, Nz_local+1)

  
    tag = 1
    call MPI_IRECV(real_pos_recv, size_sendrecv, MPI_DOUBLE_PRECISION, &
         idpz, tag, MPI_COMM_WORLD, req(1), ierr)
    if (ierr>0) then
       write (*,*) 'CalcEpotential, 3, ierr=',ierr,' myid=',myid
    end if
    call MPI_ISEND(real_neg_send, size_sendrecv, MPI_DOUBLE_PRECISION, &
         idnz, tag, MPI_COMM_WORLD, req(3), ierr)
    tag = 2
    if (ierr>0) then
       write (*,*) 'CalcEpotential, 4, ierr=',ierr,' myid=',myid
    end if
    call MPI_IRECV(real_neg_recv, size_sendrecv, MPI_DOUBLE_PRECISION, &
         idnz, tag, MPI_COMM_WORLD, req(2), ierr)
    if (ierr>0) then
       write (*,*) 'CalcEpotential, 5, ierr=',ierr,' myid=',myid
    end if
    call MPI_ISEND(real_pos_send, size_sendrecv, MPI_DOUBLE_PRECISION, &
         idpz, tag, MPI_COMM_WORLD, req(4), ierr)
    if (ierr>0) then
       write (*,*) 'CalcEpotential, 6, ierr=',ierr,' myid=',myid
    end if

! Wait to receive
    call MPI_WAITALL(2, req(1:2), MPI_STATUSES_IGNORE, ierr)
    if (ierr>0) then
       write (*,*) 'CalcEpotential, 7, ierr=',ierr,' myid=',myid
    end if

  
    if (zflow) then

      if (myid - idnz .gt. 0) then
        U(:,:,1)          = real_neg_recv(:,:)
      end if

      if (idpz - myid .gt. 0) then
        U(:,:,Nz_local+2) = real_pos_recv(:,:)
      end if
    
    else

    ! periodic
      U(:,:,1)            = real_neg_recv(:,:)
      U(:,:,Nz_local+2)   = real_pos_recv(:,:)
    
    end if
 
    call MPI_WAITALL(2, req(3:4), MPI_STATUSES_IGNORE, ierr)
    if (ierr>0) then
       write (*,*) 'CalcEpotential, 8, ierr=',ierr,' myid=',myid
    end if

    deallocate( real_neg_send )
    deallocate( real_pos_send )
    deallocate( real_neg_recv )
    deallocate( real_pos_recv )

! Y DIRECTION

    allocate( real_neg_send(Nx_local+2,  Nz_local+2) )
    allocate( real_pos_send(Nx_local+2,  Nz_local+2) )
    allocate( real_neg_recv(Nx_local+2,  Nz_local+2) )
    allocate( real_pos_recv(Nx_local+2,  Nz_local+2) )
    size_sendrecv        = (Nx_local+2)*(Nz_local+2)
    real_neg_send(:, :)  = U(:, 2, :)
    real_pos_send(:, :)  = U(:, Ny_local+1, :)

  
    tag = 1
    call MPI_IRECV(real_pos_recv, size_sendrecv, MPI_DOUBLE_PRECISION, &
         idpy, tag, MPI_COMM_WORLD, req(1), ierr)
    if (ierr>0) then
       write (*,*) 'CalcEpotential, 9, ierr=',ierr,' myid=',myid
    end if
    call MPI_ISEND(real_neg_send, size_sendrecv, MPI_DOUBLE_PRECISION, &
         idny, tag, MPI_COMM_WORLD, req(3), ierr)
    if (ierr>0) then
       write (*,*) 'CalcEpotential, 10, ierr=',ierr,' myid=',myid
    end if
    tag = 2
    call MPI_IRECV(real_neg_recv, size_sendrecv, MPI_DOUBLE_PRECISION, &
         idny, tag, MPI_COMM_WORLD, req(2), ierr)
    if (ierr>0) then
       write (*,*) 'CalcEpotential, 11, ierr=',ierr,' myid=',myid
    end if
    call MPI_ISEND(real_pos_send, size_sendrecv, MPI_DOUBLE_PRECISION, &
         idpy, tag, MPI_COMM_WORLD, req(4), ierr)
    if (ierr>0) then
       write (*,*) 'CalcEpotential, 12, ierr=',ierr,' myid=',myid
    end if

! Wait to receive
    call MPI_WAITALL(2, req(1:2), MPI_STATUSES_IGNORE, ierr)
    if (ierr>0) then
       write (*,*) 'CalcEpotential, 13, ierr=',ierr,' myid=',myid
    end if

  
    if (yflow) then

      if (myid - idny .gt. 0) then
        U(:,1,:)          = real_neg_recv(:,:)
      end if

      if (idpy - myid .gt. 0) then
        U(:,Ny_local+2,:) = real_pos_recv(:,:)
      end if
    
    else

    ! periodic
      U(:,1,:)            = real_neg_recv(:,:)
      U(:,Ny_local+2,:)   = real_pos_recv(:,:)
    
    end if
 
    call MPI_WAITALL(2, req(3:4), MPI_STATUSES_IGNORE, ierr)
    if (ierr>0) then
       write (*,*) 'CalcEpotential, 14, ierr=',ierr,' myid=',myid
    end if

    deallocate( real_neg_send )
    deallocate( real_pos_send )
    deallocate( real_neg_recv )
    deallocate( real_pos_recv )

! X DIRECTION

    allocate( real_neg_send(Ny_local+2,  Nz_local+2) )
    allocate( real_pos_send(Ny_local+2,  Nz_local+2) )
    allocate( real_neg_recv(Ny_local+2,  Nz_local+2) )
    allocate( real_pos_recv(Ny_local+2,  Nz_local+2) )
    size_sendrecv        = (Ny_local+2)*(Nz_local+2)
    real_neg_send(:, :)  = U(2, :, :)
    real_pos_send(:, :)  = U(Nx_local+1, :, :)

  
    tag = 1
    call MPI_IRECV(real_pos_recv, size_sendrecv, MPI_DOUBLE_PRECISION, &
         idpx, tag, MPI_COMM_WORLD, req(1), ierr)
    if (ierr>0) then
       write (*,*) 'CalcEpotential, 15, ierr=',ierr,' myid=',myid
    end if
    call MPI_ISEND(real_neg_send, size_sendrecv, MPI_DOUBLE_PRECISION, &
         idnx, tag, MPI_COMM_WORLD, req(3), ierr)
    if (ierr>0) then
       write (*,*) 'CalcEpotential, 16, ierr=',ierr,' myid=',myid
    end if
    tag = 2
    call MPI_IRECV(real_neg_recv, size_sendrecv, MPI_DOUBLE_PRECISION, &
         idnx, tag, MPI_COMM_WORLD, req(2), ierr)
    if (ierr>0) then
       write (*,*) 'CalcEpotential, 17, ierr=',ierr,' myid=',myid
    end if
    call MPI_ISEND(real_pos_send, size_sendrecv, MPI_DOUBLE_PRECISION, &
         idpx, tag, MPI_COMM_WORLD, req(4), ierr)
    if (ierr>0) then
       write (*,*) 'CalcEpotential, 18, ierr=',ierr,' myid=',myid
    end if

! Wait to receive
    call MPI_WAITALL(2, req(1:2), MPI_STATUSES_IGNORE, ierr)
    if (ierr>0) then
       write (*,*) 'CalcEpotential, 19, ierr=',ierr,' myid=',myid
    end if

  
    if (xflow) then

      if (myid - idnx .gt. 0) then
        U(1,:,:)          = real_neg_recv(:,:)
      end if

      if (idpx - myid .gt. 0) then
        U(Nx_local+2,:,:) = real_pos_recv(:,:)
      end if
    
    else

    ! periodic
      U(1,:,:)            = real_neg_recv(:,:)
      U(Nx_local+2,:,:)   = real_pos_recv(:,:)
    
    end if
 
    call MPI_WAITALL(2, req(3:4), MPI_STATUSES_IGNORE, ierr)
    if (ierr>0) then
       write (*,*) 'CalcEpotential, 20, ierr=',ierr,' myid=',myid
    end if

    deallocate( real_neg_send )
    deallocate( real_pos_send )
    deallocate( real_neg_recv )
    deallocate( real_pos_recv )


  end do


  return
end subroutine CalcEpotential
