
subroutine ParticlesTransfer(real_particles,int_particles,Lx_min,Lx_max,Ly_min,Ly_max,Lz_min,Lz_max,&
                         xmin,xmax,ymin,ymax,zmin,zmax,max_per_proc,num_local,iprocs,jprocs,kprocs,myid)
  use mpi
  implicit none
  


! Parameters
  integer, intent(in) :: max_per_proc, iprocs, jprocs, kprocs, myid
  real*8, intent(in) :: Lx_min,Lx_max,Ly_min,Ly_max,Lz_min,Lz_max,xmin,xmax,ymin,ymax,zmin,zmax
  integer, intent(inout) :: int_particles(max_per_proc), num_local
  real*8, intent(inout) :: real_particles(8,max_per_proc)

  
! Local variables
  integer :: ii, jj, kk, pp, neg_index, pos_index
  real*8, allocatable :: real_neg_send(:,:), real_neg_recv(:,:), real_pos_send(:,:), real_pos_recv(:,:)
  integer, allocatable :: int_neg_send(:), int_neg_recv(:), int_pos_send(:), int_pos_recv(:)
  integer :: size_neg_send, size_neg_recv, size_pos_send, size_pos_recv 
  integer :: idpx, idnx, idpy, idny, idpz, idnz

! For the implementation of non-blocking communication
  integer :: tag, ierr, req(12)

  idpx = myid + (mod(myid+iprocs+1,                           iprocs) - mod(myid+iprocs, iprocs))
  idnx = myid + (mod(myid+iprocs-1,                           iprocs) - mod(myid+iprocs, iprocs))
  idpy = myid + (mod(myid+iprocs*jprocs+iprocs,               iprocs*jprocs) - mod(myid+iprocs*jprocs, iprocs*jprocs))
  idny = myid + (mod(myid+iprocs*jprocs-iprocs,               iprocs*jprocs) - mod(myid+iprocs*jprocs, iprocs*jprocs))
  idpz = myid + (mod(myid+iprocs*jprocs*kprocs+iprocs*jprocs, iprocs*jprocs*kprocs) &
               - mod(myid+iprocs*jprocs*kprocs,               iprocs*jprocs*kprocs))
  idnz = myid + (mod(myid+iprocs*jprocs*kprocs-iprocs*jprocs, iprocs*jprocs*kprocs) &
               - mod(myid+iprocs*jprocs*kprocs,               iprocs*jprocs*kprocs))


! Clean up before the z direction

  size_neg_send = 0
  size_neg_recv = 0
  size_pos_send = 0
  size_pos_recv = 0
  neg_index = 0
  pos_index = 0


!!!! Z DIRECTION

if (idpz .ne. myid) then

  do pp = num_local, 1, -1

! Switch memory places so that the ones we want to send end up in the end

    if (real_particles(3,pp) .lt. Lz_min) then

      if (pp .ne. num_local-size_neg_send-size_pos_send) then
        real_particles(:,pp) = real_particles(:,pp) + real_particles(:,num_local-size_neg_send-size_pos_send)
        real_particles(:,num_local-size_neg_send-size_pos_send) &
                             = real_particles(:,pp) - real_particles(:,num_local-size_neg_send-size_pos_send)
        real_particles(:,pp) = real_particles(:,pp) - real_particles(:,num_local-size_neg_send-size_pos_send)

        int_particles(pp)    = int_particles(pp)    + int_particles(num_local-size_neg_send-size_pos_send)
        int_particles(num_local-size_neg_send-size_pos_send) &
                             = int_particles(pp)    - int_particles(num_local-size_neg_send-size_pos_send)
        int_particles(pp)    = int_particles(pp)    - int_particles(num_local-size_neg_send-size_pos_send)
      end if

      size_neg_send = size_neg_send + 1

    else if (real_particles(3,pp) .gt. Lz_max) then

      if (pp .ne. num_local-size_neg_send-size_pos_send) then
        real_particles(:,pp) = real_particles(:,pp) + real_particles(:,num_local-size_neg_send-size_pos_send)
        real_particles(:,num_local-size_neg_send-size_pos_send) &
                             = real_particles(:,pp) - real_particles(:,num_local-size_neg_send-size_pos_send)
        real_particles(:,pp) = real_particles(:,pp) - real_particles(:,num_local-size_neg_send-size_pos_send)

        int_particles(pp)    = int_particles(pp)    + int_particles(num_local-size_neg_send-size_pos_send)
        int_particles(num_local-size_neg_send-size_pos_send) &
                             = int_particles(pp)    - int_particles(num_local-size_neg_send-size_pos_send)
        int_particles(pp)    = int_particles(pp)    - int_particles(num_local-size_neg_send-size_pos_send)
      end if

      size_pos_send = size_pos_send + 1    

    end if
  end do

! Send and receive array sizes

  tag = 1
  call MPI_IRECV(size_pos_recv, 1, MPI_INT, idpz, tag, MPI_COMM_WORLD, req(1), ierr)
  call MPI_ISEND(size_neg_send, 1, MPI_INT, idnz, tag, MPI_COMM_WORLD, req(7), ierr)
  tag = 2
  call MPI_IRECV(size_neg_recv, 1, MPI_INT, idnz, tag, MPI_COMM_WORLD, req(2), ierr)
  call MPI_ISEND(size_pos_send, 1, MPI_INT, idpz, tag, MPI_COMM_WORLD, req(8), ierr)

! Form the send matrices

  allocate( real_neg_send(8,size_neg_send) )
  allocate( int_neg_send(size_neg_send) )
  allocate( real_pos_send(8,size_pos_send) )
  allocate( int_pos_send(size_pos_send) )

if (size_pos_send+size_neg_send .gt. 0) then
  do pp = num_local, (num_local-size_neg_send-size_pos_send+1), -1

    if (real_particles(3,pp) .lt. Lz_min) then
      neg_index                    = neg_index+1
      int_neg_send(neg_index)      = int_particles(pp)
      real_neg_send(:,neg_index)   = real_particles(:,pp)
      if (real_neg_send(3,neg_index) .lt. zmin) then
        real_neg_send(3,neg_index) = real_neg_send(3,neg_index) + (zmax-zmin)
      end if

    else if (real_particles(3,pp) .gt. Lz_max) then
      pos_index                    = pos_index+1
      int_pos_send(pos_index)      = int_particles(pp)
      real_pos_send(:,pos_index)   = real_particles(:,pp)
      if (real_pos_send(3,pos_index) .gt. zmax) then
        real_pos_send(3,pos_index) = real_pos_send(3,pos_index) - (zmax-zmin)
      end if

    end if
  end do
end if

! Change local size

  num_local = num_local-size_neg_send-size_pos_send

! Wait to receive
  call MPI_WAITALL(2, req(1:2), MPI_STATUSES_IGNORE, ierr)


  allocate( real_neg_recv(8,size_neg_recv) )
  allocate( int_neg_recv(size_neg_recv) )
  allocate( real_pos_recv(8,size_pos_recv) )
  allocate( int_pos_recv(size_pos_recv) )

! Send and receive arrays

  tag = 3
  call MPI_IRECV(real_pos_recv, 8*size_pos_recv, MPI_DOUBLE_PRECISION, idpz, tag, MPI_COMM_WORLD, req(3), ierr)
  call MPI_ISEND(real_neg_send, 8*size_neg_send, MPI_DOUBLE_PRECISION, idnz, tag, MPI_COMM_WORLD, req(9), ierr)
  tag = 4
  call MPI_IRECV(real_neg_recv, 8*size_neg_recv, MPI_DOUBLE_PRECISION, idnz, tag, MPI_COMM_WORLD, req(4), ierr)
  call MPI_ISEND(real_pos_send, 8*size_pos_send, MPI_DOUBLE_PRECISION, idpz, tag, MPI_COMM_WORLD, req(10), ierr)

  tag = 5
  call MPI_IRECV(int_pos_recv, size_pos_recv, MPI_INT, idpz, tag, MPI_COMM_WORLD, req(5), ierr)
  call MPI_ISEND(int_neg_send, size_neg_send, MPI_INT, idnz, tag, MPI_COMM_WORLD, req(11), ierr)
  tag = 6
  call MPI_IRECV(int_neg_recv, size_neg_recv, MPI_INT, idnz, tag, MPI_COMM_WORLD, req(6), ierr)
  call MPI_ISEND(int_pos_send, size_pos_send, MPI_INT, idpz, tag, MPI_COMM_WORLD, req(12), ierr)

! Wait to receive
  call MPI_WAITALL(4, req(3:6), MPI_STATUSES_IGNORE, ierr)

! Organize received arrays

  if (size_pos_recv .gt. 0) then
    int_particles(num_local+1:num_local+size_pos_recv)    = int_pos_recv(:)
    real_particles(:,num_local+1:num_local+size_pos_recv) = real_pos_recv(:,:)
    num_local = num_local+size_pos_recv
  end if

  if (size_neg_recv .gt. 0) then
    int_particles(num_local+1:num_local+size_neg_recv)    = int_neg_recv(:)
    real_particles(:,num_local+1:num_local+size_neg_recv) = real_neg_recv(:,:)
    num_local = num_local+size_neg_recv
  end if

! Clean up before the y direction

  call MPI_WAITALL(6, req(7:12), MPI_STATUSES_IGNORE, ierr)


  deallocate( real_neg_recv )
  deallocate( real_pos_recv )
  deallocate( int_neg_recv )
  deallocate( int_pos_recv )
  deallocate( real_neg_send )
  deallocate( real_pos_send )
  deallocate( int_neg_send )
  deallocate( int_pos_send )
  size_neg_send = 0
  size_neg_recv = 0
  size_pos_send = 0
  size_pos_recv = 0
  neg_index = 0
  pos_index = 0

end if

!!!! Y DIRECTION

! If outside local y domain, prepare to send!

if (idpy .ne. myid) then

  do pp = num_local, 1, -1

! Switch memory places so that the ones we want to send end up in the end

    if (real_particles(2,pp) .lt. Ly_min) then

      if (pp .ne. num_local-size_neg_send-size_pos_send) then
        real_particles(:,pp) = real_particles(:,pp) + real_particles(:,num_local-size_neg_send-size_pos_send)
        real_particles(:,num_local-size_neg_send-size_pos_send) &
                             = real_particles(:,pp) - real_particles(:,num_local-size_neg_send-size_pos_send)
        real_particles(:,pp) = real_particles(:,pp) - real_particles(:,num_local-size_neg_send-size_pos_send)

        int_particles(pp)    = int_particles(pp)    + int_particles(num_local-size_neg_send-size_pos_send)
        int_particles(num_local-size_neg_send-size_pos_send) &
                             = int_particles(pp)    - int_particles(num_local-size_neg_send-size_pos_send)
        int_particles(pp)    = int_particles(pp)    - int_particles(num_local-size_neg_send-size_pos_send)
      end if

      size_neg_send = size_neg_send + 1

    else if (real_particles(2,pp) .gt. Ly_max) then

      if (pp .ne. num_local-size_neg_send-size_pos_send) then
        real_particles(:,pp) = real_particles(:,pp) + real_particles(:,num_local-size_neg_send-size_pos_send)
        real_particles(:,num_local-size_neg_send-size_pos_send) &
                             = real_particles(:,pp) - real_particles(:,num_local-size_neg_send-size_pos_send)
        real_particles(:,pp) = real_particles(:,pp) - real_particles(:,num_local-size_neg_send-size_pos_send)

        int_particles(pp)    = int_particles(pp)    + int_particles(num_local-size_neg_send-size_pos_send)
        int_particles(num_local-size_neg_send-size_pos_send) &
                             = int_particles(pp)    - int_particles(num_local-size_neg_send-size_pos_send)
        int_particles(pp)    = int_particles(pp)    - int_particles(num_local-size_neg_send-size_pos_send)
      end if

      size_pos_send = size_pos_send + 1    

    end if
  end do

! Send and receive array sizes

  tag = 1
  call MPI_IRECV(size_pos_recv, 1, MPI_INT, idpy, tag, MPI_COMM_WORLD, req(1), ierr)
  call MPI_ISEND(size_neg_send, 1, MPI_INT, idny, tag, MPI_COMM_WORLD, req(7), ierr)
  tag = 2
  call MPI_IRECV(size_neg_recv, 1, MPI_INT, idny, tag, MPI_COMM_WORLD, req(2), ierr)
  call MPI_ISEND(size_pos_send, 1, MPI_INT, idpy, tag, MPI_COMM_WORLD, req(8), ierr)

! Form the send matrices

  allocate( real_neg_send(8,size_neg_send) )
  allocate( int_neg_send(size_neg_send) )
  allocate( real_pos_send(8,size_pos_send) )
  allocate( int_pos_send(size_pos_send) )

if (size_pos_send+size_neg_send .gt. 0) then
  do pp = num_local, (num_local-size_neg_send-size_pos_send+1), -1

    if (real_particles(2,pp) .lt. Ly_min) then
      neg_index                    = neg_index+1
      int_neg_send(neg_index)      = int_particles(pp)
      real_neg_send(:,neg_index)   = real_particles(:,pp)
      if (real_neg_send(2,neg_index) .lt. ymin) then
        real_neg_send(2,neg_index) = real_neg_send(2,neg_index) + (ymax-ymin)
      end if

    else if (real_particles(2,pp) .gt. Ly_max) then
      pos_index                    = pos_index+1
      int_pos_send(pos_index)      = int_particles(pp)
      real_pos_send(:,pos_index)   = real_particles(:,pp)
      if (real_pos_send(2,pos_index) .gt. ymax) then
        real_pos_send(2,pos_index) = real_pos_send(2,pos_index) - (ymax-ymin)
      end if

    end if
  end do
end if

! Change local size

  num_local = num_local-size_neg_send-size_pos_send

! Wait to receive
  call MPI_WAITALL(2, req(1:2), MPI_STATUSES_IGNORE, ierr)

  allocate( real_neg_recv(8,size_neg_recv) )
  allocate( real_pos_recv(8,size_pos_recv) )
  allocate( int_neg_recv(size_neg_recv) )
  allocate( int_pos_recv(size_pos_recv) )

! Send and receive arrays

  tag = 3
  call MPI_IRECV(real_pos_recv, 8*size_pos_recv, MPI_DOUBLE_PRECISION, idpy, tag, MPI_COMM_WORLD, req(3), ierr)
  call MPI_ISEND(real_neg_send, 8*size_neg_send, MPI_DOUBLE_PRECISION, idny, tag, MPI_COMM_WORLD, req(9), ierr)
  tag = 4
  call MPI_IRECV(real_neg_recv, 8*size_neg_recv, MPI_DOUBLE_PRECISION, idny, tag, MPI_COMM_WORLD, req(4), ierr)
  call MPI_ISEND(real_pos_send, 8*size_pos_send, MPI_DOUBLE_PRECISION, idpy, tag, MPI_COMM_WORLD, req(10), ierr)

  tag = 5
  call MPI_IRECV(int_pos_recv, size_pos_recv, MPI_INT, idpy, tag, MPI_COMM_WORLD, req(5), ierr)
  call MPI_ISEND(int_neg_send, size_neg_send, MPI_INT, idny, tag, MPI_COMM_WORLD, req(11), ierr)
  tag = 6
  call MPI_IRECV(int_neg_recv, size_neg_recv, MPI_INT, idny, tag, MPI_COMM_WORLD, req(6), ierr)
  call MPI_ISEND(int_pos_send, size_pos_send, MPI_INT, idpy, tag, MPI_COMM_WORLD, req(12), ierr)

! Wait to receive
  call MPI_WAITALL(4, req(3:6), MPI_STATUSES_IGNORE, ierr)

! Organize received arrays

  if (size_pos_recv .gt. 0) then
    int_particles(num_local+1:num_local+size_pos_recv)    = int_pos_recv(:)
    real_particles(:,num_local+1:num_local+size_pos_recv) = real_pos_recv(:,:)
    num_local = num_local+size_pos_recv
  end if

  if (size_neg_recv .gt. 0) then
    int_particles(num_local+1:num_local+size_neg_recv)    = int_neg_recv(:)
    real_particles(:,num_local+1:num_local+size_neg_recv) = real_neg_recv(:,:)
    num_local = num_local+size_neg_recv
  end if

! Clean up before the x direction

  call MPI_WAITALL(6, req(7:12), MPI_STATUSES_IGNORE, ierr)


  deallocate( real_neg_recv )
  deallocate( real_pos_recv )
  deallocate( int_neg_recv )
  deallocate( int_pos_recv )
  deallocate( real_neg_send )
  deallocate( real_pos_send )
  deallocate( int_neg_send )
  deallocate( int_pos_send )
  size_neg_send = 0
  size_neg_recv = 0
  size_pos_send = 0
  size_pos_recv = 0
  neg_index = 0
  pos_index = 0


end if

!!!! X DIRECTION

! If outside local x domain, prepare to send!

if (idpx .ne. myid) then

  do pp = num_local, 1, -1

! Switch memory places so that the ones we want to send end up in the end

    if (real_particles(1,pp) .lt. Lx_min) then
      
      if (pp .ne. num_local-size_neg_send-size_pos_send) then
        real_particles(:,pp) = real_particles(:,pp) + real_particles(:,num_local-size_neg_send-size_pos_send)
        real_particles(:,num_local-size_neg_send-size_pos_send) &
                             = real_particles(:,pp) - real_particles(:,num_local-size_neg_send-size_pos_send)
        real_particles(:,pp) = real_particles(:,pp) - real_particles(:,num_local-size_neg_send-size_pos_send)

        int_particles(pp)    = int_particles(pp)    + int_particles(num_local-size_neg_send-size_pos_send)
        int_particles(num_local-size_neg_send-size_pos_send) &
                             = int_particles(pp)    - int_particles(num_local-size_neg_send-size_pos_send)
        int_particles(pp)    = int_particles(pp)    - int_particles(num_local-size_neg_send-size_pos_send)
      end if

      size_neg_send = size_neg_send + 1

    else if (real_particles(1,pp) .gt. Lx_max) then

      if (pp .ne. num_local-size_neg_send-size_pos_send) then
        real_particles(:,pp) = real_particles(:,pp) + real_particles(:,num_local-size_neg_send-size_pos_send)
        real_particles(:,num_local-size_neg_send-size_pos_send) &
                             = real_particles(:,pp) - real_particles(:,num_local-size_neg_send-size_pos_send)
        real_particles(:,pp) = real_particles(:,pp) - real_particles(:,num_local-size_neg_send-size_pos_send)

        int_particles(pp)    = int_particles(pp)    + int_particles(num_local-size_neg_send-size_pos_send)
        int_particles(num_local-size_neg_send-size_pos_send) &
                             = int_particles(pp)    - int_particles(num_local-size_neg_send-size_pos_send)
        int_particles(pp)    = int_particles(pp)    - int_particles(num_local-size_neg_send-size_pos_send)
      end if

      size_pos_send = size_pos_send + 1    

    end if
  end do

! Send and receive array sizes

  tag = 1
  call MPI_IRECV(size_pos_recv, 1, MPI_INT, idpx, tag, MPI_COMM_WORLD, req(1), ierr)
  call MPI_ISEND(size_neg_send, 1, MPI_INT, idnx, tag, MPI_COMM_WORLD, req(7), ierr)
  tag = 2
  call MPI_IRECV(size_neg_recv, 1, MPI_INT, idnx, tag, MPI_COMM_WORLD, req(2), ierr)
  call MPI_ISEND(size_pos_send, 1, MPI_INT, idpx, tag, MPI_COMM_WORLD, req(8), ierr)

! Form the send matrices

  allocate( real_neg_send(8,size_neg_send) )
  allocate( int_neg_send(size_neg_send) )
  allocate( real_pos_send(8,size_pos_send) )
  allocate( int_pos_send(size_pos_send) )


if (size_pos_send+size_neg_send .gt. 0) then
    
  do pp = num_local, (num_local-size_neg_send-size_pos_send+1), -1

    if (real_particles(1,pp) .lt. Lx_min) then

      neg_index                    = neg_index+1
      int_neg_send(neg_index)      = int_particles(pp)
      real_neg_send(:,neg_index)   = real_particles(:,pp)
      if (real_neg_send(1,neg_index) .lt. xmin) then
        real_neg_send(1,neg_index) = real_neg_send(1,neg_index) + (xmax-xmin)
      end if

    else if (real_particles(1,pp) .gt. Lx_max) then
      pos_index                    = pos_index+1
      int_pos_send(pos_index)      = int_particles(pp)
      real_pos_send(:,pos_index)   = real_particles(:,pp)
      if (real_pos_send(1,pos_index) .gt. xmax) then
        real_pos_send(1,pos_index) = real_pos_send(1,pos_index) - (xmax-xmin)
      end if

    end if
  end do
end if

! Change local size

  num_local = num_local-size_neg_send-size_pos_send

! Wait to receive
  call MPI_WAITALL(2, req(1:2), MPI_STATUSES_IGNORE, ierr)

  allocate( real_neg_recv(8,size_neg_recv) )
  allocate( real_pos_recv(8,size_pos_recv) )
  allocate( int_neg_recv(size_neg_recv) )
  allocate( int_pos_recv(size_pos_recv) )

! Send and receive arrays

  tag = 3
  call MPI_IRECV(real_pos_recv, 8*size_pos_recv, MPI_DOUBLE_PRECISION, idpx, tag, MPI_COMM_WORLD, req(3), ierr)
  call MPI_ISEND(real_neg_send, 8*size_neg_send, MPI_DOUBLE_PRECISION, idnx, tag, MPI_COMM_WORLD, req(9), ierr)
  tag = 4
  call MPI_IRECV(real_neg_recv, 8*size_neg_recv, MPI_DOUBLE_PRECISION, idnx, tag, MPI_COMM_WORLD, req(4), ierr)
  call MPI_ISEND(real_pos_send, 8*size_pos_send, MPI_DOUBLE_PRECISION, idpx, tag, MPI_COMM_WORLD, req(10), ierr)

  tag = 5
  call MPI_IRECV(int_pos_recv, size_pos_recv, MPI_INT, idpx, tag, MPI_COMM_WORLD, req(5), ierr)
  call MPI_ISEND(int_neg_send, size_neg_send, MPI_INT, idnx, tag, MPI_COMM_WORLD, req(11), ierr)
  tag = 6
  call MPI_IRECV(int_neg_recv, size_neg_recv, MPI_INT, idnx, tag, MPI_COMM_WORLD, req(6), ierr)
  call MPI_ISEND(int_pos_send, size_pos_send, MPI_INT, idpx, tag, MPI_COMM_WORLD, req(12), ierr)

! Wait to receive
  call MPI_WAITALL(4, req(3:6), MPI_STATUSES_IGNORE, ierr)

! Organize received arrays

  if (size_pos_recv .gt. 0) then
    int_particles(num_local+1:num_local+size_pos_recv)    = int_pos_recv(:)
    real_particles(:,num_local+1:num_local+size_pos_recv) = real_pos_recv(:,:)
    num_local = num_local+size_pos_recv
  end if

  if (size_neg_recv .gt. 0) then
    int_particles(num_local+1:num_local+size_neg_recv)    = int_neg_recv(:)
    real_particles(:,num_local+1:num_local+size_neg_recv) = real_neg_recv(:,:)
    num_local = num_local+size_neg_recv
  end if

! WAIT FOR ALL SENDS

  call MPI_WAITALL(6, req(7:12), MPI_STATUSES_IGNORE, ierr)

end if



! If we are not using MPI in certain directions we still need to adjust the position for periodic boundary

if ( (idpz .ne. myid) .and. (idpy .ne. myid) .and. (idpx .ne. myid) ) then
  return


else if ( (idpz .eq. myid) .and. (idpy .eq. myid) .and. (idpx .eq. myid) ) then
  do pp = num_local, 1, -1
    if (real_particles(3,pp) .lt. zmin) then
      real_particles(3,pp) = real_particles(3,pp) + (zmax-zmin)
    end if

    if (real_particles(3,pp) .gt. zmax) then
      real_particles(3,pp) = real_particles(3,pp) - (zmax-zmin)
    end if

    if (real_particles(2,pp) .lt. ymin) then
      real_particles(2,pp) = real_particles(2,pp) + (ymax-ymin)
    end if

    if (real_particles(2,pp) .gt. ymax) then
      real_particles(2,pp) = real_particles(2,pp) - (ymax-ymin)
    end if

    if (real_particles(1,pp) .lt. xmin) then
      real_particles(1,pp) = real_particles(1,pp) + (xmax-xmin)
    end if

    if (real_particles(1,pp) .gt. xmax) then
      real_particles(1,pp) = real_particles(1,pp) - (xmax-xmin)
    end if
  end do  

else if ( (idpz .eq. myid) .and. (idpy .eq. myid) ) then
  do pp = num_local, 1, -1
    if (real_particles(3,pp) .lt. zmin) then
      real_particles(3,pp) = real_particles(3,pp) + (zmax-zmin)
    end if

    if (real_particles(3,pp) .gt. zmax) then
      real_particles(3,pp) = real_particles(3,pp) - (zmax-zmin)
    end if

    if (real_particles(2,pp) .lt. ymin) then
      real_particles(2,pp) = real_particles(2,pp) + (ymax-ymin)
    end if

    if (real_particles(2,pp) .gt. ymax) then
      real_particles(2,pp) = real_particles(2,pp) - (ymax-ymin)
    end if
  end do  

else if ( (idpz .eq. myid) .and. (idpx .eq. myid) ) then
  do pp = num_local, 1, -1
    if (real_particles(3,pp) .lt. zmin) then
      real_particles(3,pp) = real_particles(3,pp) + (zmax-zmin)
    end if

    if (real_particles(3,pp) .gt. zmax) then
      real_particles(3,pp) = real_particles(3,pp) - (zmax-zmin)
    end if

    if (real_particles(1,pp) .lt. xmin) then
      real_particles(1,pp) = real_particles(1,pp) + (xmax-xmin)
    end if

    if (real_particles(1,pp) .gt. xmax) then
      real_particles(1,pp) = real_particles(1,pp) - (xmax-xmin)
    end if
  end do  

else if ( (idpy .eq. myid) .and. (idpx .eq. myid) ) then
  do pp = num_local, 1, -1
    if (real_particles(2,pp) .lt. ymin) then
      real_particles(2,pp) = real_particles(2,pp) + (ymax-ymin)
    end if

    if (real_particles(2,pp) .gt. ymax) then
      real_particles(2,pp) = real_particles(2,pp) - (ymax-ymin)
    end if

    if (real_particles(1,pp) .lt. xmin) then
      real_particles(1,pp) = real_particles(1,pp) + (xmax-xmin)
    end if

    if (real_particles(1,pp) .gt. xmax) then
      real_particles(1,pp) = real_particles(1,pp) - (xmax-xmin)
    end if
  end do  

else if (idpz .eq. myid) then
  do pp = num_local, 1, -1
    if (real_particles(3,pp) .lt. zmin) then
      real_particles(3,pp) = real_particles(3,pp) + (zmax-zmin)
    end if

    if (real_particles(3,pp) .gt. zmax) then
      real_particles(3,pp) = real_particles(3,pp) - (zmax-zmin)
    end if
  end do  

else if (idpy .eq. myid) then
  do pp = num_local, 1, -1
    if (real_particles(2,pp) .lt. ymin) then
      real_particles(2,pp) = real_particles(2,pp) + (ymax-ymin)
    end if

    if (real_particles(2,pp) .gt. ymax) then
      real_particles(2,pp) = real_particles(2,pp) - (ymax-ymin)
    end if
  end do  

else if (idpx .eq. myid) then
  do pp = num_local, 1, -1
    if (real_particles(1,pp) .lt. xmin) then
      real_particles(1,pp) = real_particles(1,pp) + (xmax-xmin)
    end if

    if (real_particles(1,pp) .gt. xmax) then
      real_particles(1,pp) = real_particles(1,pp) - (xmax-xmin)
    end if
  end do  

end if

  return
end subroutine ParticlesTransfer
