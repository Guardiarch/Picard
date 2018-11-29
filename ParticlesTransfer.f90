subroutine ParticlesTransfer(particles, &
     Lx_min,Lx_max,Ly_min,Ly_max,Lz_min,Lz_max,&
     xmin,xmax,ymin,ymax,zmin,zmax,max_per_proc, &
     num_local,iprocs,jprocs,kprocs,myid)

  use SpecificTypes
  use mpi
  implicit none
  


! Parameters
  integer, intent(in) :: max_per_proc, iprocs, jprocs, kprocs, myid
  real*8, intent(in) :: Lx_min,Lx_max,Ly_min,Ly_max,Lz_min,Lz_max, &
       xmin,xmax,ymin,ymax,zmin,zmax
  integer, intent(inout) :: num_local
  type(particlearrays) particles
 
  
! Local variables
  integer :: ii, jj, kk, pp, neg_index, pos_index
  integer :: size_neg_send, size_neg_recv, size_pos_send, size_pos_recv 
  integer :: idpx, idnx, idpy, idny, idpz, idnz

! For the implementation of non-blocking communication
  integer (kind=MPI_ADDRESS_KIND) :: offsets(2)
  integer :: tag, ierr, req(12)
  type(particlearrays) neg_send, neg_recv, pos_send, pos_recv
  integer blockcounts(2), oldtypes(2), &
       b_type_neg_send, b_type_pos_send, b_type_neg_recv, b_type_pos_recv

  
  idpx = myid + &
       (mod(myid+iprocs+1,             iprocs) - mod(myid+iprocs, iprocs))
  idnx = myid + &
       (mod(myid+iprocs-1,             iprocs) - mod(myid+iprocs, iprocs))
  idpy = myid + &
       (mod(myid+iprocs*jprocs+iprocs, iprocs*jprocs) - &
       mod(myid+iprocs*jprocs, iprocs*jprocs))
  idny = myid + &
       (mod(myid+iprocs*jprocs-iprocs, iprocs*jprocs) - &
       mod(myid+iprocs*jprocs, iprocs*jprocs))
  idpz = myid + &
       (mod(myid+iprocs*jprocs*kprocs+iprocs*jprocs, iprocs*jprocs*kprocs) &
       - mod(myid+iprocs*jprocs*kprocs, iprocs*jprocs*kprocs))
  idnz = myid + &
       (mod(myid+iprocs*jprocs*kprocs-iprocs*jprocs, iprocs*jprocs*kprocs) &
       - mod(myid+iprocs*jprocs*kprocs, iprocs*jprocs*kprocs))


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

    if (particles%coordinates(3,pp) .lt. Lz_min) then

      if (pp .ne. num_local-size_neg_send-size_pos_send) then
         particles%coordinates(:,pp) = particles%coordinates(:,pp) + &
              particles%coordinates(:,num_local-size_neg_send-size_pos_send)
         particles%coordinates(:,num_local-size_neg_send-size_pos_send) &
              = particles%coordinates(:,pp) - &
              particles%coordinates(:,num_local-size_neg_send-size_pos_send)
         particles%coordinates(:,pp) = particles%coordinates(:,pp) - &
              particles%coordinates(:,num_local-size_neg_send-size_pos_send)

         particles%species(pp) = particles%species(pp) + &
              particles%species(num_local-size_neg_send-size_pos_send)
         particles%species(num_local-size_neg_send-size_pos_send) = &
              particles%species(pp) - &
              particles%species(num_local-size_neg_send-size_pos_send)
         particles%species(pp) = particles%species(pp) - &
              particles%species(num_local-size_neg_send-size_pos_send)
      end if

      size_neg_send = size_neg_send + 1

   else if (particles%coordinates(3,pp) .gt. Lz_max) then

      if (pp .ne. num_local-size_neg_send-size_pos_send) then
         particles%coordinates(:,pp) = particles%coordinates(:,pp) + &
              particles%coordinates(:,num_local-size_neg_send-size_pos_send)
         particles%coordinates(:,num_local-size_neg_send-size_pos_send) = &
              particles%coordinates(:,pp) - &
              particles%coordinates(:,num_local-size_neg_send-size_pos_send)
         particles%coordinates(:,pp) = particles%coordinates(:,pp) - &
              particles%coordinates(:,num_local-size_neg_send-size_pos_send)

         particles%species(pp) = particles%species(pp) + &
              particles%species(num_local-size_neg_send-size_pos_send)
         particles%species(num_local-size_neg_send-size_pos_send) = &
              particles%species(pp) - &
              particles%species(num_local-size_neg_send-size_pos_send)
         particles%species(pp) = particles%species(pp) - &
              particles%species(num_local-size_neg_send-size_pos_send)
      end if

      size_pos_send = size_pos_send + 1    

     end if
  end do

! Send and receive array sizes

  tag = 1
  call MPI_IRECV(size_pos_recv,1,MPI_INT,idpz,tag,MPI_COMM_WORLD,req(1),ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 1, ierr=',ierr,' myid=',myid
  end if
  call MPI_ISEND(size_neg_send,1,MPI_INT,idnz,tag,MPI_COMM_WORLD,req(7),ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 2, ierr=',ierr,' myid=',myid
  end if
  tag = 2
  call MPI_IRECV(size_neg_recv,1,MPI_INT,idnz,tag,MPI_COMM_WORLD,req(2),ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 3, ierr=',ierr,' myid=',myid
  end if
  call MPI_ISEND(size_pos_send,1,MPI_INT,idpz,tag,MPI_COMM_WORLD,req(8),ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 4, ierr=',ierr,' myid=',myid
  end if

! Allocate the send buffers
  allocate( neg_send%species(size_neg_send) )
  allocate( neg_send%coordinates(6, size_neg_send) )
  allocate( pos_send%species(size_pos_send) )
  allocate( pos_send%coordinates(6, size_pos_send) )

if (size_pos_send+size_neg_send .gt. 0) then
  do pp = num_local, (num_local-size_neg_send-size_pos_send+1), -1

    if (particles%coordinates(3,pp) .lt. Lz_min) then
      neg_index                    = neg_index+1
      neg_send%species(neg_index)  = particles%species(pp)
      neg_send%coordinates(:,neg_index)   = particles%coordinates(:,pp)
      if (neg_send%coordinates(3,neg_index) .lt. zmin) then
         neg_send%coordinates(3,neg_index)=neg_send%coordinates(3,neg_index)+ &
              (zmax-zmin)
      end if

    else if (particles%coordinates(3,pp) .gt. Lz_max) then
      pos_index                    = pos_index+1
      pos_send%species(pos_index)   = particles%species(pp)
      pos_send%coordinates(:,pos_index) = particles%coordinates(:,pp)
      if (pos_send%coordinates(3,pos_index) .gt. zmax) then
         pos_send%coordinates(3,pos_index)=pos_send%coordinates(3,pos_index)- &
              (zmax-zmin)
      end if

    end if
  end do
end if

! Change local size

  num_local = num_local-size_neg_send-size_pos_send

! Wait to receive
  call MPI_WAITALL(2, req(1:2), MPI_STATUSES_IGNORE, ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 5, ierr=',ierr,' myid=',myid
  end if

! Allocate receive buffers
  allocate( neg_recv%species(size_neg_recv) )
  allocate( neg_recv%coordinates(6, size_neg_recv) )
  allocate( pos_recv%species(size_pos_recv) )
  allocate( pos_recv%coordinates(6, size_pos_recv) )

! Create new MPI datatypes for the buffers
  oldtypes(1) = MPI_INTEGER
  oldtypes(2) = MPI_DOUBLE_PRECISION
  ! negative receive
  blockcounts(1) = size_neg_recv
  blockcounts(2) = size_neg_recv*6
  call MPI_GET_ADDRESS(neg_recv%species, offsets(1), ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 6, ierr=',ierr,' myid=',myid
  end if
  call MPI_GET_ADDRESS(neg_recv%coordinates, offsets(2), ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 7, ierr=',ierr,' myid=',myid
  end if
  offsets = offsets-offsets(1)
  call MPI_TYPE_CREATE_STRUCT &
       (2,blockcounts,offsets,oldtypes,b_type_neg_recv,ierr)
    if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 8, ierr=',ierr,' myid=',myid
  end if
  call MPI_TYPE_COMMIT(b_type_neg_recv, ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 9, ierr=',ierr,' myid=',myid
  end if
  ! positive receive
  blockcounts(1) = size_pos_recv
  blockcounts(2) = size_pos_recv*6
  call MPI_GET_ADDRESS(pos_recv%species, offsets(1), ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 10, ierr=',ierr,' myid=',myid
  end if
  call MPI_GET_ADDRESS(pos_recv%coordinates, offsets(2), ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 11, ierr=',ierr,' myid=',myid
  end if
  offsets = offsets-offsets(1)
  call MPI_TYPE_CREATE_STRUCT &
       (2,blockcounts,offsets,oldtypes,b_type_pos_recv,ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 12, ierr=',ierr,' myid=',myid
  end if
  call MPI_TYPE_COMMIT(b_type_pos_recv, ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 13, ierr=',ierr,' myid=',myid
  end if
  ! negative send
  blockcounts(1) = size_neg_send
  blockcounts(2) = size_neg_send*6
  call MPI_GET_ADDRESS(neg_send%species, offsets(1), ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 14, ierr=',ierr,' myid=',myid
  end if
  call MPI_GET_ADDRESS(neg_send%coordinates, offsets(2), ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 15, ierr=',ierr,' myid=',myid
  end if
  offsets = offsets-offsets(1)
  call MPI_TYPE_CREATE_STRUCT &
       (2,blockcounts,offsets,oldtypes,b_type_neg_send,ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 16, ierr=',ierr,' myid=',myid
  end if
  call MPI_TYPE_COMMIT(b_type_neg_send, ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 17, ierr=',ierr,' myid=',myid
  end if
  ! positive send
  blockcounts(1) = size_pos_send
  blockcounts(2) = size_pos_send*6
  call MPI_GET_ADDRESS(pos_send%species, offsets(1), ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 18, ierr=',ierr,' myid=',myid
  end if
  call MPI_GET_ADDRESS(pos_send%coordinates, offsets(2), ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 19, ierr=',ierr,' myid=',myid
  end if
  offsets = offsets-offsets(1)
  call MPI_TYPE_CREATE_STRUCT &
       (2,blockcounts,offsets,oldtypes,b_type_pos_send,ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 20, ierr=',ierr,' myid=',myid
  end if
  call MPI_TYPE_COMMIT(b_type_pos_send, ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 21, ierr=',ierr,' myid=',myid
  end if
  

! Send and receive the buffers
  tag = 3
  call MPI_IRECV(pos_recv%species, 1, b_type_pos_recv, idpz, &
       tag, MPI_COMM_WORLD, req(3), ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 22, ierr=',ierr,' myid=',myid
  end if
  call MPI_ISEND(neg_send%species,1, b_type_neg_send, idnz, &
       tag, MPI_COMM_WORLD, req(9), ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 23, ierr=',ierr,' myid=',myid
  end if
  tag = 4
  call MPI_IRECV(neg_recv%species, 1, b_type_neg_recv, idnz, &
       tag, MPI_COMM_WORLD, req(4), ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 24, ierr=',ierr,' myid=',myid
  end if
  call MPI_ISEND(pos_send%species,1, b_type_pos_send, idpz, &
       tag, MPI_COMM_WORLD, req(10), ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 25, ierr=',ierr,' myid=',myid
  end if

! Wait to receive
  call MPI_WAITALL(2, req(3:4), MPI_STATUSES_IGNORE, ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 26, ierr=',ierr,' myid=',myid
  end if

! Organize received arrays

  if (size_pos_recv .gt. 0) then
     particles%species(num_local+1:num_local+size_pos_recv)=pos_recv%species(:)
     particles%coordinates(:,num_local+1:num_local+size_pos_recv) = &
          pos_recv%coordinates(:,:)
     num_local = num_local+size_pos_recv
  end if

  if (size_neg_recv .gt. 0) then
     particles%species(num_local+1:num_local+size_neg_recv)=neg_recv%species(:)
     particles%coordinates(:,num_local+1:num_local+size_neg_recv) = &
          neg_recv%coordinates(:,:)
     num_local = num_local+size_neg_recv
  end if

! Clean up before the y direction
  call MPI_WAITALL(4, req(7:10), MPI_STATUSES_IGNORE, ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 27, ierr=',ierr,' myid=',myid
  end if
  call MPI_TYPE_FREE(b_type_neg_recv, ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 28, ierr=',ierr,' myid=',myid
  end if
  call MPI_TYPE_FREE(b_type_pos_recv, ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 29, ierr=',ierr,' myid=',myid
  end if
  call MPI_TYPE_FREE(b_type_neg_send, ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 30, ierr=',ierr,' myid=',myid
  end if
  call MPI_TYPE_FREE(b_type_pos_send, ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 31, ierr=',ierr,' myid=',myid
  end if
  deallocate( neg_recv%species )
  deallocate( neg_recv%coordinates )
  deallocate( pos_recv%species )
  deallocate( pos_recv%coordinates )
  deallocate( neg_send%species )
  deallocate( neg_send%coordinates )
  deallocate( pos_send%species )
  deallocate( pos_send%coordinates )
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

    if (particles%coordinates(2,pp) .lt. Ly_min) then

      if (pp .ne. num_local-size_neg_send-size_pos_send) then
         particles%coordinates(:,pp) = particles%coordinates(:,pp) + &
              particles%coordinates(:,num_local-size_neg_send-size_pos_send)
         particles%coordinates(:,num_local-size_neg_send-size_pos_send) = &
              particles%coordinates(:,pp) - &
              particles%coordinates(:,num_local-size_neg_send-size_pos_send)
         particles%coordinates(:,pp) = particles%coordinates(:,pp) - &
              particles%coordinates(:,num_local-size_neg_send-size_pos_send)

         particles%species(pp) = particles%species(pp) + &
              particles%species(num_local-size_neg_send-size_pos_send)
         particles%species(num_local-size_neg_send-size_pos_send) = &
              particles%species(pp) - &
              particles%species(num_local-size_neg_send-size_pos_send)
         particles%species(pp) = particles%species(pp) - &
              particles%species(num_local-size_neg_send-size_pos_send)
      end if

      size_neg_send = size_neg_send + 1

    else if (particles%coordinates(2,pp) .gt. Ly_max) then

       if (pp .ne. num_local-size_neg_send-size_pos_send) then
          particles%coordinates(:,pp) = particles%coordinates(:,pp) + &
               particles%coordinates(:,num_local-size_neg_send-size_pos_send)
          particles%coordinates(:,num_local-size_neg_send-size_pos_send) = &
               particles%coordinates(:,pp) - &
               particles%coordinates(:,num_local-size_neg_send-size_pos_send)
          particles%coordinates(:,pp) = particles%coordinates(:,pp) - &
               particles%coordinates(:,num_local-size_neg_send-size_pos_send)

          particles%species(pp) = particles%species(pp) + &
               particles%species(num_local-size_neg_send-size_pos_send)
          particles%species(num_local-size_neg_send-size_pos_send) = &
               particles%species(pp) - &
               particles%species(num_local-size_neg_send-size_pos_send)
          particles%species(pp) = particles%species(pp) - &
               particles%species(num_local-size_neg_send-size_pos_send)
      end if

      size_pos_send = size_pos_send + 1    

    end if
  end do

! Send and receive array sizes

  tag = 1
  call MPI_IRECV(size_pos_recv,1,MPI_INT,idpy,tag,MPI_COMM_WORLD,req(1),ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 32, ierr=',ierr,' myid=',myid
  end if
  call MPI_ISEND(size_neg_send,1,MPI_INT,idny,tag,MPI_COMM_WORLD,req(7),ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 33, ierr=',ierr,' myid=',myid
  end if
  tag = 2
  call MPI_IRECV(size_neg_recv,1,MPI_INT,idny,tag,MPI_COMM_WORLD,req(2),ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 34, ierr=',ierr,' myid=',myid
  end if
  call MPI_ISEND(size_pos_send,1,MPI_INT,idpy,tag,MPI_COMM_WORLD,req(8),ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 35, ierr=',ierr,' myid=',myid
  end if

! Allocate the send buffers
  allocate( neg_send%species(size_neg_send) )
  allocate( neg_send%coordinates(6, size_neg_send) )
  allocate( pos_send%species(size_pos_send) )
  allocate( pos_send%coordinates(6, size_pos_send) )

if (size_pos_send+size_neg_send .gt. 0) then
  do pp = num_local, (num_local-size_neg_send-size_pos_send+1), -1

    if (particles%coordinates(2,pp) .lt. Ly_min) then
      neg_index                    = neg_index+1
      neg_send%species(neg_index)  = particles%species(pp)
      neg_send%coordinates(:,neg_index)   = particles%coordinates(:,pp)
      if (neg_send%coordinates(2,neg_index) .lt. ymin) then
         neg_send%coordinates(2,neg_index)=neg_send%coordinates(2,neg_index)+ &
              (ymax-ymin)
      end if

    else if (particles%coordinates(2,pp) .gt. Ly_max) then
      pos_index                    = pos_index+1
      pos_send%species(pos_index)  = particles%species(pp)
      pos_send%coordinates(:,pos_index)   = particles%coordinates(:,pp)
      if (pos_send%coordinates(2,pos_index) .gt. ymax) then
         pos_send%coordinates(2,pos_index)=pos_send%coordinates(2,pos_index)- &
              (ymax-ymin)
      end if

    end if
  end do
end if

! Change local size

  num_local = num_local-size_neg_send-size_pos_send

! Wait to receive
  call MPI_WAITALL(2, req(1:2), MPI_STATUSES_IGNORE, ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 36, ierr=',ierr,' myid=',myid
  end if

! Allocate receive buffers
  allocate( neg_recv%species(size_neg_recv) )
  allocate( neg_recv%coordinates(6, size_neg_recv) )
  allocate( pos_recv%species(size_pos_recv) )
  allocate( pos_recv%coordinates(6, size_pos_recv) )

! Create new MPI datatypes for the buffers
  oldtypes(1) = MPI_INTEGER
  oldtypes(2) = MPI_DOUBLE_PRECISION
  ! negative receive
  blockcounts(1) = size_neg_recv
  blockcounts(2) = size_neg_recv*6
  call MPI_GET_ADDRESS(neg_recv%species, offsets(1), ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 37, ierr=',ierr,' myid=',myid
  end if
  call MPI_GET_ADDRESS(neg_recv%coordinates, offsets(2), ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 38, ierr=',ierr,' myid=',myid
  end if
  offsets = offsets-offsets(1)
  call MPI_TYPE_CREATE_STRUCT &
       (2,blockcounts,offsets,oldtypes,b_type_neg_recv,ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 39, ierr=',ierr,' myid=',myid
  end if
  call MPI_TYPE_COMMIT(b_type_neg_recv, ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 40, ierr=',ierr,' myid=',myid
  end if
  ! positive receive
  blockcounts(1) = size_pos_recv
  blockcounts(2) = size_pos_recv*6
  call MPI_GET_ADDRESS(pos_recv%species, offsets(1), ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 41, ierr=',ierr,' myid=',myid
  end if
  call MPI_GET_ADDRESS(pos_recv%coordinates, offsets(2), ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 42, ierr=',ierr,' myid=',myid
  end if
  offsets = offsets-offsets(1)
  call MPI_TYPE_CREATE_STRUCT &
       (2,blockcounts,offsets,oldtypes,b_type_pos_recv,ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 43, ierr=',ierr,' myid=',myid
  end if
  call MPI_TYPE_COMMIT(b_type_pos_recv, ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 44, ierr=',ierr,' myid=',myid
  end if
  ! negative send
  blockcounts(1) = size_neg_send
  blockcounts(2) = size_neg_send*6
  call MPI_GET_ADDRESS(neg_send%species, offsets(1), ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 45, ierr=',ierr,' myid=',myid
  end if
  call MPI_GET_ADDRESS(neg_send%coordinates, offsets(2), ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 46, ierr=',ierr,' myid=',myid
  end if
  offsets = offsets-offsets(1)
  call MPI_TYPE_CREATE_STRUCT &
       (2,blockcounts,offsets,oldtypes,b_type_neg_send,ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 47, ierr=',ierr,' myid=',myid
  end if
  call MPI_TYPE_COMMIT(b_type_neg_send, ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 48, ierr=',ierr,' myid=',myid
  end if
  ! positive send
  blockcounts(1) = size_pos_send
  blockcounts(2) = size_pos_send*6
  call MPI_GET_ADDRESS(pos_send%species, offsets(1), ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 49, ierr=',ierr,' myid=',myid
  end if
  call MPI_GET_ADDRESS(pos_send%coordinates, offsets(2), ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 50, ierr=',ierr,' myid=',myid
  end if
  offsets = offsets-offsets(1)
  call MPI_TYPE_CREATE_STRUCT &
       (2,blockcounts,offsets,oldtypes,b_type_pos_send,ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 51, ierr=',ierr,' myid=',myid
  end if
  call MPI_TYPE_COMMIT(b_type_pos_send, ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 52, ierr=',ierr,' myid=',myid
  end if


! Send and receive the buffers
  tag = 3
  call MPI_IRECV(pos_recv%species, 1, b_type_pos_recv, idpy, &
       tag, MPI_COMM_WORLD, req(3), ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 53, ierr=',ierr,' myid=',myid
  end if
  call MPI_ISEND(neg_send%species,1, b_type_neg_send, idny, &
       tag, MPI_COMM_WORLD, req(9), ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 54, ierr=',ierr,' myid=',myid
  end if
  tag = 4
  call MPI_IRECV(neg_recv%species, 1, b_type_neg_recv, idny, &
       tag, MPI_COMM_WORLD, req(4), ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 55, ierr=',ierr,' myid=',myid
  end if
  call MPI_ISEND(pos_send%species,1, b_type_pos_send, idpy, &
       tag, MPI_COMM_WORLD, req(10), ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 56, ierr=',ierr,' myid=',myid
  end if

! Wait to receive
  call MPI_WAITALL(2, req(3:4), MPI_STATUSES_IGNORE, ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 57, ierr=',ierr,' myid=',myid
  end if

! Organize received arrays
  if (size_pos_recv .gt. 0) then
    particles%species(num_local+1:num_local+size_pos_recv)=pos_recv%species(:)
    particles%coordinates(:,num_local+1:num_local+size_pos_recv) = &
         pos_recv%coordinates(:,:)
    num_local = num_local+size_pos_recv
  end if

  if (size_neg_recv .gt. 0) then
    particles%species(num_local+1:num_local+size_neg_recv)=neg_recv%species(:)
    particles%coordinates(:,num_local+1:num_local+size_neg_recv) = &
         neg_recv%coordinates(:,:)
    num_local = num_local+size_neg_recv
  end if

! Clean up before the x direction
  call MPI_WAITALL(4, req(7:10), MPI_STATUSES_IGNORE, ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 58, ierr=',ierr,' myid=',myid
  end if
  call MPI_TYPE_FREE(b_type_neg_recv, ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 59, ierr=',ierr,' myid=',myid
  end if
  call MPI_TYPE_FREE(b_type_pos_recv, ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 60, ierr=',ierr,' myid=',myid
  end if
  call MPI_TYPE_FREE(b_type_neg_send, ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 61, ierr=',ierr,' myid=',myid
  end if
  call MPI_TYPE_FREE(b_type_pos_send, ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 62, ierr=',ierr,' myid=',myid
  end if
  deallocate( neg_recv%species )
  deallocate( neg_recv%coordinates )
  deallocate( pos_recv%species )
  deallocate( pos_recv%coordinates )
  deallocate( neg_send%species )
  deallocate( neg_send%coordinates )
  deallocate( pos_send%species )
  deallocate( pos_send%coordinates )
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

    if (particles%coordinates(1,pp) .lt. Lx_min) then
      
      if (pp .ne. num_local-size_neg_send-size_pos_send) then
         particles%coordinates(:,pp) = particles%coordinates(:,pp) + &
              particles%coordinates(:,num_local-size_neg_send-size_pos_send)
         particles%coordinates(:,num_local-size_neg_send-size_pos_send) = &
              particles%coordinates(:,pp) - &
              particles%coordinates(:,num_local-size_neg_send-size_pos_send)
         particles%coordinates(:,pp) = particles%coordinates(:,pp) - &
              particles%coordinates(:,num_local-size_neg_send-size_pos_send)

         particles%species(pp) = particles%species(pp) + &
              particles%species(num_local-size_neg_send-size_pos_send)
        particles%species(num_local-size_neg_send-size_pos_send) = &
             particles%species(pp) - &
             particles%species(num_local-size_neg_send-size_pos_send)
        particles%species(pp) = particles%species(pp) - &
             particles%species(num_local-size_neg_send-size_pos_send)
      end if

      size_neg_send = size_neg_send + 1

    else if (particles%coordinates(1,pp) .gt. Lx_max) then

      if (pp .ne. num_local-size_neg_send-size_pos_send) then
         particles%coordinates(:,pp) = particles%coordinates(:,pp) + &
              particles%coordinates(:,num_local-size_neg_send-size_pos_send)
         particles%coordinates(:,num_local-size_neg_send-size_pos_send) = &
              particles%coordinates(:,pp) - &
              particles%coordinates(:,num_local-size_neg_send-size_pos_send)
         particles%coordinates(:,pp) = particles%coordinates(:,pp) - &
              particles%coordinates(:,num_local-size_neg_send-size_pos_send)

         particles%species(pp) = particles%species(pp) + &
              particles%species(num_local-size_neg_send-size_pos_send)
         particles%species(num_local-size_neg_send-size_pos_send) = &
              particles%species(pp) - &
              particles%species(num_local-size_neg_send-size_pos_send)
         particles%species(pp) = particles%species(pp) - &
              particles%species(num_local-size_neg_send-size_pos_send)
      end if

      size_pos_send = size_pos_send + 1    

    end if
  end do

! Send and receive array sizes
  tag = 1
  call MPI_IRECV(size_pos_recv,1,MPI_INT,idpx,tag,MPI_COMM_WORLD,req(1),ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 63, ierr=',ierr,' myid=',myid
  end if
  call MPI_ISEND(size_neg_send,1,MPI_INT,idnx,tag,MPI_COMM_WORLD,req(7),ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 64, ierr=',ierr,' myid=',myid
  end if
  tag = 2
  call MPI_IRECV(size_neg_recv,1,MPI_INT,idnx,tag,MPI_COMM_WORLD,req(2),ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 65, ierr=',ierr,' myid=',myid
  end if
  call MPI_ISEND(size_pos_send,1,MPI_INT,idpx,tag,MPI_COMM_WORLD,req(8),ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 66, ierr=',ierr,' myid=',myid
  end if

! Allocate the send buffers
  allocate( neg_send%species(size_neg_send) )
  allocate( neg_send%coordinates(6, size_neg_send) )
  allocate( pos_send%species(size_pos_send) )
  allocate( pos_send%coordinates(6, size_pos_send) )

if (size_pos_send+size_neg_send .gt. 0) then
    
  do pp = num_local, (num_local-size_neg_send-size_pos_send+1), -1

    if (particles%coordinates(1,pp) .lt. Lx_min) then

      neg_index                    = neg_index+1
      neg_send%species(neg_index)  = particles%species(pp)
      neg_send%coordinates(:,neg_index) = particles%coordinates(:,pp)
      if (neg_send%coordinates(1,neg_index) .lt. xmin) then
         neg_send%coordinates(1,neg_index)=neg_send%coordinates(1,neg_index)+ &
              (xmax-xmin)
      end if

    else if (particles%coordinates(1,pp) .gt. Lx_max) then
      pos_index                    = pos_index+1
      pos_send%species(pos_index)  = particles%species(pp)
      pos_send%coordinates(:,pos_index) = particles%coordinates(:,pp)
      if (pos_send%coordinates(1,pos_index) .gt. xmax) then
         pos_send%coordinates(1,pos_index)=pos_send%coordinates(1,pos_index)- &
              (xmax-xmin)
      end if

    end if
  end do
end if

! Change local size
  num_local = num_local-size_neg_send-size_pos_send

! Wait to receive
  call MPI_WAITALL(2, req(1:2), MPI_STATUSES_IGNORE, ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 67, ierr=',ierr,' myid=',myid
  end if

! Allocate receive buffers
  allocate( neg_recv%species(size_neg_recv) )
  allocate( neg_recv%coordinates(6, size_neg_recv) )
  allocate( pos_recv%species(size_pos_recv) )
  allocate( pos_recv%coordinates(6, size_pos_recv) )

! Create new MPI datatypes for the buffers
  oldtypes(1) = MPI_INTEGER
  oldtypes(2) = MPI_DOUBLE_PRECISION
  ! negative receive
  blockcounts(1) = size_neg_recv
  blockcounts(2) = size_neg_recv*6
  call MPI_GET_ADDRESS(neg_recv%species, offsets(1), ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 68, ierr=',ierr,' myid=',myid
  end if
  call MPI_GET_ADDRESS(neg_recv%coordinates, offsets(2), ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 69, ierr=',ierr,' myid=',myid
  end if
  offsets = offsets-offsets(1)
  call MPI_TYPE_CREATE_STRUCT &
       (2,blockcounts,offsets,oldtypes,b_type_neg_recv,ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 70, ierr=',ierr,' myid=',myid
  end if
  call MPI_TYPE_COMMIT(b_type_neg_recv, ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 71, ierr=',ierr,' myid=',myid
  end if
  ! positive receive
  blockcounts(1) = size_pos_recv
  blockcounts(2) = size_pos_recv*6
  call MPI_GET_ADDRESS(pos_recv%species, offsets(1), ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 72, ierr=',ierr,' myid=',myid
  end if
  call MPI_GET_ADDRESS(pos_recv%coordinates, offsets(2), ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 73, ierr=',ierr,' myid=',myid
  end if
  offsets = offsets-offsets(1)
  call MPI_TYPE_CREATE_STRUCT &
       (2,blockcounts,offsets,oldtypes,b_type_pos_recv,ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 74, ierr=',ierr,' myid=',myid
  end if
  call MPI_TYPE_COMMIT(b_type_pos_recv, ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 75, ierr=',ierr,' myid=',myid
  end if
  ! negative send
  blockcounts(1) = size_neg_send
  blockcounts(2) = size_neg_send*6
  call MPI_GET_ADDRESS(neg_send%species, offsets(1), ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 76, ierr=',ierr,' myid=',myid
  end if
  call MPI_GET_ADDRESS(neg_send%coordinates, offsets(2), ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 77, ierr=',ierr,' myid=',myid
  end if
  offsets = offsets-offsets(1)
  call MPI_TYPE_CREATE_STRUCT &
       (2,blockcounts,offsets,oldtypes,b_type_neg_send,ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 78, ierr=',ierr,' myid=',myid
  end if
  call MPI_TYPE_COMMIT(b_type_neg_send, ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 79, ierr=',ierr,' myid=',myid
  end if
  ! positive send
  blockcounts(1) = size_pos_send
  blockcounts(2) = size_pos_send*6
  call MPI_GET_ADDRESS(pos_send%species, offsets(1), ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 80, ierr=',ierr,' myid=',myid
  end if
  call MPI_GET_ADDRESS(pos_send%coordinates, offsets(2), ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 81, ierr=',ierr,' myid=',myid
  end if
  offsets = offsets-offsets(1)
  call MPI_TYPE_CREATE_STRUCT &
       (2,blockcounts,offsets,oldtypes,b_type_pos_send,ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 82, ierr=',ierr,' myid=',myid
  end if
  call MPI_TYPE_COMMIT(b_type_pos_send, ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 83, ierr=',ierr,' myid=',myid
  end if


! Send and receive the buffers
  tag = 3
  call MPI_IRECV(pos_recv%species, 1, b_type_pos_recv, idpx, &
       tag, MPI_COMM_WORLD, req(3), ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 84, ierr=',ierr,' myid=',myid
  end if
  call MPI_ISEND(neg_send%species,1, b_type_neg_send, idnx, &
       tag, MPI_COMM_WORLD, req(9), ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 85, ierr=',ierr,' myid=',myid
  end if
  tag = 4
  call MPI_IRECV(neg_recv%species, 1, b_type_neg_recv, idnx, &
       tag, MPI_COMM_WORLD, req(4), ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 86, ierr=',ierr,' myid=',myid
  end if
  call MPI_ISEND(pos_send%species,1, b_type_pos_send, idpx, &
       tag, MPI_COMM_WORLD, req(10), ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 87, ierr=',ierr,' myid=',myid
  end if

! Wait to receive
  call MPI_WAITALL(2, req(3:4), MPI_STATUSES_IGNORE, ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 88, ierr=',ierr,' myid=',myid
  end if

! Organize received arrays

  if (size_pos_recv .gt. 0) then
     particles%species(num_local+1:num_local+size_pos_recv)=pos_recv%species(:)
     particles%coordinates(:,num_local+1:num_local+size_pos_recv) = &
          pos_recv%coordinates(:,:)
     num_local = num_local+size_pos_recv
  end if

  if (size_neg_recv .gt. 0) then
     particles%species(num_local+1:num_local+size_neg_recv)=neg_recv%species(:)
     particles%coordinates(:,num_local+1:num_local+size_neg_recv) = &
          neg_recv%coordinates(:,:)
     num_local = num_local+size_neg_recv
  end if

! WAIT FOR ALL SENDS
  call MPI_WAITALL(4, req(7:10), MPI_STATUSES_IGNORE, ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 89, ierr=',ierr,' myid=',myid
  end if
  call MPI_TYPE_FREE(b_type_neg_recv, ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 90, ierr=',ierr,' myid=',myid
  end if
  call MPI_TYPE_FREE(b_type_pos_recv, ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 91, ierr=',ierr,' myid=',myid
  end if
  call MPI_TYPE_FREE(b_type_neg_send, ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 92, ierr=',ierr,' myid=',myid
  end if
  call MPI_TYPE_FREE(b_type_pos_send, ierr)
  if (ierr>0) then
     write (*,*) 'ParticlesTransfer, 93, ierr=',ierr,' myid=',myid
  end if
  deallocate( neg_recv%species )
  deallocate( neg_recv%coordinates )
  deallocate( pos_recv%species )
  deallocate( pos_recv%coordinates )
  deallocate( neg_send%species )
  deallocate( neg_send%coordinates )
  deallocate( pos_send%species )
  deallocate( pos_send%coordinates )

end if



! If we are not using MPI in certain directions we still need to
! adjust the position for periodic boundary

if ( (idpz .ne. myid) .and. (idpy .ne. myid) .and. (idpx .ne. myid) ) then
  return


else if ( (idpz .eq. myid) .and. (idpy .eq. myid) .and. (idpx .eq. myid) ) then
  do pp = num_local, 1, -1
    if (particles%coordinates(3,pp) .lt. zmin) then
      particles%coordinates(3,pp) = particles%coordinates(3,pp) + (zmax-zmin)
    end if

    if (particles%coordinates(3,pp) .gt. zmax) then
      particles%coordinates(3,pp) = particles%coordinates(3,pp) - (zmax-zmin)
    end if

    if (particles%coordinates(2,pp) .lt. ymin) then
      particles%coordinates(2,pp) = particles%coordinates(2,pp) + (ymax-ymin)
    end if

    if (particles%coordinates(2,pp) .gt. ymax) then
      particles%coordinates(2,pp) = particles%coordinates(2,pp) - (ymax-ymin)
    end if

    if (particles%coordinates(1,pp) .lt. xmin) then
      particles%coordinates(1,pp) = particles%coordinates(1,pp) + (xmax-xmin)
    end if

    if (particles%coordinates(1,pp) .gt. xmax) then
      particles%coordinates(1,pp) = particles%coordinates(1,pp) - (xmax-xmin)
    end if
  end do  

else if ( (idpz .eq. myid) .and. (idpy .eq. myid) ) then
  do pp = num_local, 1, -1
    if (particles%coordinates(3,pp) .lt. zmin) then
      particles%coordinates(3,pp) = particles%coordinates(3,pp) + (zmax-zmin)
    end if

    if (particles%coordinates(3,pp) .gt. zmax) then
      particles%coordinates(3,pp) = particles%coordinates(3,pp) - (zmax-zmin)
    end if

    if (particles%coordinates(2,pp) .lt. ymin) then
      particles%coordinates(2,pp) = particles%coordinates(2,pp) + (ymax-ymin)
    end if

    if (particles%coordinates(2,pp) .gt. ymax) then
      particles%coordinates(2,pp) = particles%coordinates(2,pp) - (ymax-ymin)
    end if
  end do  

else if ( (idpz .eq. myid) .and. (idpx .eq. myid) ) then
  do pp = num_local, 1, -1
    if (particles%coordinates(3,pp) .lt. zmin) then
      particles%coordinates(3,pp) = particles%coordinates(3,pp) + (zmax-zmin)
    end if

    if (particles%coordinates(3,pp) .gt. zmax) then
      particles%coordinates(3,pp) = particles%coordinates(3,pp) - (zmax-zmin)
    end if

    if (particles%coordinates(1,pp) .lt. xmin) then
      particles%coordinates(1,pp) = particles%coordinates(1,pp) + (xmax-xmin)
    end if

    if (particles%coordinates(1,pp) .gt. xmax) then
      particles%coordinates(1,pp) = particles%coordinates(1,pp) - (xmax-xmin)
    end if
  end do  

else if ( (idpy .eq. myid) .and. (idpx .eq. myid) ) then
  do pp = num_local, 1, -1
    if (particles%coordinates(2,pp) .lt. ymin) then
      particles%coordinates(2,pp) = particles%coordinates(2,pp) + (ymax-ymin)
    end if

    if (particles%coordinates(2,pp) .gt. ymax) then
      particles%coordinates(2,pp) = particles%coordinates(2,pp) - (ymax-ymin)
    end if

    if (particles%coordinates(1,pp) .lt. xmin) then
      particles%coordinates(1,pp) = particles%coordinates(1,pp) + (xmax-xmin)
    end if

    if (particles%coordinates(1,pp) .gt. xmax) then
      particles%coordinates(1,pp) = particles%coordinates(1,pp) - (xmax-xmin)
    end if
  end do  

else if (idpz .eq. myid) then
  do pp = num_local, 1, -1
    if (particles%coordinates(3,pp) .lt. zmin) then
      particles%coordinates(3,pp) = particles%coordinates(3,pp) + (zmax-zmin)
    end if

    if (particles%coordinates(3,pp) .gt. zmax) then
      particles%coordinates(3,pp) = particles%coordinates(3,pp) - (zmax-zmin)
    end if
  end do  

else if (idpy .eq. myid) then
  do pp = num_local, 1, -1
    if (particles%coordinates(2,pp) .lt. ymin) then
      particles%coordinates(2,pp) = particles%coordinates(2,pp) + (ymax-ymin)
    end if

    if (particles%coordinates(2,pp) .gt. ymax) then
      particles%coordinates(2,pp) = particles%coordinates(2,pp) - (ymax-ymin)
    end if
  end do  

else if (idpx .eq. myid) then
  do pp = num_local, 1, -1
    if (particles%coordinates(1,pp) .lt. xmin) then
      particles%coordinates(1,pp) = particles%coordinates(1,pp) + (xmax-xmin)
    end if

    if (particles%coordinates(1,pp) .gt. xmax) then
      particles%coordinates(1,pp) = particles%coordinates(1,pp) - (xmax-xmin)
    end if
  end do  

end if

  return
end subroutine ParticlesTransfer
