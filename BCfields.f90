subroutine BCfields(E,Nx_local,Ny_local,Nz_local, &
     xflow,yflow,zflow,iprocs,jprocs,kprocs,myid)

  implicit none

  logical, intent(in) :: xflow,yflow,zflow
  integer, intent(in) :: Nx_local,Ny_local,Nz_local,iprocs,jprocs,kprocs,myid
  real*8, intent(inout) :: E(3,Nx_local+2,Ny_local+2,Nz_local+2)

  real*8, parameter :: mu0  = 1.25663706143591729538505735331180d-6  ! Magnetic constant [ SI ]

! locals
  integer :: idpx, idnx, idpy, idny, idpz, idnz
  real*8 :: Sx(Ny_local+2,Nz_local+2)


  idpx = myid + (mod(myid+iprocs+1,                           iprocs) - mod(myid+iprocs, iprocs))
  idnx = myid + (mod(myid+iprocs-1,                           iprocs) - mod(myid+iprocs, iprocs))
  idpy = myid + (mod(myid+iprocs*jprocs+iprocs,               iprocs*jprocs) - mod(myid+iprocs*jprocs, iprocs*jprocs))
  idny = myid + (mod(myid+iprocs*jprocs-iprocs,               iprocs*jprocs) - mod(myid+iprocs*jprocs, iprocs*jprocs))
  idpz = myid + (mod(myid+iprocs*jprocs*kprocs+iprocs*jprocs, iprocs*jprocs*kprocs) &
               - mod(myid+iprocs*jprocs*kprocs,               iprocs*jprocs*kprocs))
  idnz = myid + (mod(myid+iprocs*jprocs*kprocs-iprocs*jprocs, iprocs*jprocs*kprocs) &
               - mod(myid+iprocs*jprocs*kprocs,               iprocs*jprocs*kprocs))

! if I have the boundary, make it into the boundary field.

  if (zflow) then

    if (idpz - myid .lt. 1) then
!      E(:,:,:,Nz_local+1:Nz_local+2) = 0.0d0
      E(:,:,:,Nz_local+2) = 0.0d0
    end if
    
    if (myid - idnz .lt. 1) then
!      E(:,:,:,1:2)                   = 0.0d0
      E(:,:,:,1)                   = 0.0d0
    end if

  end if

  if (yflow) then

    if (idpy - myid .lt. 1) then
!      E(:,:,Ny_local+1:Ny_local+2,:) = 0.0d0
      E(:,:,Ny_local+2,:) = 0.0d0
    end if
    
    if (myid - idny .lt. 1) then
!      E(:,:,1:2,:)                   = 0.0d0
      E(:,:,1,:)                   = 0.0d0
    end if

  end if

  if (xflow) then

! INFLOW
    if (idpx - myid .lt. 1) then
!      E(:,Nx_local+1:Nx_local+2,:,:) = 0.0d0
      E(:,Nx_local+2,:,:) = 0.0d0
    end if

! OUTFLOW
    if (myid - idnx .lt. 1) then
       E(:,1,:,:)                     = 0.0d0
    end if
    
  end if

  return

end subroutine BCfields
