subroutine BCpotentials(U,C,D,xmin,xmax,ymin,ymax,zmin,zmax, &
     Nx_local,Ny_local,Nz_local,Nx,Ny,Nz, &
     xflow,yflow,zflow,Nspecies,iprocs,jprocs,kprocs,myid)
  use mpi
  implicit none


  logical, intent(in) :: xflow, yflow, zflow
  integer, intent(in) :: Nx_local, Ny_local, Nz_local, Nx, Ny, Nz, Nspecies, iprocs, jprocs, kprocs, myid
  real*8, intent(in) :: xmin, xmax, ymin, ymax, zmin, zmax
  real*8, intent(in) :: D(2*Nx+1,2*Ny+1,2*Nz+1)
  real*8, intent(in) :: C(Nx_local+2, Ny_local+2, Nz_local+2)
  real*8, intent(inout) :: U(Nx_local+2, Ny_local+2, Nz_local+2)

! locals
  integer :: ierr, idpx, idnx, idpy, idny, idpz, idnz, ii, jj, kk, bi, bj, bk, di, dj, dk
  real*8 :: Ux0_loc(Ny,Nz), Ux1_loc(Ny,Nz), Uy0_loc(Nx,Nz), Uy1_loc(Nx,Nz), Uz0_loc(Nx,Ny), Uz1_loc(Nx,Ny)
  real*8 :: Ux0_glo(Ny,Nz), Ux1_glo(Ny,Nz), Uy0_glo(Nx,Nz), Uy1_glo(Nx,Nz), Uz0_glo(Nx,Ny), Uz1_glo(Nx,Ny)

  real*8, parameter :: eps0 = 8.854187817d-12  ! Electric constant [ SI ]
  real*8, parameter :: tau  = 6.2831853071795864769252867665590d0  !  2*pi

  idpx = myid + (mod(myid+iprocs+1,                           iprocs) - mod(myid+iprocs, iprocs))
  idnx = myid + (mod(myid+iprocs-1,                           iprocs) - mod(myid+iprocs, iprocs))
  idpy = myid + (mod(myid+iprocs*jprocs+iprocs,               iprocs*jprocs) - mod(myid+iprocs*jprocs, iprocs*jprocs))
  idny = myid + (mod(myid+iprocs*jprocs-iprocs,               iprocs*jprocs) - mod(myid+iprocs*jprocs, iprocs*jprocs))
  idpz = myid + (mod(myid+iprocs*jprocs*kprocs+iprocs*jprocs, iprocs*jprocs*kprocs) &
               - mod(myid+iprocs*jprocs*kprocs,               iprocs*jprocs*kprocs))
  idnz = myid + (mod(myid+iprocs*jprocs*kprocs-iprocs*jprocs, iprocs*jprocs*kprocs) &
               - mod(myid+iprocs*jprocs*kprocs,               iprocs*jprocs*kprocs))

! Block number
  bi = mod(myid+iprocs, iprocs)
  bj = mod(myid/iprocs+jprocs, jprocs)
  bk = mod(myid/(iprocs*jprocs)+kprocs, kprocs)

  Ux0_loc = 0.0d0
  Ux1_loc = 0.0d0
  Uy0_loc = 0.0d0
  Uy1_loc = 0.0d0
  Uz0_loc = 0.0d0
  Uz1_loc = 0.0d0
  Ux0_glo = 0.0d0
  Ux1_glo = 0.0d0
  Uy0_glo = 0.0d0
  Uy1_glo = 0.0d0
  Uz0_glo = 0.0d0
  Uz1_glo = 0.0d0


  if (xflow) then

    do kk = 1, Nz
      do jj = 1, Ny
        ! Distance to block start from boundary in terms of cells
        dj = Ny+1+(1-jj+Ny_local*bj)
        dk = Nz+1+(1-kk+Nz_local*bk)

        di = Nx+1+(1+Nx_local*bi)
        Ux0_loc(jj,kk) = sum( C(2:Nx_local+1,2:Ny_local+1,2:Nz_local+1)*D(di:di+Nx_local-1,dj:dj+Ny_local-1,dk:dk+Nz_local-1) )

        di = Nx+1-(Nx_local*(iprocs-bi))
        Ux1_loc(jj,kk) = sum( C(2:Nx_local+1,2:Ny_local+1,2:Nz_local+1)*D(di:di+Nx_local-1,dj:dj+Ny_local-1,dk:dk+Nz_local-1) )
      end do
    end do

    call MPI_allreduce( Ux0_loc, Ux0_glo, Ny*Nz, MPI_DOUBLE_PRECISION, &
         MPI_SUM, MPI_comm_world, ierr )
    if (ierr>0) then
       write (*,*) 'BCpotentials, 1, ierr=',ierr,' myid=',myid
    end if
    call MPI_allreduce( Ux1_loc, Ux1_glo, Ny*Nz, MPI_DOUBLE_PRECISION, &
         MPI_SUM, MPI_comm_world, ierr )
    if (ierr>0) then
       write (*,*) 'BCpotentials, 2, ierr=',ierr,' myid=',myid
    end if

    if (myid - idnx .le. 0) then
      U(1,          2:Ny_local+1, 2:Nz_local+1) = &
      Ux0_glo(1+bj*Ny_local:Ny_local*(bj+1),1+bk*Nz_local:Nz_local*(bk+1))*(2.0d0*tau*eps0)**(-1)
    end if
    if (idpx - myid .le. 0) then
      U(Nx_local+2, 2:Ny_local+1, 2:Nz_local+1) = &
      Ux1_glo(1+bj*Ny_local:Ny_local*(bj+1),1+bk*Nz_local:Nz_local*(bk+1))*(2.0d0*tau*eps0)**(-1)
    end if

  end if


  if (yflow) then

    do kk = 1, Nz
      do ii = 1, Nx
        ! Distance to block start from boundary in terms of cells
        di = Nx+1+(1-ii+Nx_local*bi)
        dk = Nz+1+(1-kk+Nz_local*bk)

        dj = Ny+1+(1+Ny_local*bj)
        Uy0_loc(ii,kk) = sum( C(2:Nx_local+1,2:Ny_local+1,2:Nz_local+1)*D(di:di+Nx_local-1,dj:dj+Ny_local-1,dk:dk+Nz_local-1) )

        dj = Ny+1-(Ny_local*(jprocs-bj))
        Uy1_loc(ii,kk) = sum( C(2:Nx_local+1,2:Ny_local+1,2:Nz_local+1)*D(di:di+Nx_local-1,dj:dj+Ny_local-1,dk:dk+Nz_local-1) )
      end do
    end do

    call MPI_allreduce( Uy0_loc, Uy0_glo, Nx*Nz, MPI_DOUBLE_PRECISION, &
         MPI_SUM, MPI_comm_world, ierr )
    if (ierr>0) then
       write (*,*) 'BCpotentials, 3, ierr=',ierr,' myid=',myid
    end if
    call MPI_allreduce( Uy1_loc, Uy1_glo, Nx*Nz, MPI_DOUBLE_PRECISION, &
         MPI_SUM, MPI_comm_world, ierr )
    if (ierr>0) then
       write (*,*) 'BCpotentials, 4, ierr=',ierr,' myid=',myid
    end if

    if (myid - idny .le. 0) then
      U(2:Nx_local+1, 1,          2:Nz_local+1) = &
      Uy0_glo(1+bi*Nx_local:Nx_local*(bi+1),1+bk*Nz_local:Nz_local*(bk+1))*(2.0d0*tau*eps0)**(-1)
    end if
    if (idpy - myid .le. 0) then
      U(2:Nx_local+1, Ny_local+2, 2:Nz_local+1) = &
      Uy1_glo(1+bi*Nx_local:Nx_local*(bi+1),1+bk*Nz_local:Nz_local*(bk+1))*(2.0d0*tau*eps0)**(-1)
    end if

  end if


  if (zflow) then

    do jj = 1, Ny
      do ii = 1, Nx
        ! Distance to block start from boundary in terms of cells
        di = Nx+1+(1-ii+Nx_local*bi)
        dj = Ny+1+(1-jj+Ny_local*bj)

        dk = Nz+1+(1+Nz_local*bk)
        Uz0_loc(ii,jj) = sum( C(2:Nx_local+1,2:Ny_local+1,2:Nz_local+1)*D(di:di+Nx_local-1,dj:dj+Ny_local-1,dk:dk+Nz_local-1) )

        dk = Nz+1-(Nz_local*(kprocs-bk))
        Uz1_loc(ii,jj) = sum( C(2:Nx_local+1,2:Ny_local+1,2:Nz_local+1)*D(di:di+Nx_local-1,dj:dj+Ny_local-1,dk:dk+Nz_local-1) )
      end do
    end do

    call MPI_allreduce( Uz0_loc, Uz0_glo, Nx*Ny, MPI_DOUBLE_PRECISION, &
         MPI_SUM, MPI_comm_world, ierr )
    if (ierr>0) then
       write (*,*) 'BCpotentials, 5, ierr=',ierr,' myid=',myid
    end if
    call MPI_allreduce( Uz1_loc, Uz1_glo, Nx*Ny, MPI_DOUBLE_PRECISION, &
         MPI_SUM, MPI_comm_world, ierr )
    if (ierr>0) then
       write (*,*) 'BCpotentials, 6, ierr=',ierr,' myid=',myid
    end if

    if (myid - idnz .le. 0) then
      U(2:Nx_local+1, 2:Ny_local+1, 1         ) = &
      Uz0_glo(1+bi*Nx_local:Nx_local*(bi+1),1+bj*Ny_local:Ny_local*(bj+1))*(2.0d0*tau*eps0)**(-1)
    end if
    if (idpz - myid .le. 0) then
      U(2:Nx_local+1, 2:Ny_local+1, Nz_local+2) = &
      Uz1_glo(1+bi*Nx_local:Nx_local*(bi+1),1+bj*Ny_local:Ny_local*(bj+1))*(2.0d0*tau*eps0)**(-1)
    end if

  end if


  return
end subroutine BCpotentials
