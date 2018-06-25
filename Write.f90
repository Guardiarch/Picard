
subroutine DumpParticles( real_particles, int_particles, num_local, max_per_proc, writeParticles_iteration, myid )
  implicit none

  integer, intent(in) :: num_local, max_per_proc, myid
  integer, intent(in) :: int_particles(max_per_proc)
  real*8, intent(inout) :: real_particles(8,max_per_proc)

  integer, intent(inout) :: writeParticles_iteration

  ! locals
  
  character :: itername*5, myidname*5, filename*36

  real_particles(7,1:num_local) = dble( int_particles(1:num_local) )

  write(itername,fmt='(I5.5)') writeParticles_iteration
  write(myidname,fmt='(I5.5)') myid

  filename = 'save/particles/particles_'//itername//'_'//myidname

  open(unit=1,file=filename)

  write(1,fmt='(7E10.4)') real_particles(1:7,1:num_local)

  close(1)

  writeParticles_iteration = writeParticles_iteration + 1

  return
end subroutine DumpParticles


subroutine DumpFields( U, N, F, E, Nx_local, Ny_local, Nz_local, writeFields_iteration, species, myid )
  implicit none

  integer, intent(in) :: Nx_local, Ny_local, Nz_local, species, myid
  real*8, intent(in) :: E(3,Nx_local+2, Ny_local+2, Nz_local+2)
  real*8, intent(in) :: F(species,3,Nx_local+2, Ny_local+2, Nz_local+2)
  real*8, intent(in) :: N(species,Nx_local+2, Ny_local+2, Nz_local+2)
  real*8, intent(in) :: U(Nx_local+2, Ny_local+2, Nz_local+2)

  integer, intent(inout) :: writeFields_iteration

  ! locals
  integer :: ii, jj, kk
  character :: itername*5, myidname*5, filename*26


  write(itername,fmt='(I5.5)') writeFields_iteration
  write(myidname,fmt='(I5.5)') myid

  filename = 'save/fields/UE_'//itername//'_'//myidname
  open(unit=1,file=filename)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
        write(1,fmt='(6E10.4)') U(ii,jj,kk)
      end do
    end do
  end do
  close(1)

  filename = 'save/fields/Ex_'//itername//'_'//myidname
  open(unit=1,file=filename)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
        write(1,fmt='(6E10.4)') E(1,ii,jj,kk)
      end do
    end do
  end do
  close(1)

  filename = 'save/fields/Ey_'//itername//'_'//myidname
  open(unit=1,file=filename)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
        write(1,fmt='(6E10.4)') E(2,ii,jj,kk)
      end do
    end do
  end do
  close(1)

  filename = 'save/fields/Ez_'//itername//'_'//myidname
  open(unit=1,file=filename)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
        write(1,fmt='(6E10.4)') E(3,ii,jj,kk)
      end do
    end do
  end do
  close(1)


  if (species .ge. 1) then
  filename = 'save/fields/N1_'//itername//'_'//myidname
  open(unit=1,file=filename)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
        write(1,fmt='(6E10.4)') N(1,ii,jj,kk)
      end do
    end do
  end do
  close(1)

  if (species .ge. 2) then
  filename = 'save/fields/N2_'//itername//'_'//myidname
  open(unit=1,file=filename)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
        write(1,fmt='(6E10.4)') N(2,ii,jj,kk)
      end do
    end do
  end do
  close(1)

  if (species .ge. 3) then
  filename = 'save/fields/N3_'//itername//'_'//myidname
  open(unit=1,file=filename)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
        write(1,fmt='(6E10.4)') N(3,ii,jj,kk)
      end do
    end do
  end do
  close(1)

  if (species .ge. 4) then
  filename = 'save/fields/N4_'//itername//'_'//myidname
  open(unit=1,file=filename)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
        write(1,fmt='(6E10.4)') N(4,ii,jj,kk)
      end do
    end do
  end do
  close(1)

  if (species .ge. 5) then
  filename = 'save/fields/N5_'//itername//'_'//myidname
  open(unit=1,file=filename)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
        write(1,fmt='(6E10.4)') N(5,ii,jj,kk)
      end do
    end do
  end do
  close(1)

  if (species .ge. 6) then
  filename = 'save/fields/N6_'//itername//'_'//myidname
  open(unit=1,file=filename)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
        write(1,fmt='(6E10.4)') N(6,ii,jj,kk)
      end do
    end do
  end do
  close(1)

  if (species .ge. 7) then
  filename = 'save/fields/N7_'//itername//'_'//myidname
  open(unit=1,file=filename)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
        write(1,fmt='(6E10.4)') N(7,ii,jj,kk)
      end do
    end do
  end do
  close(1)

  if (species .ge. 8) then
  filename = 'save/fields/N8_'//itername//'_'//myidname
  open(unit=1,file=filename)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
        write(1,fmt='(6E10.4)') N(8,ii,jj,kk)
      end do
    end do
  end do
  close(1)

  end if
  end if
  end if
  end if
  end if
  end if
  end if
  end if

  if (species .ge. 1) then
  filename = 'save/fields/X1_'//itername//'_'//myidname
  open(unit=1,file=filename)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
        write(1,fmt='(6E10.4)') F(1,1,ii,jj,kk)
      end do
    end do
  end do
  close(1)

  if (species .ge. 2) then
  filename = 'save/fields/X2_'//itername//'_'//myidname
  open(unit=1,file=filename)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
        write(1,fmt='(6E10.4)') F(2,1,ii,jj,kk)
      end do
    end do
  end do
  close(1)

  if (species .ge. 3) then
  filename = 'save/fields/X3_'//itername//'_'//myidname
  open(unit=1,file=filename)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
        write(1,fmt='(6E10.4)') F(3,1,ii,jj,kk)
      end do
    end do
  end do
  close(1)

  if (species .ge. 4) then
  filename = 'save/fields/X4_'//itername//'_'//myidname
  open(unit=1,file=filename)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
        write(1,fmt='(6E10.4)') F(4,1,ii,jj,kk)
      end do
    end do
  end do
  close(1)

  if (species .ge. 5) then
  filename = 'save/fields/X5_'//itername//'_'//myidname
  open(unit=1,file=filename)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
        write(1,fmt='(6E10.4)') F(5,1,ii,jj,kk)
      end do
    end do
  end do
  close(1)

  if (species .ge. 6) then
  filename = 'save/fields/X6_'//itername//'_'//myidname
  open(unit=1,file=filename)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
        write(1,fmt='(6E10.4)') F(6,1,ii,jj,kk)
      end do
    end do
  end do
  close(1)

  if (species .ge. 7) then
  filename = 'save/fields/X7_'//itername//'_'//myidname
  open(unit=1,file=filename)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
        write(1,fmt='(6E10.4)') F(7,1,ii,jj,kk)
      end do
    end do
  end do
  close(1)

  if (species .ge. 8) then
  filename = 'save/fields/X8_'//itername//'_'//myidname
  open(unit=1,file=filename)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
        write(1,fmt='(6E10.4)') F(8,1,ii,jj,kk)
      end do
    end do
  end do
  close(1)

  end if
  end if
  end if
  end if
  end if
  end if
  end if
  end if

  if (species .ge. 1) then
  filename = 'save/fields/Y1_'//itername//'_'//myidname
  open(unit=1,file=filename)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
        write(1,fmt='(6E10.4)') F(1,2,ii,jj,kk)
      end do
    end do
  end do
  close(1)

  if (species .ge. 2) then
  filename = 'save/fields/Y2_'//itername//'_'//myidname
  open(unit=1,file=filename)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
        write(1,fmt='(6E10.4)') F(2,2,ii,jj,kk)
      end do
    end do
  end do
  close(1)

  if (species .ge. 3) then
  filename = 'save/fields/Y3_'//itername//'_'//myidname
  open(unit=1,file=filename)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
        write(1,fmt='(6E10.4)') F(3,2,ii,jj,kk)
      end do
    end do
  end do
  close(1)

  if (species .ge. 4) then
  filename = 'save/fields/Y4_'//itername//'_'//myidname
  open(unit=1,file=filename)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
        write(1,fmt='(6E10.4)') F(4,2,ii,jj,kk)
      end do
    end do
  end do
  close(1)

  if (species .ge. 5) then
  filename = 'save/fields/Y5_'//itername//'_'//myidname
  open(unit=1,file=filename)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
        write(1,fmt='(6E10.4)') F(5,2,ii,jj,kk)
      end do
    end do
  end do
  close(1)

  if (species .ge. 6) then
  filename = 'save/fields/Y6_'//itername//'_'//myidname
  open(unit=1,file=filename)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
        write(1,fmt='(6E10.4)') F(6,2,ii,jj,kk)
      end do
    end do
  end do
  close(1)

  if (species .ge. 7) then
  filename = 'save/fields/Y7_'//itername//'_'//myidname
  open(unit=1,file=filename)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
        write(1,fmt='(6E10.4)') F(7,2,ii,jj,kk)
      end do
    end do
  end do
  close(1)

  if (species .ge. 8) then
  filename = 'save/fields/Y8_'//itername//'_'//myidname
  open(unit=1,file=filename)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
        write(1,fmt='(6E10.4)') F(8,2,ii,jj,kk)
      end do
    end do
  end do
  close(1)

  end if
  end if
  end if
  end if
  end if
  end if
  end if
  end if


  if (species .ge. 1) then
  filename = 'save/fields/Z1_'//itername//'_'//myidname
  open(unit=1,file=filename)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
        write(1,fmt='(6E10.4)') F(1,3,ii,jj,kk)
      end do
    end do
  end do
  close(1)

  if (species .ge. 2) then
  filename = 'save/fields/Z2_'//itername//'_'//myidname
  open(unit=1,file=filename)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
        write(1,fmt='(6E10.4)') F(2,3,ii,jj,kk)
      end do
    end do
  end do
  close(1)

  if (species .ge. 3) then
  filename = 'save/fields/Z3_'//itername//'_'//myidname
  open(unit=1,file=filename)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
        write(1,fmt='(6E10.4)') F(3,3,ii,jj,kk)
      end do
    end do
  end do
  close(1)

  if (species .ge. 4) then
  filename = 'save/fields/Z4_'//itername//'_'//myidname
  open(unit=1,file=filename)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
        write(1,fmt='(6E10.4)') F(4,3,ii,jj,kk)
      end do
    end do
  end do
  close(1)

  if (species .ge. 5) then
  filename = 'save/fields/Z5_'//itername//'_'//myidname
  open(unit=1,file=filename)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
        write(1,fmt='(6E10.4)') F(5,3,ii,jj,kk)
      end do
    end do
  end do
  close(1)

  if (species .ge. 6) then
  filename = 'save/fields/Z6_'//itername//'_'//myidname
  open(unit=1,file=filename)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
        write(1,fmt='(6E10.4)') F(6,3,ii,jj,kk)
      end do
    end do
  end do
  close(1)

  if (species .ge. 7) then
  filename = 'save/fields/Z7_'//itername//'_'//myidname
  open(unit=1,file=filename)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
        write(1,fmt='(6E10.4)') F(7,3,ii,jj,kk)
      end do
    end do
  end do
  close(1)

  if (species .ge. 8) then
  filename = 'save/fields/Z8_'//itername//'_'//myidname
  open(unit=1,file=filename)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
        write(1,fmt='(6E10.4)') F(8,3,ii,jj,kk)
      end do
    end do
  end do
  close(1)

  end if
  end if
  end if
  end if
  end if
  end if
  end if
  end if


  writeFields_iteration = writeFields_iteration + 1
  return

end subroutine DumpFields

