recursive subroutine  DumpDump(Nspecies, num_local, max_per_proc, &
     int_particles, real_particles, myid, iteration, attemptno, Nretries)

  use mpi
  use ifport

  implicit none

  integer, intent(in) :: Nspecies, max_per_proc, iteration, myid, Nretries, &
       int_particles(max_per_proc), num_local
  real*8, intent(in) :: real_particles(8,max_per_proc)
  integer, intent(inout) :: attemptno

  integer ierr, ii, jj
  character filename*42, tmpfile*44, pnumber*4

  attemptno = attemptno + 1

! Process 0 writes the iteration number
  if (myid==0) then
     tmpfile = 'dumps/XXdump_iteration.picard.dump'
     filename = 'dumps/dump_iteration.picard.dump'
     open(unit=1,file=tmpfile,status='replace',err=97)
     write (1,fmt='(I7)',err=98) iteration
     close(1,err=99)
     ierr = rename(tmpfile,filename)

  end if

! Everybody writes the particles
  write(pnumber,fmt='(i4.4)') myid
  tmpfile='dumps/XXdump_particles.p'//pnumber//'.picard.dump' 
  filename='dumps/dump_particles.p'//pnumber//'.picard.dump' 
  open(unit=1,form='unformatted',file=tmpfile,status='replace',err=97)
  write (1,err=98) num_local
  write (1,err=98) (int_particles(ii), ii=1, max_per_proc)
  do jj = 1, max_per_proc
     write (1,err=98) (real_particles(ii,jj), ii=1,8)
  end do
  close (1,err=99)

  ! Make sure all processes have written their temporary files 
  ! and then rename them.
  call MPI_Barrier(MPI_comm_world, ierr)
  ierr = rename(tmpfile,filename)

  return

! Error handling section
97 write (*,*) 'DumpDump: error in open statement, myid=', myid
  goto 100
98 write (*,*) 'DumpDump: error in write statement, myid=', myid
  close(1,err=99)
  goto 100
99 write (*,*) 'DumpDump: error in close statement, myid=', myid
100 if (attemptno<=Nretries) then
     call DumpDump(Nspecies, num_local, max_per_proc, &
          int_particles, real_particles, myid, iteration, attemptno, Nretries)
  end if

end subroutine DumpDump

!------------------------------------------------------------------------

recursive subroutine  LoadDump(Nspecies, num_local, max_per_proc, &
     int_particles, real_particles, myid, iteration, attemptno, Nretries)

  use mpi

  implicit none

  integer, intent(in) :: Nspecies, max_per_proc,  myid, Nretries
  real*8, intent(inout) :: real_particles(8,max_per_proc)
  integer, intent(inout) :: num_local, iteration, attemptno, &
       int_particles(max_per_proc)

  integer ierr, ii, jj
  character filename*42, pnumber*4

  attemptno = attemptno + 1

! Read the iteration number
  filename = 'dumps/dump_iteration.picard.dump'
  open(unit=1,file=filename,status='old',err=97)
  read (1,fmt='(I7)',err=98) iteration
  close(1,err=99)


! Read the particles
  write(pnumber,fmt='(i4.4)') myid
  filename='dumps/dump_particles.p'//pnumber//'.picard.dump' 
  open(unit=1,form='unformatted',file=filename,status='old',err=97)
  read (1,err=98) num_local
  read (1,err=98) (int_particles(ii), ii=1, max_per_proc)
  do jj = 1, max_per_proc
     read (1,err=98) (real_particles(ii,jj), ii=1,8)
  end do
  close (1,err=99)

  return

! Error handling section
97 write (*,*) 'LoadDump: error in open statement, myid=', myid
  goto 100
98 write (*,*) 'LoadDump: error in write statement, myid=', myid
  close(1,err=99)
  goto 100
99 write (*,*) 'LoadDump: error in close statement, myid=', myid
100 if (attemptno<=Nretries) then
     call LoadDump(Nspecies, num_local, max_per_proc, &
          int_particles, real_particles, myid, iteration, attemptno, Nretries)
  else
     write (*,*) 'LoadDump: disk access failed ',Nretries,' times, myid=', myid
     stop
  end if

end subroutine LoadDump

!------------------------------------------------------------------------

subroutine DumpParticles( real_particles, int_particles, num_local, &
     max_per_proc, iteration, myid )
  implicit none

  integer, intent(in) :: num_local, max_per_proc, myid
  integer, intent(in) :: int_particles(max_per_proc)
  real*8, intent(inout) :: real_particles(8,max_per_proc)

  integer, intent(in) :: iteration

  ! locals
  
  character :: itername*6, myidname*5, filename*50
  integer ii, jj
  
  real_particles(7,1:num_local) = dble( int_particles(1:num_local) )

  write(itername,fmt='(I6.6)') iteration
  write(myidname,fmt='(I5.5)') myid

  filename = 'outp/datfiles/particles/particles_'// &
       itername//'p'//myidname//'.dat'

  open(unit=1,file=filename)
  do ii = 1, num_local
     write(1,fmt='(7(E11.4E3,a))') (real_particles(jj,ii),' ',jj=1,7)
  end do

  close(1)

  return
end subroutine DumpParticles

!------------------------------------------------------------------------

subroutine DumpFields( U, N, F, E, Nx_local, Ny_local, Nz_local, &
     iteration, Nspecies, dV, myid )
  implicit none

  integer, intent(in) :: Nx_local, Ny_local, Nz_local, Nspecies, myid
  real*8, intent(in) :: E(3,Nx_local+2, Ny_local+2, Nz_local+2)
  real*8, intent(in) :: F(Nspecies,3,Nx_local+2, Ny_local+2, Nz_local+2)
  real*8, intent(in) :: N(Nspecies,Nx_local+2, Ny_local+2, Nz_local+2)
  real*8, intent(in) :: U(Nx_local+2, Ny_local+2, Nz_local+2)
  real*8, intent(in) :: dV

  integer, intent(in) :: iteration

  ! locals
  integer :: ii, jj, kk, ll
  character :: itername*6, myidname*5, filename*50, form*20
 

  write(itername,fmt='(I6.6)') iteration
  write(myidname,fmt='(I5.5)') myid

  ! Write the potential
  filename = 'outp/datfiles/potential/UE_'//itername//'p'//myidname//'.dat'
  open(unit=1,file=filename)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
        write(1,fmt='(E11.4E3)') U(ii,jj,kk)
      end do
    end do
  end do
  close(1)

  ! Write the E-field: Ex Ey Ez
  filename = 'outp/datfiles/Efield/E_'//itername//'p'//myidname//'.dat'
  open(unit=1,file=filename)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
       do ii = 1, Nx_local+2
          write(1,fmt='(3(E11.4E3,a))') (E(ll,ii,jj,kk),' ',ll=1,3)
      end do
    end do
  end do
  close(1)

  ! Density: n1 n2 n3...
  ! We dump density, not particles per cell. Hence the division by dV below.
  write (form,fmt='(a,i5.5,a)') '(', Nspecies, '(E11.4E3,a))'
  filename = 'outp/datfiles/density/n_'//itername//'p'//myidname//'.dat'
  open(unit=1,file=filename)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
        write(1,fmt=form) (N(ll,ii,jj,kk)/dV,' ',ll=1,Nspecies)
      end do
    end do
  end do
  close(1)


  ! Write the fluxes: Fx Fy Fz, one file per species
  if (Nspecies .ge. 1) then
  filename = 'outp/datfiles/flux/F1_'//itername//'p'//myidname//'.dat'
  open(unit=1,file=filename)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
         write(1,fmt='(E11.4E3,a,E11.4E3,a,E11.4E3)') &
              F(1,1,ii,jj,kk), ' ', F(1,2,ii,jj,kk), ' ', F(1,3,ii,jj,kk)
      end do
    end do
  end do
  close(1)

  if (Nspecies .ge. 2) then
  filename = 'outp/datfiles/flux/F2_'//itername//'p'//myidname//'.dat'
  open(unit=1,file=filename)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
         write(1,fmt='(E11.4E3,a,E11.4E3,a,E11.4E3)') &
              F(2,1,ii,jj,kk), ' ', F(2,2,ii,jj,kk), ' ', F(2,3,ii,jj,kk)
      end do
    end do
  end do
  close(1)

  if (Nspecies .ge. 3) then
  filename = 'outp/datfiles/flux/F3_'//itername//'p'//myidname//'.dat'
  open(unit=1,file=filename)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
         write(1,fmt='(E11.4E3,a,E11.4E3,a,E11.4E3)') &
              F(3,1,ii,jj,kk), ' ', F(3,2,ii,jj,kk), ' ', F(3,3,ii,jj,kk)
      end do
    end do
  end do
  close(1)

  if (Nspecies .ge. 4) then
  filename = 'outp/datfiles/flux/F4_'//itername//'p'//myidname//'.dat'
  open(unit=1,file=filename)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
         write(1,fmt='(E11.4E3,a,E11.4E3,a,E11.4E3)') &
              F(4,1,ii,jj,kk), ' ', F(4,2,ii,jj,kk), ' ', F(4,3,ii,jj,kk)
      end do
    end do
  end do
  close(1)

  if (Nspecies .ge. 5) then
  filename = 'outp/datfiles/flux/F5_'//itername//'p'//myidname//'.dat'
  open(unit=1,file=filename)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
         write(1,fmt='(E11.4E3,a,E11.4E3,a,E11.4E3)') &
              F(5,1,ii,jj,kk), ' ', F(5,2,ii,jj,kk), ' ', F(5,3,ii,jj,kk)
      end do
    end do
  end do
  close(1)

  if (Nspecies .ge. 6) then
  filename = 'outp/datfiles/flux/F6_'//itername//'p'//myidname//'.dat'
  open(unit=1,file=filename)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
         write(1,fmt='(E11.4E3,a,E11.4E3,a,E11.4E3)') &
              F(6,1,ii,jj,kk), ' ', F(6,2,ii,jj,kk), ' ', F(6,3,ii,jj,kk)
      end do
    end do
  end do
  close(1)

  if (Nspecies .ge. 7) then
  filename = 'outp/datfiles/flux/F7_'//itername//'p'//myidname//'.dat'
  open(unit=1,file=filename)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
         write(1,fmt='(E11.4E3,a,E11.4E3,a,E11.4E3)') &
              F(7,1,ii,jj,kk), ' ', F(7,2,ii,jj,kk), ' ', F(7,3,ii,jj,kk)
      end do
    end do
  end do
  close(1)

  if (Nspecies .ge. 8) then
  filename = 'outp/datfiles/flux/F8_'//itername//'p'//myidname//'.dat'
  open(unit=1,file=filename)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
         write(1,fmt='(E11.4E3,a,E11.4E3,a,E11.4E3)') &
              F(8,1,ii,jj,kk), ' ', F(8,2,ii,jj,kk), ' ', F(8,3,ii,jj,kk)
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


  return

end subroutine DumpFields

!------------------------------------------------------------------------

! Here I can put my new dump probe function. HG 2018-06-26
