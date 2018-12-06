recursive subroutine  DumpDump(Nspecies, num_local, max_per_proc, &
     particles, myid, iteration, U, Nx_local, Ny_local, Nz_local, &
     attemptno, Nretries)

  use SpecificTypes
  use mpi
  use ifport

  implicit none

  integer, intent(in) :: Nspecies, max_per_proc, iteration, myid, Nretries, &
       num_local, Nx_local, Ny_local, Nz_local
  integer, intent(inout) :: attemptno
  real*8, intent(in) :: U(Nx_local+2,Ny_local+2,Nz_local+2)
  type(particlearrays) particles

  integer ierr, ii, jj
  character filename*43, tmpfile*45, pnumber*5

  attemptno = attemptno + 1

! Everybody writes the particles
  write(pnumber,fmt='(i5.5)') myid
  tmpfile='dumps/XXdump_particles.p'//pnumber//'.picard.dump' 
  filename='dumps/dump_particles.p'//pnumber//'.picard.dump' 
  open(unit=1,form='unformatted',file=tmpfile,status='replace',err=97)
  write (1,err=98) num_local
  write (1,err=98) (particles%species(ii), ii=1, num_local)
  do jj = 1, num_local
     write (1,err=98) (particles%coordinates(ii,jj), ii=1,6)
  end do
  ! and the potential
  write (1,err=98) U
  close (1,err=99)

  ! Make sure all processes have written their temporary files 
  ! and then rename them.
  call MPI_Barrier(MPI_comm_world, ierr)
  if (ierr>0) then
     write (*,*) 'DumpDump, 1, ierr=',ierr,' myid=',myid
  end if
  ierr = rename(tmpfile,filename)
  if (ierr>0) then
     write (*,*) 'DumpDump, 2, ierr=',ierr,' myid=',myid
  end if

! Process 0 writes the iteration number
  if (myid==0) then
     tmpfile = 'dumps/XXdump_iteration.picard.dump'
     filename = 'dumps/dump_iteration.picard.dump'
     open(unit=1,file=tmpfile,status='replace',err=97)
     write (1,fmt='(I7)',err=98) iteration
     close(1,err=99)
     ierr = rename(tmpfile,filename)
     if (ierr>0) then
        write (*,*) 'DumpDump, 3, ierr=',ierr,' myid=',myid
     end if
  end if

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
          particles, myid, iteration, U, Nx_local, Ny_local, Nz_local, &
          attemptno, Nretries)
  end if

end subroutine DumpDump

!------------------------------------------------------------------------

recursive subroutine  LoadDump(Nspecies, num_local, max_per_proc, &
     particles, myid, iteration, U, Nx_local, Ny_local, Nz_local, &
     attemptno, Nretries)

  use SpecificTypes

  implicit none

  integer, intent(in) :: Nspecies, max_per_proc,  myid, Nretries
  integer, intent(inout) :: num_local, Nx_local, Ny_local, Nz_local, &
       iteration, attemptno
  real*8, intent(inout) :: U(Nx_local+2,Ny_local+2,Nz_local+2)
  type(particlearrays) particles

  integer ierr, ii, jj
  character filename*43, pnumber*5

  U = 0.0d0

  attemptno = attemptno + 1

  
! Read the iteration number
  filename = 'dumps/dump_iteration.picard.dump'
  open(unit=1,file=filename,status='old',err=97)
  read (1,fmt='(I7)',err=98) iteration
  close(1,err=99)


! Read the particles
  write(pnumber,fmt='(i5.5)') myid
  filename='dumps/dump_particles.p'//pnumber//'.picard.dump' 
  open(unit=1,form='unformatted',file=filename,status='old',err=97)
  read (1,err=98) num_local
  read (1,err=98) (particles%species(ii), ii=1, num_local)
  do jj = 1, num_local
     read (1,err=98) (particles%coordinates(ii,jj), ii=1,6)
  end do
  ! and the potential
  read (1,err=98) U  
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
          particles, myid, iteration, U, Nx_local, Ny_local, Nz_local, &
          attemptno, Nretries)
  else
     write (*,*) 'LoadDump: disk access failed ',Nretries,' times, myid=', myid
     stop
  end if

end subroutine LoadDump

!------------------------------------------------------------------------

subroutine DumpParticles( particles, num_local, iteration, myid )

  use SpecificTypes
  implicit none

  integer, intent(in) :: num_local, myid
  type(particlearrays) particles

  integer, intent(in) :: iteration

  ! locals
  
  character :: itername*6, myidname*5, filename*57
  integer ii, jj
  
  write(itername,fmt='(I6.6)') iteration
  write(myidname,fmt='(I5.5)') myid

  filename = 'outp/datfiles/particles/particles_'// &
       itername//'p'//myidname//'.picard.dat'

  open(unit=1,file=filename)
  do ii = 1, num_local
     write(1,fmt='(6(E11.4E3,a),I3.3)') &
          (particles%coordinates(jj,ii),' ',jj=1,6), particles%species(ii)
  end do

  close(1)

  return
end subroutine DumpParticles

!------------------------------------------------------------------------

recursive subroutine DumpEfield( E, Nx_local, Ny_local, Nz_local, &
     iteration, Nspecies, myid, attemptno, Nretries)
  implicit none

  integer, intent(in) :: Nx_local, Ny_local, Nz_local, Nspecies, myid
  real*8, intent(in) :: E(3,Nx_local+2, Ny_local+2, Nz_local+2)

  integer, intent(in) :: iteration, Nretries
  integer attemptno

  ! locals
  integer :: ii, jj, kk, ll
  character :: itername*7, myidname*5, filename*58, form*20, snumber*2
 
  attemptno = attemptno + 1
  
  write(itername,fmt='(I7.7)') iteration
  write(myidname,fmt='(I5.5)') myid

  ! Write the E-field: Ex Ey Ez
  filename = 'outp/datfiles/Efield/E_'//itername//'p'//myidname//'.picard.dat'
  open(unit=1,file=filename,err=97)
  do kk = 1, Nz_local+2
     do jj = 1, Ny_local+2
        do ii = 1, Nx_local+2
           write(1,fmt='(3(E11.4E3,a))',err=98) (E(ll,ii,jj,kk),' ',ll=1,3)
        end do
     end do
  end do
 
  close(1,err=99)

  return

  ! Error handling section
97 write (*,*) 'DumpEfield: error in open statement, file ', filename
  goto 100
98 write (*,*) 'DumpEfield: error in write statement, file ', filename
  close(1,err=99)
  goto 100
99 write (*,*) 'DumpEfield: error in close statement, file ', filename
100 if (attemptno<=Nretries) then
     call DumpEfield( E, Nx_local, Ny_local, Nz_local, &
          iteration, Nspecies, myid, attemptno, Nretries)
  end if

end subroutine DumpEfield

!------------------------------------------------------------------------
recursive subroutine DumpPotential( U, Nx_local, Ny_local, Nz_local, &
     iteration, Nspecies, myid, attemptno, Nretries)
  implicit none

  integer, intent(in) :: Nx_local, Ny_local, Nz_local, Nspecies, myid
  real*8, intent(in) :: U(Nx_local+2, Ny_local+2, Nz_local+2)

  integer, intent(in) :: iteration, Nretries
  integer attemptno

  ! locals
  integer :: ii, jj, kk, ll
  character :: itername*7, myidname*5, filename*58, form*20, snumber*2

  attemptno = attemptno +1

  write(itername,fmt='(I7.7)') iteration
  write(myidname,fmt='(I5.5)') myid

  ! Write the potential
  filename = 'outp/datfiles/potential/UE_'//itername//'p'//myidname// &
       '.picard.dat'
  open(unit=1,file=filename,err=97)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
        write(1,fmt='(E11.4E3)',err=98) U(ii,jj,kk)
      end do
    end do
 end do
 
  close(1,err=99)

  return

! Error handling section
97 write (*,*) 'DumpPotential: error in open statement, file ', filename
  goto 100
98 write (*,*) 'DumpPotential: error in write statement, file ', filename
  close(1,err=99)
  goto 100
99 write (*,*) 'DumpPotential: error in close statement, file ', filename
100 if (attemptno<=Nretries) then
     call DumpPotential( U, Nx_local, Ny_local, Nz_local, &
          iteration, Nspecies, myid, attemptno, Nretries)
  end if

end subroutine DumpPotential

!------------------------------------------------------------------------
recursive subroutine DumpDensity( N, Nx_local, Ny_local, Nz_local, &
     iteration, Nspecies, dV, myid, attemptno, Nretries )
  implicit none

  integer, intent(in) :: Nx_local, Ny_local, Nz_local, Nspecies, myid
  real*8, intent(in) :: N(Nspecies,Nx_local+2, Ny_local+2, Nz_local+2)
  real*8, intent(in) :: dV

  integer, intent(in) :: iteration, Nretries
  integer attemptno
  
  ! locals
  integer :: ii, jj, kk, ll
  character :: itername*7, myidname*5, filename*58, form*20, snumber*2
 
  attemptno = attemptno + 1

  write(itername,fmt='(I7.7)') iteration
  write(myidname,fmt='(I5.5)') myid


  ! Density: n1 n2 n3...
  ! We dump density, not particles per cell. Hence the division by dV below.
  write (form,fmt='(a,i5.5,a)') '(', Nspecies, '(E11.4E3,a))'
  filename = 'outp/datfiles/density/n_'//itername//'p'//myidname// &
       '.picard.dat'
  open(unit=1,file=filename,err=97)
  do kk = 1, Nz_local+2
    do jj = 1, Ny_local+2
      do ii = 1, Nx_local+2
        write(1,fmt=form,err=98) (N(ll,ii,jj,kk)/dV,' ',ll=1,Nspecies)
      end do
    end do
  end do
  close(1,err=99)

  return

! Error handling section
97 write (*,*) 'DumpDensity: error in open statement, file ', filename
  goto 100
98 write (*,*) 'DumpDensity: error in write statement, file ', filename
  close(1,err=99)
  goto 100
99 write (*,*) 'DumpDensity: error in close statement, file ', filename
100 if (attemptno<=Nretries) then
     call DumpDensity( N, Nx_local, Ny_local, Nz_local, &
          iteration, Nspecies, dV, myid, attemptno, Nretries )
  end if

end subroutine DumpDensity

!------------------------------------------------------------------------
recursive subroutine DumpFlux( F, Nx_local, Ny_local, Nz_local, &
     iteration, Nspecies, myid, attemptno, Nretries)
  implicit none

  integer, intent(in) :: Nx_local, Ny_local, Nz_local, Nspecies, myid
  real*8, intent(in) :: F(Nspecies,3,Nx_local+2, Ny_local+2, Nz_local+2)

  integer, intent(in) :: iteration
  integer attemptno, Nretries
  
  ! locals
  integer :: ii, jj, kk, ll
  character :: itername*7, myidname*5, filename*58, form*20, snumber*2
 
  attemptno = attemptno + 1
  
  write(itername,fmt='(I7.7)') iteration
  write(myidname,fmt='(I5.5)') myid

  ! Write the fluxes: Fx Fy Fz, one file per species
  do ll = 1, Nspecies
     write (snumber,fmt='(i2.2)') ll
     filename = 'outp/datfiles/flux/F'//snumber//'_'//itername//'p'// &
          myidname//'.picard.dat'
     open(unit=1,file=filename,err=97)
     do kk = 1, Nz_local+2
        do jj = 1, Ny_local+2
           do ii = 1, Nx_local+2
              write(1,fmt='(E11.4E3,a,E11.4E3,a,E11.4E3)',err=98) &
                   F(ll,1,ii,jj,kk), ' ', F(ll,2,ii,jj,kk), ' ', &
                   F(ll,3,ii,jj,kk)
           end do
        end do
     end do
  end do
  close(1,err=99)

  return

! Error handling section
97 write (*,*) 'DumpFlux: error in open statement, file ', filename
  goto 100
98 write (*,*) 'DumpFlux: error in write statement, file ', filename
  close(1,err=99)
  goto 100
99 write (*,*) 'DumpFlux: error in close statement, file ', filename
100 if (attemptno<=Nretries) then
     call DumpFlux( F, Nx_local, Ny_local, Nz_local, &
          iteration, Nspecies, myid, attemptno, Nretries)
  end if

end subroutine DumpFlux

!------------------------------------------------------------------------

recursive subroutine DumpEprobe(Nprobes, thisprobe, fields_collected, &
     attemptno, Nretries)

  use SpecificTypes
  use ifport

  implicit none

  integer Nprobes, fields_collected, attemptno, Nretries
  type(EandPprobe) thisprobe
 

  character filename*69, tmpfile*62, fnumber*7, tnumber*7,  pnumber*4, form*20
  integer ifi, ierr

  attemptno = attemptno + 1

! Electric field
  ! Assemble filename
  write (fnumber,fmt='(i7.7)') thisprobe%iterations(1)
  write (tnumber,fmt='(i7.7)') thisprobe%iterations(fields_collected)
  write (pnumber,fmt='(i4.4)') thisprobe%probeno
  tmpfile = 'outp/datfiles/Eprobe/XXEprobe' &
       //'xxxxxxx'//'p'//pnumber//'.picard.dat'
  filename = 'outp/datfiles/Eprobe/Eprobe' &
       //fnumber//'To'//tnumber//'p'//pnumber//'.picard.dat'


  ! open file and write Ebuffer
  open(unit=1,file=tmpfile,status='replace',err=97)
  do ifi = 1, fields_collected
     write (1,fmt='(i7.7,a,E13.6E3,a,E13.6E3,a,E13.6E3)',err=98) &
          thisprobe%iterations(ifi), ' ', thisprobe%E(1,ifi), ' ', &
          thisprobe%E(2,ifi), ' ', thisprobe%E(3,ifi)
  end do
  close(1,err=99)
  ierr = rename(tmpfile,filename)
  if ( ierr>0 ) then
     write (*,*) 'DumpEprobe: error in rename statement'
     goto 100
  end if

  return

! Error handling section
97 write (*,*) 'DumpEprobe: error in open statement, file ', tmpfile
  goto 100
98 write (*,*) 'DumpEprobe: error in write statement, file ', tmpfile
  close(1,err=99)
  goto 100
99 write (*,*) 'DumpEprobe: error in close statement, file ', tmpfile
100 if (attemptno<=Nretries) then
     call DumpEprobe(Nprobes, thisprobe, fields_collected, &
     attemptno, Nretries)
  end if

end subroutine DumpEprobe

!------------------------------------------------------------------------

recursive subroutine DumpPprobe(Nprobes, thisprobe, particles, num_local, &
     iteration, myid, attemptno, Nretries)

  use SpecificTypes
  use ifport

  implicit none

  integer Nprobes, num_local, iteration, myid, attemptno, Nretries
  type(EandPprobe) thisprobe
  type(particlearrays) particles

  real*8 rp2
  character filename*69, tmpfile*62, tnumber*7,  pnumber*4, form*20
  integer ii, jj, ierr

  attemptno = attemptno + 1

  ! Assemble filename
  write (tnumber,fmt='(i7.7)') iteration
  write (pnumber,fmt='(i4.4)') thisprobe%probeno
  tmpfile = 'outp/datfiles/Pprobe/XXPprobe' &
       //'_iter'//tnumber//'p'//pnumber//'.picard.dat'
  filename = 'outp/datfiles/Pprobe/Pprobe' &
       //'_iter'//tnumber//'p'//pnumber//'.picard.dat'


  ! open file and write particles close enough to the probe centre
  rp2 = thisprobe%rprobe**2.0d0
  open(unit=1,file=tmpfile,status='replace',err=97)
  do ii = 1, num_local
     if (sum((particles%coordinates(1:3,ii)-thisprobe%rc)**2.0d0) <= rp2) then
        write(1,fmt='(6(E11.4E3,a),I3.3)',err=98) &
             (particles%coordinates(jj,ii),' ',jj=1,6), particles%species(ii)
     end if
  end do

  close(1,err=99)
  ierr = rename(tmpfile,filename)
  if ( ierr>0 ) then
     write (*,*) 'DumpPprobe: error in rename statement'
     goto 100
  end if
  
  return

! Error handling section
97 write (*,*) 'DumpPprobe: error in open statement, file ', &
        tmpfile, ' myid=', myid
  goto 100
98 write (*,*) 'DumpPprobe: error in write statement, file ', &
        tmpfile, ' myid=', myid
  close(1,err=99)
  goto 100
99 write (*,*) 'DumpPprobe: error in close statement, file ', &
        tmpfile, 'myid=', myid
100 if (attemptno<=Nretries) then
     call DumpPprobe(Nprobes, thisprobe, particles, num_local, &
     iteration, myid, attemptno, Nretries)
  end if

end subroutine DumpPprobe

!------------------------------------------------------------------------
