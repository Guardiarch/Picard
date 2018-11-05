subroutine GetGeneralInput(startfromdumpfile, dump_period_dump)

  implicit none

  logical startfromdumpfile
  integer dump_period_dump

  integer i, j, k
  character indata*132, filename*242

  ! Defaults
  startfromdumpfile = .false.
  dump_period_dump = 1000

  filename = 'inputpicarda1.m'
  open(unit=1,file=filename,status='old')

  read (1,'(a)') indata
  do while (index(indata,'%END')+index(indata,'%end') == 0)
     i=index(indata,'=')
     if (i>0) then
        do j=1,i-1
           if(indata(j:j) /= ' ') exit
        enddo
        ! If there is a semicolon on the line, ignore it and everything beyond.
        k = index(indata(i+1:132),';')
        if (k==0) then
           k = 133-i
        end if
        select case (indata(j:i-1))
        case ('startfromdumpfile')
           if (index(indata(i+1:i+k-1),'''yes''') > 0) then
              startfromdumpfile = .true.
           end if
        case ('dump_period_dump')
           read (indata(i+1:i+k-1),*) dump_period_dump
        case default
           if (index(indata(j:j+1),'%')==0) then
              write (*,*) 'Input quantity ', indata(j:i-1), ' is unknown.'
           end if
        end select
     end if
     read (1,'(a)') indata
  end do

  close(1)

end subroutine GetGeneralInput

!------------------------------------------------------------------------

subroutine ParticlesInitialize(real_particles,int_particles, &
     Lx_min,Ly_min,Lz_min,dxyz,kelvin,mass,v0,ve,tB,max_per_proc,&
     num_local,ppc,epp,hpp,ipp,Nspecies,Nx_local,Ny_local,Nz_local)

  implicit none

! Parameters
  integer, intent(in) :: max_per_proc, ppc(8), epp(8), hpp(8), ipp(8), Nspecies, Nx_local, Ny_local, Nz_local
  real*8, intent(in) :: Lx_min, Ly_min, Lz_min, dxyz(3), kelvin(8), mass(8), v0(4,8), ve(4), tB(4,8)

  real*8, intent(inout) :: real_particles(8,max_per_proc)
  integer, intent(inout) :: int_particles(max_per_proc), num_local

! Local variables
  integer ee, pp, ss, ii, jj, kk
  real*8 rndp(3), x0, y0, z0, vmin(3), vmid(3)


!  Graphics
!
!     xmin                  Plasma                    xmax
!      |-----------------------------------------------|
! 
! Proc N:   
!
!                Lx_local   
!      |<--------------------->|  
!  |_______________________________|
!  |   |   |   |   |   |   |   |   |
!            
!      |                       | 
!    Lx_min                  Lx_max
!
! 
!                            Proc N+1:
!                                       Lx_local                 
!                              |<--------------------->|
!                          |_______________________________|
!                          |   |   |   |   |   |   |   |   |




! Ok, let's try something here: A common problem is that we initialize
! the plasma with fully developed fluctuation levels, i.e. there will be
! a lot of electrostatic energy beacuse electrons and ions do not cancel
! out completely within each cell. Let's try to put electrons and ions
! at the exact same locations to begin with. This will give us a plasma
! with exactly the energy we want (if the chosen # particles/cell for
! each species scales with the charge density).

! Loop through each cell which is not a guard cell
do kk = 2, Nz_local+1
  do jj = 2, Ny_local+1
    do ii = 2, Nx_local+1
      ! Lower edge of cells
      x0 = Lx_min + dxyz(1)*(ii-2)
      y0 = Ly_min + dxyz(2)*(jj-2)
      z0 = Lz_min + dxyz(3)*(kk-2)
      ! Loop through all ions
      do ss = 1, Nspecies

        ! Loop through all particles to be created
        do pp = 1, ipp(ss)
          int_particles(num_local+1)        = ss
          call random_number(rndp)
          real_particles(1,num_local+1)     = x0 + dxyz(1)*rndp(1)
          real_particles(2,num_local+1)     = y0 + dxyz(2)*rndp(2)
          real_particles(3,num_local+1)     = z0 + dxyz(3)*rndp(3)
          call Gaussian( real_particles(4:6,num_local+1),v0(1:3,ss),kelvin(ss),mass(ss) )
          ! Advance velocity dt/2
          vmin(1:3) = real_particles(4:6,num_local+1) - ve(1:3)
          vmid(1) = vmin(1) + vmin(2)*tB(3,ss)-vmin(3)*tB(2,ss)
          vmid(2) = vmin(2) + vmin(3)*tB(1,ss)-vmin(1)*tB(3,ss)
          vmid(3) = vmin(3) + vmin(1)*tB(2,ss)-vmin(2)*tB(1,ss)
          if ( sum(vmid(:)**2) .gt. 0.0d0 ) then
            vmid = vmid * dsqrt( sum(vmin(:)**2)/sum(vmid(:)**2) )
          end if
          real_particles(4:6,num_local+1) = vmid(1:3) + ve(1:3)

          ! For each ion, create as many electrons as its charged state.
          do ee = 1, epp(ss)
            int_particles(num_local+1+ee)         = 1
            real_particles(:,num_local+1+ee)      = real_particles(:,num_local+1)
            call Gaussian( real_particles(4:6,num_local+1+ee),v0(1:3,ss),kelvin(1),mass(1) )
            ! Advance velocity dt/2
            vmin(1:3) = real_particles(4:6,num_local+1+ee) - ve(1:3)
            vmid(1) = vmin(1) + vmin(2)*tB(3,1)-vmin(3)*tB(2,1)
            vmid(2) = vmin(2) + vmin(3)*tB(1,1)-vmin(1)*tB(3,1)
            vmid(3) = vmin(3) + vmin(1)*tB(2,1)-vmin(2)*tB(1,1)
            if ( sum(vmid(:)**2) .gt. 0.0d0 ) then
              vmid = vmid * dsqrt( sum(vmin(:)**2)/sum(vmid(:)**2) )
            end if
            real_particles(4:6,num_local+1+ee) = vmid(1:3) + ve(1:3)
          end do
          
          do ee = 1, hpp(ss)
            int_particles(num_local+1+epp(ss)+ee)    = 2
            real_particles(:,num_local+1+epp(ss)+ee) = real_particles(:,num_local+1)
            call Gaussian( real_particles(4:6,num_local+1+epp(ss)+ee),v0(1:3,ss),kelvin(2),mass(2) )
            ! Advance velocity dt/2
            vmin(1:3) = real_particles(4:6,num_local+1+epp(ss)+ee) - ve(1:3)
            vmid(1) = vmin(1) + vmin(2)*tB(3,2)-vmin(3)*tB(2,2)
            vmid(2) = vmin(2) + vmin(3)*tB(1,2)-vmin(1)*tB(3,2)
            vmid(3) = vmin(3) + vmin(1)*tB(2,2)-vmin(2)*tB(1,2)
            if ( sum(vmid(:)**2) .gt. 0.0d0 ) then
              vmid = vmid * dsqrt( sum(vmin(:)**2)/sum(vmid(:)**2) )
            end if
            real_particles(4:6,num_local+1+epp(ss)+ee) = vmid(1:3) + ve(1:3)
          end do

          num_local = num_local+1+epp(ss)+hpp(ss)

        end do
      end do
    end do
  end do
end do

! Now we have a neutral plasma with no excess energy...
    

return
end subroutine ParticlesInitialize

