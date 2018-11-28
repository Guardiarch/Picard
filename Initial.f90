subroutine GetGeneralInput(iter_end, startfromdumpfile, dump_period_dump, &
     dump_period_particles, dump_period_pprobes, &
     dump_period_fields, write_period_screen, &
     nucleusradius, flatradius, Galandradius, fadeoutradius, &
     iprocs, jprocs, kprocs, Nx_local, Ny_local, Nz_local, Nspecies, &
     xmin, xmax, ymin, ymax, zmin, zmax, dt, B0, Nprobes)

  implicit none

  logical startfromdumpfile
  integer iter_end, dump_period_dump, dump_period_particles, &
       dump_period_pprobes, dump_period_fields, write_period_screen, &
       iprocs, jprocs, kprocs, &
       Nx_local, Ny_local, Nz_local, Nspecies, Nprobes
  real*8  nucleusradius, flatradius, Galandradius, fadeoutradius, &
       xmin, xmax, ymin, ymax, zmin, zmax, dt, B0(4)
  
  integer i, j, k
  character indata*132, filename*242

  ! Defaults
  iter_end = 10000000  ! Bitter end
  startfromdumpfile = .false. ! start from dump or not
  dump_period_dump = 1000  ! iterations between dumps one can start from
  dump_period_particles=100000 ! number of iterations between particle dumps
  dump_period_pprobes=5000 ! number of iterations between particle probe dumps
  dump_period_fields=10   ! number of iterations between field dumps
  write_period_screen=10  ! number of iterations between showing onscreen
                          ! life signs.

  nucleusradius = 5.0d0;   ! radius of the nucleus [m], 2.0d3/400.0d0
  flatradius = 1.0d1       ! radius of flat density part [m], 4.0d3/400.0d0
  Galandradius = 8.0d2     ! outer radius of Galand model
  fadeoutradius = 8.8d2    ! radius at which the cometary ion density is zero
  iprocs = 2      ! processes in x dimension, should be multiple of 2
  jprocs = 2      ! processes in y dimension, should be multiple of 2
  kprocs = 2      ! processes in z dimension, should be multiple of 2
  Nx_local = 30   ! should be equal to INTEGER*8-2, e.g., 6, 30, 54, 78, 102
  Ny_local = 30   ! should be equal to INTEGER*8-2, e.g., 6, 30, 54, 78, 102
  Nz_local = 30   ! should be equal to INTEGER*8-2, e.g., 6, 30, 54, 78, 102
  Nspecies = 3    ! number of species included in the simulation
  xmin = -240.0d0  ! physical domain
  ymin = -240.0d0  ! physical domain
  zmin = -240.0d0  ! physical domain
  xmax =  240.0d0  ! physical domain
  ymax =  240.0d0  ! physical domain
  zmax =  240.0d0  ! physical domain
  dt   =  4.0d-7   ! timestep
  B0   = 0.0d0     ! Magnetic flux density
  Nprobes = 0      ! Number of electric field probes

  filename = 'inputpicarda1.m'
  open(unit=1,file=filename,status='old',err=98)
  goto 100
  
98 write (*,*) 'GetGeneralInput: Error opening input file ', filename
  stop
  
99 write (*,*) 'GetGeneralInput: Error reading input file ', filename
  stop
  
100 read (1,'(a)',err=99) indata
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
        case ('iter_end')
           read (indata(i+1:i+k-1),*) iter_end
        case ('startfromdumpfile')
           if (index(indata(i+1:i+k-1),'''yes''') > 0) then
              startfromdumpfile = .true.
           end if
        case ('dump_period_dump')
           read (indata(i+1:i+k-1),*) dump_period_dump
        case ('dump_period_particles')
           read (indata(i+1:i+k-1),*) dump_period_particles
        case ('dump_period_pprobes')
           read (indata(i+1:i+k-1),*) dump_period_pprobes
        case ('dump_period_fields')
           read (indata(i+1:i+k-1),*) dump_period_fields
        case ('write_period_screen')
           read (indata(i+1:i+k-1),*) write_period_screen
        case ('nucleusradius')
           read (indata(i+1:i+k-1),*) nucleusradius
        case ('flatradius')
           read (indata(i+1:i+k-1),*) flatradius
        case ('Galandradius', 'galandradius')
           read (indata(i+1:i+k-1),*) Galandradius
        case ('fadeoutradius')
           read (indata(i+1:i+k-1),*) fadeoutradius
        case ('iprocs')
           read (indata(i+1:i+k-1),*) iprocs
        case ('jprocs')
           read (indata(i+1:i+k-1),*) jprocs
        case ('kprocs')
           read (indata(i+1:i+k-1),*) kprocs
        case ('Nx_local', 'nx_local')
           read (indata(i+1:i+k-1),*) Nx_local
        case ('Ny_local', 'ny_local')
           read (indata(i+1:i+k-1),*) Ny_local
        case ('Nz_local', 'nz_local')
           read (indata(i+1:i+k-1),*) Nz_local
        case ('Nspecies', 'nspecies')
           read (indata(i+1:i+k-1),*) Nspecies
        case ('xmin')
           read (indata(i+1:i+k-1),*) xmin
        case ('ymin')
           read (indata(i+1:i+k-1),*) ymin
        case ('zmin')
           read (indata(i+1:i+k-1),*) zmin
        case ('xmax')
           read (indata(i+1:i+k-1),*) xmax
        case ('ymax')
           read (indata(i+1:i+k-1),*) ymax
        case ('zmax')
           read (indata(i+1:i+k-1),*) zmax
        case ('dt')
           read (indata(i+1:i+k-1),*) dt
        case ('B0x', 'b0x')
           read (indata(i+1:i+k-1),*) B0(1)
        case ('B0y', 'b0y')
           read (indata(i+1:i+k-1),*) B0(2)
        case ('B0z', 'b0z')
           read (indata(i+1:i+k-1),*) B0(3)
        case ('Nprobes', 'nprobes')
           read (indata(i+1:i+k-1),*) Nprobes
        case default
           if (index(indata(j:j+1),'%')==0) then
              write (*,*) 'Input quantity ', indata(j:i-1), ' is unknown.'
           end if
        end select
     end if
     read (1,'(a)') indata
  end do

  close(1,err=101)
  return

101 write (*,*) 'GetGeneralInput: Error closing input file ', filename
  stop

end subroutine GetGeneralInput

!------------------------------------------------------------------------

subroutine GetSpecies(Nspecies, species)

  use SpecificTypes

  implicit none

  integer Nspecies
  type(particlespecies) species(Nspecies)

  character indata*132, filename*242
  integer ii, i, j, k

  real*8, parameter :: elc  = 1.602176634d-19  ! elementary charge [C]
  
  filename='inputpicarda1.m'
  open(unit=1,file=filename,status='old',err=98)
  goto 100
  
98 write (*,*) 'GetSpecies: Error opening input file ', filename
  stop
  
99 write (*,*) 'GetSpecies: Error reading input file ', filename
  stop
  
100 do ii = 1, Nspecies
     ! Defaults
     species(ii)%gyration        = .true.  ! ZERO B-FIELD FOR CHOSEN SPECIES.
                                           ! TRUE: F=q(E+vxB). FALSE: F=qE
     ! Only Jesper has tested the fluid option, better keep it .false.
     species(ii)%fluid           = .false. ! FLUID APPROXIMATION.
                                           ! TRUE: v_perp=FxB/(q*B^2),
                                           ! a_para = F dot B...
                                           ! FALSE: Lorentz force
     species(ii)%ppc             = 0
     species(ii)%mass            = 9.10938356d-31
     species(ii)%charge          = -1.602176634d-19 
     species(ii)%upstreamdensity = 0.0d0
     species(ii)%upstreamkelvin  = 0.0d0
     species(ii)%v0              = 0.0d0
     species(ii)%vn              = 0.0d0
     species(ii)%Qn              = 0.0d0
     species(ii)%Ek              = 0.0d0
     species(ii)%nu_i            = 0.0d0
     species(ii)%n2p             = 0
     species(ii)%tB              = 0.0d0
     species(ii)%sB              = 0.0d0
     species(ii)%ipp             = 0
     species(ii)%epp             = 0
     species(ii)%hpp             = 0
     species(ii)%eta             = 1.7588200380926763d11
     species(ii)%productspecies  = 0
     species(ii)%cometion        = .false.
     read (1,'(a)',err=99) indata
     do while (index(indata,'%SPEC')+index(indata,'%spec') == 0)
        read (1,'(a)',err=99) indata
     end do
     do while (index(indata,'%END')+index(indata,'%end') == 0)
        i=index(indata,'=')
        if (i>0) then
           do j=1,i-1
              if(indata(j:j) /= ' ') exit
           enddo
           ! If there is a semicolon on the line, ignore it and 
           ! everything beyond.
           k = index(indata(i+1:132),';')
           if (k==0) then
              k = 133-i
           end if
           select case (indata(j:i-1))
           case ('gyration')
              if (index(indata(i+1:i+k-1),'''no''') > 0) then
                 species(ii)%gyration = .false.
              end if
           case('fluid')
              if (index(indata(i+1:i+k-1),'''yes''') > 0) then
                 species(ii)%fluid = .true.
              end if
           case ('ppc')
              read (indata(i+1:i+k-1),*) species(ii)%ppc
           case ('mass')
              read (indata(i+1:i+k-1),*) species(ii)%mass
           case ('charge')
              read (indata(i+1:i+k-1),*) species(ii)%charge
           case ('upstreamdensity')
              read (indata(i+1:i+k-1),*) species(ii)%upstreamdensity
           case ('upstreamkelvin')
              read (indata(i+1:i+k-1),*) species(ii)%upstreamkelvin
           case ('v0x')
              read (indata(i+1:i+k-1),*) species(ii)%v0(1)
           case ('v0y')
              read (indata(i+1:i+k-1),*) species(ii)%v0(2)
           case ('v0z')
              read (indata(i+1:i+k-1),*) species(ii)%v0(3)
           case ('vn')
              read (indata(i+1:i+k-1),*) species(ii)%vn
           case ('Qn')
              read (indata(i+1:i+k-1),*) species(ii)%Qn
           case ('Ek')
              read (indata(i+1:i+k-1),*) species(ii)%Ek
              species(ii)%Ek = species(ii)%Ek * elc ! convert eV -> J
           case ('nu_i')
              read (indata(i+1:i+k-1),*) species(ii)%nu_i
           case ('cometion')
              if (index(indata(i+1:i+k-1),'''yes''') > 0) then
                 species(ii)%cometion = .true.
              end if
           case ('productspecies')
              read (indata(i+1:i+k-1),*) species(ii)%productspecies
           case default
              if (index(indata(j:j+1),'%')==0) then
                 write (*,*) 'Input quantity ', indata(j:i-1), ' is unknown.'
              end if
           end select
        end if
        read (1,'(a)') indata
     end do
  end do

  close(1,err=101)
  return

101 write (*,*) 'GetSpecies: Error closing input file ', filename
  stop

end subroutine GetSpecies

!------------------------------------------------------------------------

recursive subroutine GetProbes(Nprobes, probe, dxyz, &
     Lx_min, Lx_max, Ly_min, Ly_max, Lz_min, Lz_max, attemptno, Nretries)

  use SpecificTypes

  implicit none

  integer Nprobes, attemptno, Nretries
  type(EandPprobe) probe(Nprobes)
  real*8 dxyz(3), Lx_min, Lx_max, Ly_min, Ly_max, Lz_min, Lz_max
  
  character indata*132, filename*242
  integer ii, i, j, k

  attemptno = attemptno + 1

  filename='inputpicarda1.m'
  open(unit=1,file=filename,status='old',err=97)

  do ii = 1, Nprobes
     ! Defaults
     probe(ii)%xc = 1.0d0/0.0d0 ! The default should be outside 
     probe(ii)%yc = 1.0d0/0.0d0 ! the simulation box, and I think
     probe(ii)%zc = 1.0d0/0.0d0 ! this is.
     probe(ii)%mine = .false.
     probe(ii)%ic = 0
     probe(ii)%rc = 1.0d0/0.0d0
     probe(ii)%rprobe = -1.0d0
     probe(ii)%E = 0.0d0
     probe(ii)%iterations = 0
     read (1,'(a)',err=98) indata
     do while (index(indata,'%PROB')+index(indata,'%prob') == 0)
        read (1,'(a)') indata
     end do
     do while (index(indata,'%END')+index(indata,'%end') == 0)
        i=index(indata,'=')
        if (i>0) then
           do j=1,i-1
              if(indata(j:j) /= ' ') exit
           enddo
           ! If there is a semicolon on the line, ignore it and 
           ! everything beyond.
           k = index(indata(i+1:132),';')
           if (k==0) then
              k = 133-i
           end if
           select case (indata(j:i-1))
           case ('xc')
              read (indata(i+1:i+k-1),*) probe(ii)%xc
           case ('yc')
              read (indata(i+1:i+k-1),*) probe(ii)%yc
           case ('zc')
              read (indata(i+1:i+k-1),*) probe(ii)%zc
           case ('rprobe')
              read (indata(i+1:i+k-1),*) probe(ii)%rprobe
           case default
              if (index(indata(j:j+1),'%')==0) then
                 write (*,*) 'Input quantity ', indata(j:i-1), ' is unknown.'
              end if
           end select
        end if
        read (1,'(a)',err=98) indata
     end do
  end do

  close(1,err=99)

  ! Find our place in the grid here.
  do ii = 1, Nprobes
     if ( (probe(ii)%xc>=Lx_min) .and. (probe(ii)%xc<Lx_max) .and. &
          (probe(ii)%yc>=Ly_min) .and. (probe(ii)%yc<Ly_max) .and. &
          (probe(ii)%zc>=Lz_min) .and. (probe(ii)%zc<Lz_max) ) then
        ! This probe is mine
        probe(ii)%mine = .true.
        ! Indeces of nearest cell centre
        probe(ii)%ic(1) = 1 + nint((probe(ii)%xc-Lx_min)/dxyz(1)+0.5d0 )
        probe(ii)%ic(2) = 1 + nint((probe(ii)%yc-Ly_min)/dxyz(2)+0.5d0 )
        probe(ii)%ic(3) = 1 + nint((probe(ii)%zc-Lz_min)/dxyz(3)+0.5d0 )
        ! real position to which this corresponds
        probe(ii)%rc(1) = (dble(probe(ii)%ic(1)-1)-0.5d0)*dxyz(1)
        probe(ii)%rc(2) = (dble(probe(ii)%ic(2)-1)-0.5d0)*dxyz(2)
        probe(ii)%rc(3) = (dble(probe(ii)%ic(3)-1)-0.5d0)*dxyz(3)
     else
        ! This probe isn't mine
        probe(ii)%mine = .false.
        deallocate( probe(ii)%E )
        deallocate( probe(ii)%iterations )
     end if
  end do

  return

! Error handling section
97 write (*,*) 'GetProbes: error in open statement'
  goto 100
98 write (*,*) 'GetProbes: error in read statement'
  close(1,err=99)
  goto 100
99 write (*,*) 'GetProbes: error in close statement'
100 if (attemptno<=Nretries) then
     call GetProbes(Nprobes, probe, dxyz, &
          Lx_min, Lx_max, Ly_min, Ly_max, Lz_min, Lz_max, attemptno, Nretries)
  end if

end subroutine GetProbes

!------------------------------------------------------------------------

subroutine ParticlesInitialize(particles, &
     Lx_min,Ly_min,Lz_min,dxyz,ve,max_per_proc,&
     num_local,Nspecies,Nx_local,Ny_local,Nz_local,species)

  use SpecificTypes
  implicit none

! Parameters
  integer, intent(in) :: max_per_proc, Nspecies, Nx_local, Ny_local, Nz_local
  real*8, intent(in) :: Lx_min, Ly_min, Lz_min, dxyz(3), ve(4)

  integer, intent(inout) :: num_local
  type(particlearrays) particles
  type(particlespecies) species(Nspecies)

! Local variables
  integer ee, pp, ss, ii, jj, kk
  real*8 rndp(3), x0, y0, z0, vmin(3), vmid(3)

! This subroutine is designed to create a uniform plasma of the particle
! species that continuously come in from the boundary during the simulation. 
! In other words, it initialises the simulation domain with a solar wind plama.

! Following Hurtig et al. we put electrons and ions at the exact same locations.

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
        do pp = 1, species(ss)%ipp
          particles%species(num_local+1)        = ss
          call random_number(rndp)
          particles%coordinates(1,num_local+1)     = x0 + dxyz(1)*rndp(1)
          particles%coordinates(2,num_local+1)     = y0 + dxyz(2)*rndp(2)
          particles%coordinates(3,num_local+1)     = z0 + dxyz(3)*rndp(3)
          call Gaussian( particles%coordinates(4:6,num_local+1),species(ss)%v0(1:3), &
               species(ss)%upstreamkelvin,species(ss)%mass )
          ! Advance velocity dt/2
          vmin(1:3) = particles%coordinates(4:6,num_local+1) - ve(1:3)
          vmid(1) = vmin(1)+vmin(2)*species(ss)%tB(3)-vmin(3)*species(ss)%tB(2)
          vmid(2) = vmin(2)+vmin(3)*species(ss)%tB(1)-vmin(1)*species(ss)%tB(3)
          vmid(3) = vmin(3)+vmin(1)*species(ss)%tB(2)-vmin(2)*species(ss)%tB(1)
          if ( sum(vmid(:)**2) .gt. 0.0d0 ) then
            vmid = vmid * dsqrt( sum(vmin(:)**2)/sum(vmid(:)**2) )
          end if
          particles%coordinates(4:6,num_local+1) = vmid(1:3) + ve(1:3)

          ! For each ion, create as many electrons as its charged state.
          do ee = 1, species(ss)%epp
            particles%species(num_local+1+ee)    = 1
            particles%coordinates(:,num_local+1+ee) = particles%coordinates(:,num_local+1)
            call Gaussian( particles%coordinates(4:6,num_local+1+ee), &
                 species(ss)%v0(1:3),species(1)%upstreamkelvin,species(1)%mass )
            ! Advance velocity dt/2
            vmin(1:3) = particles%coordinates(4:6,num_local+1+ee) - ve(1:3)
            vmid(1) = vmin(1)+vmin(2)*species(1)%tB(3)-vmin(3)*species(1)%tB(2)
            vmid(2) = vmin(2)+vmin(3)*species(1)%tB(1)-vmin(1)*species(1)%tB(3)
            vmid(3) = vmin(3)+vmin(1)*species(1)%tB(2)-vmin(2)*species(1)%tB(1)
            if ( sum(vmid(:)**2) .gt. 0.0d0 ) then
              vmid = vmid * dsqrt( sum(vmin(:)**2)/sum(vmid(:)**2) )
            end if
            particles%coordinates(4:6,num_local+1+ee) = vmid(1:3) + ve(1:3)
          end do
          
          do ee = 1, species(ss)%hpp
            particles%species(num_local+1+species(ss)%epp+ee)    = 2
            particles%coordinates(:,num_local+1+species(ss)%epp+ee) = &
                 particles%coordinates(:,num_local+1)
            call Gaussian(particles%coordinates(4:6,num_local+1+species(ss)%epp+ee), &
                 species(ss)%v0(1:3),species(2)%upstreamkelvin,species(2)%mass )
            ! Advance velocity dt/2
            vmin(1:3)=particles%coordinates(4:6,num_local+1+species(ss)%epp+ee) - ve(1:3)
            vmid(1) = vmin(1)+vmin(2)*species(2)%tB(3)-vmin(3)*species(2)%tB(2)
            vmid(2) = vmin(2)+vmin(3)*species(2)%tB(1)-vmin(1)*species(2)%tB(3)
            vmid(3) = vmin(3)+vmin(1)*species(2)%tB(2)-vmin(2)*species(2)%tB(1)
            if ( sum(vmid(:)**2) .gt. 0.0d0 ) then
               vmid = vmid * dsqrt( sum(vmin(:)**2)/sum(vmid(:)**2) )
            end if
            particles%coordinates(4:6,num_local+1+species(ss)%epp+ee) = &
                 vmid(1:3) + ve(1:3)
          end do

          num_local = num_local+1+species(ss)%epp+species(ss)%hpp

        end do
      end do
    end do
  end do
end do

return
end subroutine ParticlesInitialize

!------------------------------------------------------------------------

subroutine ParticlesProfile(particles, &
     Lx_min,Ly_min,Lz_min,dxyz,ve,&
     max_per_proc,num_local,Nspecies, &
     Nx_local,Ny_local,Nz_local,dt, &
     nucleusradius,flatradius,Galandradius,fadeoutradius,species)

  use SpecificTypes
  implicit none

  interface
     real*8 function Galand(nu, Q, r, nucleusradius, u)
       real*8, intent(in) :: nu, Q, r, nucleusradius, u
     end function Galand
  end interface

  ! Parameters
  integer, intent(in) :: max_per_proc, Nspecies, Nx_local, Ny_local, Nz_local
  real*8, intent(in) :: Lx_min, Ly_min, Lz_min, dxyz(3), ve(4), &
       dt, nucleusradius, flatradius,Galandradius,fadeoutradius

  integer, intent(inout) :: num_local
  type(particlearrays) particles
  type(particlespecies) species(Nspecies)

  real*8, parameter :: pi  = 3.141592653589793d0

  ! Local variables
  integer ee, pp, ss, ii, jj, kk, num_create
  real*8 rnd, rndp(7), x0, y0, z0, vmin(3), vmid(3), r1, nden, dV
  real*8, allocatable :: ui(:), ue(:), uh(:)

  allocate (ui(Nspecies) )
  allocate (ue(Nspecies) )
  allocate (uh(Nspecies) )
  dV = dxyz(1)*dxyz(2)*dxyz(3)
  ue = 0.0d0
  uh = 0.0d0
  ui = 0.0d0

  do ss = 1, Nspecies
     ue(ss)=dsqrt(2.0d0*species(ss)%mass*species(ss)%Ek* &
          (species(1)%mass*(species(ss)%mass+species(1)%mass))**(-1.0d0))
     uh(ss)=dsqrt(2.0d0*species(ss)%mass*species(ss)%Ek* &
          (species(2)%mass*(species(ss)%mass+species(2)%mass))**(-1.0d0))
     if (species(ss)%epp .gt. 0) then
        ui(ss) = dsqrt( 2.0d0*species(1)%mass*species(ss)%Ek* &
             (species(ss)%mass*(species(ss)%mass+species(1)%mass))**(-1.0d0) )
     else if (species(ss)%hpp .gt. 0) then
        ui(ss) = dsqrt( 2.0d0*species(2)%mass*species(ss)%Ek * &
             (species(ss)%mass*(species(ss)%mass+species(2)%mass))**(-1.0d0) )
     else
        ui(ss)  = 0.0d0
     end if
  end do
  
  ! Loop through each cell which is not a guard cell
  do kk = 2, Nz_local+1
     do jj = 2, Ny_local+1
        do ii = 2, Nx_local+1
           ! Lower edge of cells
           x0 = Lx_min + dxyz(1)*(ii-2)
           y0 = Ly_min + dxyz(2)*(jj-2)
           z0 = Lz_min + dxyz(3)*(kk-2)
           r1 =dsqrt((x0+0.5d0*dxyz(1))**2.0d0 + (y0+0.5d0*dxyz(2))**2.0d0 + &
                (z0+0.5d0*dxyz(3))**2.0d0)
           if ( (4.0d0*r1)**3.0d0 .lt. dV) then
              r1 = ( dsqrt((x0+0.25d0*dxyz(1))**2.0d0 + &
                   (y0+0.25d0*dxyz(2))**2.0d0 + &
                   (z0+0.25d0*dxyz(3))**2.0d0) &
                   + dsqrt((x0+0.25d0*dxyz(1))**2.0d0 + &
                   (y0+0.25d0*dxyz(2))**2.0d0 + &
                   (z0+0.75d0*dxyz(3))**2) &
                   + dsqrt((x0+0.25d0*dxyz(1))**2.0d0 + &
                   (y0+0.75d0*dxyz(2))**2.0d0 + &
                   (z0+0.25d0*dxyz(3))**2) + &
                   dsqrt((x0+0.75d0*dxyz(1))**2.0d0 + &
                   (y0+0.25d0*dxyz(2))**2.0d0 + &
                   (z0+0.25d0*dxyz(3))**2) + &
                   dsqrt((x0+0.25d0*dxyz(1))**2.0d0 + &
                   (y0+0.75d0*dxyz(2))**2.0d0 + &
                   (z0+0.75d0*dxyz(3))**2) + &
                   dsqrt((x0+0.75d0*dxyz(1))**2.0d0 + &
                   (y0+0.25d0*dxyz(2))**2.0d0 + &
                   (z0+0.75d0*dxyz(3))**2) + &
                   dsqrt((x0+0.75d0*dxyz(1))**2.0d0 + &
                   (y0+0.75d0*dxyz(2))**2.0d0 + &
                   (z0+0.25d0*dxyz(3))**2) + &
                   dsqrt((x0+0.75d0*dxyz(1))**2.0d0 + &
                   (y0+0.75d0*dxyz(2))**2.0d0 + &
                   (z0+0.75d0*dxyz(3))**2) )*0.125d0
           end if

           do ss = 1, Nspecies
              if ((species(ss)%cometion) .and. (species(ss)%vn .gt. 0.0d0)) then
                 ! Create a density following an approximate 1/r profile
                 ! (Galand et al., 2016), except inside flatradius, where
                 ! it is flat, and outside galandradius, where it fades
                 ! to zero which it reaches at fadeoutradius
                 if ( r1 > fadeoutradius ) then
                    nden = 0.0d0
                 elseif ( r1 > Galandradius ) then
                    nden = ( 1.0d0 - (r1-Galandradius) / &
                         (fadeoutradius-Galandradius) ) * &
                         Galand(species(ss)%nu_i, species(ss)%Qn, &
                         Galandradius, nucleusradius, species(ss)%vn)
                 elseif ( r1 > flatradius ) then
                    nden = Galand(species(ss)%nu_i, species(ss)%Qn, r1, &
                         nucleusradius, species(ss)%vn)
                 else
                    nden=Galand(species(ss)%nu_i,species(ss)%Qn, &
                         flatradius,nucleusradius,species(ss)%vn)
                 end if
              else
                 nden = 0.0d0
              end if

              if ( species(ss)%n2p .gt. 0.0d0 ) then
                 call random_number(rnd)
                 num_create = idint( nden*dV/species(ss)%n2p + rnd )
              else
                 num_create = 0
              end if
              
              if (num_local+(1+species(ss)%epp+species(ss)%hpp)*num_create > &
                   max_per_proc) then
                 write (*,*) 'Cannot create more particles than what fills ', &
                      'the allocated memory, sorry.'
                 write (*,*) 'num_local+(1+species(ss)%epp+species(ss)%hpp)', &
                      '*num_create=', num_local+(1+species(ss)%epp+ &
                      species(ss)%hpp)*num_create, ' nden=',nden, &
                      ' num_local=', num_local,' num_create=', num_create
                 stop
              end if
              
              do pp = 1, num_create
                 call random_number(rndp)
                 particles%species(num_local+1)      = ss
                 particles%coordinates(1,num_local+1)   = x0 + dxyz(1)*rndp(1)
                 particles%coordinates(2,num_local+1)   = y0 + dxyz(2)*rndp(2)
                 particles%coordinates(3,num_local+1)   = z0 + dxyz(3)*rndp(3)
                 particles%coordinates(4:6,num_local+1) = species(ss)%vn * &
                      particles%coordinates(1:3,num_local+1) / dsqrt(sum( &
                      particles%coordinates(1:3,num_local+1)**2.0d0)) + &
                      ui(ss) * rndp(4:6)/dsqrt(sum(rndp(4:6)**2.0d0))
                 particles%coordinates(1:3,num_local+1) = &
                      particles%coordinates(1:3,num_local+1) + &
                      particles%coordinates(4:6,num_local+1)*dt*rndp(7)
                 ! For each ion, create as many electrons as its charged state.
                 do ee = 1, species(ss)%epp
                    particles%species(num_local+1+ee) = &
                         species(ss)%productspecies
                    particles%coordinates(4:6,num_local+1+ee)   = &
                         particles%coordinates(4:6,num_local+1) - &
                         (ue(ss)+ui(ss))*rndp(4:6)/dsqrt(sum(rndp(4:6)**2.0d0))
                    particles%coordinates(1:3,num_local+1+ee)   = &
                         particles%coordinates(1:3,num_local+1) &
                         - particles%coordinates(4:6,num_local+1)*dt*rndp(7) &
                         + particles%coordinates(4:6,num_local+1+ee)*dt*rndp(7)
                 end do
                 do ee = 1, species(ss)%hpp
                    particles%species(num_local+1+species(ss)%epp+ee)       = 2
                    particles%coordinates(4:6,num_local+1+species(ss)%epp+ee) &
                         = particles%coordinates(4:6,num_local+1) - &
                         (uh(ss)+ui(ss))*rndp(4:6)/dsqrt(sum(rndp(4:6)**2.0d0))
                    particles%coordinates(1:3,num_local+1+species(ss)%epp+ee) &
                         = particles%coordinates(1:3,num_local+1) - &
                         particles%coordinates(4:6,num_local+1)*dt*rndp(7) + &
                         particles%coordinates(4:6,num_local+1+ &
                         species(ss)%epp+ee)*dt*rndp(7)
                 end do
                 num_local = num_local+1+species(ss)%epp+species(ss)%hpp
              end do

           end do
        end do
     end do
  end do

  return
end subroutine ParticlesProfile

!------------------------------------------------------------------------
