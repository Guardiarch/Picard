!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! The program 'picard' was developed by Jesper Lindkvist in
! 2016-2018 using resources provided by the Swedish National
! Infrastructure for Computing (SNIC) at the High Performance
! Computing Center North (HPC2N), Umeå University, Sweden.
! Jesper Lindkvist was funded by the Swedish National Space
! Board (SNSB project 201/15).
! @author    :  Jesper Lindkvist
! Email      :  jesper.lindkvist@umu.se
!
! Some updates were made by Herbert Gunell in 2018.
! Email      :  herbert.gunell@physics.org
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program picard

  use SpecificTypes
  use mpi
  implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1.0 Initialization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  logical :: writeParticles, writeFields, xflow, yflow, zflow, &
       updateEfield, updateBfield, startfromdumpfile
  integer :: Nx, Ny, Nz, iteration, iter_end, iter_start, &
       ii, jj, kk, ss, num_local, num_global, maxit, &
       dump_period_dump, dump_period_particles, dump_period_pprobes, &
       dump_period_fields, write_period_screen,Nretries, attempts, Nprobes, &
       fields_collected
  real*8  :: Lx, Ly, Lz, Lx_min, Lx_max, Ly_min, Ly_max, Lz_min, Lz_max, & 
       Lx_local, Ly_local, Lz_local, xmin, xmax, ymin, ymax, &
       zmin, zmax, dxyz(3), &
       starttime, endtime, time, dt, maxerr, E0(4), B0(4), &
       uE_local, uE_global, uB_local, uB_global, uP_local, uP_global, &
       uEB0_global, uE0_global, uB0_global, uP0_global, u0_global, &
       dV, ve(4), nucleusradius, flatradius, Galandradius, fadeoutradius

  type(particlearrays) particles
  type(particlespecies), allocatable :: species(:)
  type(EandPprobe), allocatable :: probe(:)


! For MPI handling purposes
  integer :: myid, Nprocs, ierr
  integer, allocatable :: seed(:)



  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1.1 DOMAIN PARAMETERS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! should be multiple of 2
  integer iprocs, jprocs, kprocs
  ! should be equal to INTEGER*8-2, e.g., 6, 30, 54, 78, 102
  integer Nx_local, Ny_local, Nz_local
  ! number of species included in the simulation
  integer Nspecies
  ! should be a multiple of 8
  integer, parameter :: max_per_proc = 2**24


! Physical and mathematical constants in SI units
  ! Magnetic constant [ SI ]
  real*8, parameter :: mu0  = 1.25663706143591729538505735331180d-6 
  real*8, parameter :: eps0 = 8.854187817d-12 ! Electric constant [ SI ]
  real*8, parameter :: kB   = 1.380649d-23  ! Bolzmann's constant [J/K]


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1.2 Domain Definitions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Inverse distance matrix
  real*8, allocatable :: D(:,:,:)

! Electrostatic potential
  real*8, allocatable :: U(:,:,:)

! Electric field matrices
  real*8, allocatable :: E(:,:,:,:)

! Magnetic field matrices
  real*8, allocatable :: B(:,:,:,:)

! Number of particles for each species.
  real*8, allocatable :: N(:,:,:,:)

! Particle flux for each species.
  real*8, allocatable :: F(:,:,:,:,:)

! Charge in C for all species.
  real*8, allocatable :: C(:,:,:)

! Number of retries for disk access 
  parameter(Nretries = 10)

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1.3 Get parameters from file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Eventually all parameters should come from the file, but
  ! for now this is all that has been implemented
  call GetGeneralInput(iter_end, startfromdumpfile, dump_period_dump, &
       dump_period_particles, dump_period_pprobes, &
       dump_period_fields, write_period_screen, &
       nucleusradius, flatradius, Galandradius, fadeoutradius, &
       iprocs, jprocs, kprocs, Nx_local, Ny_local, Nz_local, Nspecies, &
       xmin, xmax, ymin, ymax, zmin, zmax, dt, B0, Nprobes)

  ! Allocate the array for the species number of all particles
  allocate( particles%species(max_per_proc) )
  ! the array with positions and velocities
  allocate( particles%coordinates(6, max_per_proc) )
  ! and the species specific information.
  allocate( species(Nspecies) )

  ! Now, read the input file again to get the species specific parameters.
  call GetSpecies(Nspecies, species)
  
  ! Allocation of memory that depends on parameters just read
  allocate( D(2*Nx_local*iprocs+1,2*Ny_local*jprocs+1,2*Nz_local*kprocs+1) )
  allocate( U(Nx_local+2,Ny_local+2,Nz_local+2) )
  allocate( E(3,Nx_local+2,Ny_local+2,Nz_local+2) )
  allocate( B(3,Nx_local+2,Ny_local+2,Nz_local+2) )
  allocate( N(Nspecies,Nx_local+2,Ny_local+2,Nz_local+2) )
  allocate( F(Nspecies,3,Nx_local+2,Ny_local+2,Nz_local+2) )
  allocate( C(Nx_local+2,Ny_local+2,Nz_local+2) )

  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2.0 Initialize MPI
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call MPI_init( ierr )
  call MPI_comm_rank( MPI_comm_world, myid, ierr )
  call MPI_comm_size( MPI_comm_world, Nprocs, ierr )

  if (Nprocs .ne. iprocs*jprocs*kprocs) then
    if (myid .eq. 0) then
       write(*,*) 'Wrong Number of Processes!'
       write (*,*) 'Nprocs=', Nprocs
       write (*,*) 'iprocs=', iprocs
       write (*,*) 'jprocs=', jprocs
       write (*,*) 'kprocs=', kprocs
    end if
    call MPI_barrier(MPI_comm_world, ierr)
  end if


  call random_seed( size = Nprocs )
  allocate( seed(Nprocs) )
  call random_seed( put = seed )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2.1.0 SIMULATION PARAMETERS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! LOGICAL
  writeParticles = .true.  ! should we write particles to disk?
  writeFields = .true.  ! should we write fields to disk?
  xflow = .true.  ! should we have non-periodic in x?
  yflow = .true.  ! should we have non-periodic in y?
  zflow = .true.  ! should we have non-periodic in z?
  updateEfield = .true.  ! should we update E-field?
  updateBfield = .false.  ! should we update B-field? NOT IMPLEMENTED

! INTEGER
  maxit = 2**20  ! maximum number of iterations in Poisson solver

! REAL
  maxerr = 2.0d0**(-16.0d0)           ! maximum error in Poisson solver


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2.2 Simulation Definitions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! The layout of the processes for Nx_local = 6 looks like:
!
! 
! Proc N...
! Proc 'myid':   
!
!              Lx_local   
!      |<--------------------->|  
!  |_______________________________|
!  |   |   |   |   |   |   |   |   |
!            
!      |                       | 
!    Lx_min                  Lx_max
!
!
! Proc N+1... 
! Proc 'myid + (mod(myid+1, iprocs) - mod(myid, iprocs))':
!
!                                      Lx_local                 
!                              |<--------------------->|
!                          |_______________________________|
!                          |   |   |   |   |   |   |   |   |



! Number of cells of domain without guard cells in x,y,z
  Nx = Nx_local*iprocs
  Ny = Ny_local*jprocs
  Nz = Nz_local*kprocs

! Length of domain in x,y,z
  Lx = xmax - xmin
  Ly = ymax - ymin
  Lz = zmax - zmin

! Local length without guard cells of x,y,z
  Lx_local = Lx/dble(iprocs)
  Ly_local = Ly/dble(jprocs)
  Lz_local = Lz/dble(kprocs)

! Local min/max. Numbering through x, then y, and last z.

  Lx_min = xmin   + Lx_local*mod(myid+iprocs, iprocs)
  Lx_max = Lx_min + Lx_local
  Ly_min = ymin   + Ly_local*mod(myid/iprocs+jprocs, jprocs)
  Ly_max = Ly_min + Ly_local
  Lz_min = zmin   + Lz_local*mod(myid/(iprocs*jprocs)+kprocs, kprocs)
  Lz_max = Lz_min + Lz_local

! Cell sizes
  dxyz(1) = Lx_local/dble(Nx_local)
  dxyz(2) = Ly_local/dble(Ny_local)
  dxyz(3) = Lz_local/dble(Nz_local)
  dV = dxyz(1)*dxyz(2)*dxyz(3)

! Let's allocate the probe buffers, and get the probe-specific information
  ! We create buffers for dump_period_dump time steps, which means that
  ! we save the probe data every time we dump files we can start from.
  ! This is meant to ensure there are no gaps in the time series.
  allocate( probe(Nprobes) )
  do ii = 1, Nprobes
     probe(ii)%probeno=ii
     allocate(probe(ii)%E(3,dump_period_dump))
     allocate(probe(ii)%iterations(dump_period_dump))
  end do
  attempts = 0
  call GetProbes(Nprobes, probe, dxyz, &
       Lx_min, Lx_max, Ly_min, Ly_max, Lz_min, Lz_max, attempts, Nretries)
  
  ! The electrons in species 1 quasi-neutralise the solar wind ions
  ! I removed the computation of species 1 density. If the user
  ! wants quasi-neutrality, the user has to make it so in the
  ! input file.
  ! I still don't like the ipp, epp, ppc jungle, but now that it
  ! works for what I want it to do, I'll leave it this way. HG.
  
  species%epp   = max(idint(0.5d0-species%charge/species(1)%charge ), 0 )
! I'm killing off the hpp functionality here, by setting hpp=0 for all species,
! and it won't be resurrected unless Jesper can convince me it is useful. HG.
!!$  species%hpp   = max(idint(0.5d0-species%charge/species(2)%charge ), 0 )
!!$  species(1)%hpp= 0
  species%hpp = 0

  do ss = 3, Nspecies
     if (species(ss)%upstreamdensity > 0.0d0) then
        species(ss)%ipp = species(ss)%ppc
     else
        species(ss)%ipp = 0
     end if
  end do

  species(2)%ipp  = species(2)%ppc - &
       sum(species(3:)%ipp*species(3:)%hpp)
  species(1)%ipp  = 0

  species(1)%ppc  = sum(species(2:)%ipp*species(2:)%epp)

  species%n2p = 0.0d0
  do ss = 1, Nspecies
     if (species(ss)%upstreamdensity>0.0d0) then
        if (species(ss)%ppc .gt. 0) then
           species(ss)%n2p = dV*species(ss)%upstreamdensity * &
                dble(species(ss)%ppc)**(-1.0d0)
        else
           species(ss)%n2p = species(1)%n2p
        end if
     else
        species(ss)%n2p = species(1)%n2p * species(1)%ppc / species(ss)%ppc
     end if
  end do


! q/m [#]
  species%eta = species%charge / species%mass

! Velocity and B-field
  ve(1) = sum( species%ipp*(species%n2p*species%mass + &
       species%epp*species(1)%n2p*species(1)%mass + &
       species%hpp*species(2)%n2p*species(2)%mass)*species%v0(1) ) &
       / sum( species%ipp*(species%n2p*species%mass+species%epp* &
       species(1)%n2p*species(1)%mass + &
       species%hpp*species(2)%n2p*species(2)%mass) )
  ve(2) = sum( species%ipp*(species%n2p*species%mass + &
       species%epp*species(1)%n2p*species(1)%mass + &
       species%hpp*species(2)%n2p*species(2)%mass)*species%v0(2) ) &
       / sum( species%ipp*(species%n2p*species%mass+species%epp* &
       species(1)%n2p*species(1)%mass + &
       species%hpp*species(2)%n2p*species(2)%mass) )
  ve(3) = sum( species%ipp*(species%n2p*species%mass + &
       species%epp*species(1)%n2p*species(1)%mass + &
       species%hpp*species(2)%n2p*species(2)%mass)*species%v0(3) ) &
       / sum( species%ipp*(species%n2p*species%mass+species%epp* &
       species(1)%n2p*species(1)%mass + &
       species%hpp*species(2)%n2p*species(2)%mass) )

  do ss= 1, Nspecies
    species(ss)%v0(4) = dsqrt( sum( species(ss)%v0(1:3)**2.0d0 ) )
  end do
  ve(4)   = dsqrt( sum( ve(1:3)**2.0d0 ) )
  B0(4)   = dsqrt( sum( B0(1:3)**2.0d0 ) )
  B(1,:,:,:) = B0(1)
  B(2,:,:,:) = B0(2)
  B(3,:,:,:) = B0(3)

  E0(1) = ve(3)*B0(2)-ve(2)*B0(3)
  E0(2) = ve(1)*B0(3)-ve(3)*B0(1)
  E0(3) = ve(2)*B0(1)-ve(1)*B0(2)
  E0(4) = dsqrt( sum( E0(1:3)**2.0d0 ) )

! Initialize Distance matrix.
  call CalcDistance(D,Nx,Ny,Nz,dxyz,xflow,yflow,zflow)   


! Initialize fields and potentials to zero.
  U                     = 0.0d0
  N                     = 0.0d0
  F                     = 0.0d0
  E                     = 0.0d0
  particles%coordinates = 0.0d0
  particles%species     = 0


! Initialize scalars to zero.
  num_local = 0
  num_global = 0
  fields_collected = 0   ! Counter of iterations saved

! For the Boris algorithm from Birdsall and Langdon for constant
! homogeneous magnetic field
  do ss = 1, 8
    if ( species(ss)%gyration ) then
      species(ss)%tB(1) = 0.5d0*dt*B0(1)*species(ss)%eta
      species(ss)%tB(2) = 0.5d0*dt*B0(2)*species(ss)%eta
      species(ss)%tB(3) = 0.5d0*dt*B0(3)*species(ss)%eta
      species(ss)%tB(4) = 0.5d0*dt*B0(4)*species(ss)%eta
    else
      species(ss)%tB = 0.0d0
    end if
  end do

  species%sB(1) = 2.0d0*species%tB(1)/(1.0d0+species%tB(4)**2.0d0)
  species%sB(2) = 2.0d0*species%tB(2)/(1.0d0+species%tB(4)**2.0d0)
  species%sB(3) = 2.0d0*species%tB(3)/(1.0d0+species%tB(4)**2.0d0)
  species%sB(4) = 2.0d0*species%tB(4)/(1.0d0+species%tB(4)**2.0d0)


! Diagnostics
  uB_local      = 0.0d0
  uB_global     = 0.0d0
  uE_local      = 0.0d0
  uE_global     = 0.0d0
  uP_local      = 0.0d0
  uP_global     = 0.0d0
  uE0_global    = 0.5d0*( E0(4)**2.0d0*eps0 ) * Lx*Ly*Lz
  uB0_global    = 0.5d0*( B0(4)**2.0d0/mu0  ) * Lx*Ly*Lz
  uEB0_global   = uE0_global + uB0_global
  uP0_global    = 0.5d0*sum(species%upstreamdensity * &
       ( species%mass*species%v0(4)**2.0d0 + &
       3.0d0*kB*species%upstreamkelvin))* Lx*Ly*Lz ! Addition of both bulk and
                                                  ! temperature (not integrated)
  u0_global     = uEB0_global + uP0_global

  call MPI_barrier( MPI_comm_world, ierr)

! Write particle status
  if (myid .eq. 0) then
    write(*,*)'****************************************************************'
    write(*,*) 'PARAM SPECIES_1   SPECIES_2   SPECIES_3   SPECIES_4   SPECIES_5'
    write(*,*) 'n2p', int(species%n2p)
    write(*,*) 'ppc', species%ppc
!    write(*,*) 'ipp', species%ipp
!    write(*,*) 'epp', species%epp
!    write(*,*) 'hpp', species%hpp
    write(*,*)'****************************************************************'
 end if


! Loading or initializing the particles in the domain and writing to disk
 if (startfromdumpfile) then
    attempts = 0
    call LoadDump(Nspecies, num_local, max_per_proc, &
         particles, myid, iter_start, attempts, Nretries)
 else
    iter_start = 0
    ! Initialise the domain with a solar wind plasma
    call ParticlesInitialize(particles, &
         Lx_min,Ly_min,Lz_min,dxyz,ve,max_per_proc,&
         num_local,Nspecies,Nx_local,Ny_local,Nz_local,species)
    ! Initialise the domain with a comet ion plasma
    call ParticlesProfile(particles, &
         Lx_min,Ly_min,Lz_min,dxyz,ve,&
         max_per_proc,num_local,Nspecies, &
         Nx_local,Ny_local,Nz_local,dt, &
         nucleusradius,flatradius,Galandradius,fadeoutradius,species)
 end if

 ! Initialise time
 time = dble(iter_start)*dt

  if (writeFields) then
     call ParticlesGridifyAll(F,N,C,particles, &
          Lx_min,Lx_max,Ly_min,Ly_max,Lz_min,Lz_max,dxyz,&
          xmin,xmax,ymin,ymax,zmin,zmax,Nspecies, &
          Nx_local,Ny_local,Nz_local,max_per_proc,num_local,&
          xflow,yflow,zflow,iprocs,jprocs,kprocs,species,myid)
     call BCpotentials(U,C,D,xmin,xmax,ymin,ymax,zmin,zmax, &
          Nx_local,Ny_local,Nz_local,Nx,Ny,Nz,&
          xflow,yflow,zflow,Nspecies,iprocs,jprocs,kprocs,myid)
     call CalcEpotential(U,C,Nx_local,Ny_local,Nz_local, &
          Lx_min,Lx_max,Ly_min,Ly_max,Lz_min,Lz_max,xflow,yflow,zflow,&
          xmin,xmax,ymin,ymax,zmin,zmax,dxyz,maxerr,maxit, &
          iprocs,jprocs,kprocs,Nspecies,myid)
     call CalcEfield(E,U,Nx_local,Ny_local,Nz_local,dxyz, &
          xflow,yflow,zflow,iprocs,jprocs,kprocs,myid)
     call BCfields(E,U,B0,E0,Nx_local,Ny_local,Nz_local, &
          xflow,yflow,zflow,iprocs,jprocs,kprocs,myid)
     call MPI_barrier( MPI_comm_world, ierr)
  end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3.0 Start of Main loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (myid .eq. 0) then
    write(*,*) 'Engaging main loop!'
  end if

  do iteration = iter_start+1, iter_end

    time = time + dt   

!!$    call ParticlesModifyExp( ) ! Not implemented. MC model could be put here.
    call BCinflow(particles, &
         Lx_min,Ly_min,Lz_min,Lx_max,Ly_max,Lz_max,dxyz,ve,&
         max_per_proc,num_local,Nspecies, &
         Nx_local,Ny_local,Nz_local,&
         xmin,xmax,ymin,ymax,zmin,zmax,xflow,yflow,zflow,species)

    call AdvancePosition(particles, dt, max_per_proc, num_local) 

    call BCoutflow(particles, &
         Lx_min,Lx_max,Ly_min,Ly_max,Lz_min,Lz_max,&
         xmin,xmax,ymin,ymax,zmin,zmax,max_per_proc,num_local,xflow,yflow,zflow)

    call ParticlesTransfer(particles, &
         Lx_min,Lx_max,Ly_min,Ly_max,Lz_min,Lz_max,&
         xmin,xmax,ymin,ymax,zmin,zmax,max_per_proc, &
         num_local,iprocs,jprocs,kprocs,myid)

    if (updateEfield) then

       call ParticlesGridify(N,C,particles, &
            Lx_min,Lx_max,Ly_min,Ly_max,Lz_min,Lz_max,dxyz,&
            xmin,xmax,ymin,ymax,zmin,zmax,Nspecies, &
            Nx_local,Ny_local,Nz_local,max_per_proc,num_local,&
            xflow,yflow,zflow,iprocs,jprocs,kprocs,species,myid)

       call BCpotentials(U,C,D,xmin,xmax,ymin,ymax,zmin,zmax, &
            Nx_local,Ny_local,Nz_local,Nx,Ny,Nz,&
            xflow,yflow,zflow,Nspecies,iprocs,jprocs,kprocs,myid)

       call CalcEpotential(U,C,Nx_local,Ny_local,Nz_local, &
            Lx_min,Lx_max,Ly_min,Ly_max,Lz_min,Lz_max,xflow,yflow,zflow,&
            xmin,xmax,ymin,ymax,zmin,zmax,dxyz,maxerr,maxit, &
            iprocs,jprocs,kprocs,Nspecies,myid)

       call CalcEfield(E,U,Nx_local,Ny_local,Nz_local,dxyz, &
            xflow,yflow,zflow,iprocs,jprocs,kprocs,myid)

       call BCfields(E,U,B0,E0,Nx_local,Ny_local,Nz_local, &
            xflow,yflow,zflow,iprocs,jprocs,kprocs,myid)

    end if

        ! collect E fields in the probe buffers
        fields_collected = fields_collected+1
        do ii = 1, Nprobes
           if (probe(ii)%mine) then
              probe(ii)%E(:,fields_collected) = &
                   E(:,probe(ii)%ic(1), probe(ii)%ic(2), probe(ii)%ic(3))
              probe(ii)%iterations(fields_collected) = iteration
           end if
        end do
     ! if the buffers are full, then dump
     if (fields_collected >= dump_period_dump) then
        do ii = 1, Nprobes
           if (probe(ii)%mine) then
              attempts = 0
              call DumpEprobe(Nprobes, probe(ii), fields_collected, &
                   attempts, Nretries)
           end if
        end do
        fields_collected = 0
     end if

    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3.1 Advance velocities and write to disk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    if (writeFields .or. writeParticles) then
       if ( ( mod(iteration, dump_period_fields) .eq. 0 ) .or. &
            ( mod(iteration, dump_period_particles) .eq. 0 ) ) then

          call AdvanceVelocityImp(particles,E,B0, &
               Lx_min,Ly_min,Lz_min,dxyz,0.5d0*dt,dt,ve,&
               num_local,max_per_proc,Nspecies, &
               Nx_local,Ny_local,Nz_local,species)

          if (writeParticles) then
             if ( mod(iteration, dump_period_particles) .eq. 0)  then
                call DumpParticles( particles, num_local, iteration, myid )
             end if
             if ( mod(iteration, dump_period_pprobes) .eq. 0)  then
                do ii = 1, Nprobes
                   if (probe(ii)%mine .and. probe(ii)%rprobe>0.0d0) then
                      attempts = 0
                      call DumpPprobe( Nprobes, probe(ii), particles, &
                           num_local, iteration, myid, attempts, Nretries)
                   end if
                end do
             end if
          end if
       
          if ( (writeFields) .and. &
               ( mod(iteration, dump_period_fields) .eq. 0 ) ) then
             call ParticlesGridifyAll(F,N,C,particles, &
                  Lx_min,Lx_max,Ly_min,Ly_max,Lz_min,Lz_max,dxyz,&
                  xmin,xmax,ymin,ymax,zmin,zmax,Nspecies, &
                  Nx_local,Ny_local,Nz_local,max_per_proc,num_local,&
                  xflow,yflow,zflow,iprocs,jprocs,kprocs,species,myid)
             attempts = 0
             call DumpEfield( E, Nx_local, Ny_local, Nz_local, &
                  iteration, Nspecies, myid, attempts, Nretries)
             attempts = 0
             call DumpPotential( U, Nx_local, Ny_local, Nz_local, &
                  iteration, Nspecies, myid, attempts, Nretries)
             attempts = 0
             call DumpDensity( N, Nx_local, Ny_local, Nz_local, &
                  iteration, Nspecies, dV, myid, attempts, Nretries)
             attempts = 0
             call DumpFlux( F, Nx_local, Ny_local, Nz_local, &
                  iteration, Nspecies, myid, attempts, Nretries)
          end if

          call AdvanceVelocityExp(particles,E,B0, &
               Lx_min,Ly_min,Lz_min,dxyz,0.5d0*dt,dt,ve,&
               num_local,max_per_proc,Nspecies, &
               Nx_local,Ny_local,Nz_local,species)

       else

          call AdvanceVelocity(particles,E,B0, &
               Lx_min,Ly_min,Lz_min,dxyz,dt,ve,&
               num_local,max_per_proc,Nspecies, &
               Nx_local,Ny_local,Nz_local,species)
    
      end if
    else

       call AdvanceVelocity(particles,E,B0, &
            Lx_min,Ly_min,Lz_min,dxyz,dt,ve, &
            num_local,max_per_proc,Nspecies, &
            Nx_local,Ny_local,Nz_local,species)

    end if
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3.2 Write information to default output
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (iteration .eq. 1) then
      if (myid .eq. 0) then
        write (*,*) 'The first iteration was a success!'
        write (*,*) ' '
        write (*,*) '****************************************************************'
        write (*,*) '   ITER           TIME   MPI_RANK  NUM_LOCAL   NUM_GLOBAL    PART_ENERGY       E_ENERGY       B_ENERGY'
      end if
    end if


    if (mod(iteration,write_period_screen) .eq. 0 ) then

       call MPI_reduce( num_local, num_global, 1, MPI_INT, MPI_SUM, 0, &
            MPI_comm_world, ierr )

       uB_local = sum( B(:,2:Nx_local+1,2:Ny_local+1,2:Nz_local+1 )**2.0d0 ) &
            *0.5d0/mu0*dV
       call MPI_reduce( uB_local, uB_global, 1, MPI_DOUBLE_PRECISION, &
            MPI_SUM, 0, MPI_comm_world, ierr )
      
       uE_local=sum((E(1,2:Nx_local+1,2:Ny_local+1,2:Nz_local+1)+E0(1))**2.0d0 &
            + (E(2,2:Nx_local+1,2:Ny_local+1,2:Nz_local+1)+E0(2))**2.0d0 &
            + (E(3,2:Nx_local+1,2:Ny_local+1,2:Nz_local+1)+E0(3))**2.0d0 &
                      )*0.5d0*eps0*dV
       call MPI_reduce( uE_local, uE_global, 1, MPI_DOUBLE_PRECISION, &
            MPI_SUM, 0, MPI_comm_world, ierr )

       uP_local = sum( species(particles%species(1:num_local))%n2p* &
            species(particles%species(1:num_local))%mass*&
            ( particles%coordinates(4,1:num_local)**2.0d0 &
            + particles%coordinates(5,1:num_local)**2.0d0 &
            + particles%coordinates(6,1:num_local)**2.0d0 )  )*0.5d0
       call MPI_reduce( uP_local, uP_global, 1, MPI_DOUBLE_PRECISION, &
            MPI_SUM, 0, MPI_comm_world, ierr )

       if (myid .eq. 0) then
          write(*,fmt='(I8,a,E12.6,a,I8,a,I8,a,I10,a,E12.6,a,E12.6,a,E12.6)') &
               iteration,'   ',time,'   ',myid,'   ',num_local,'   ', &
               num_global,'   ',uP_global/u0_global,'   ',&
               uE_global/u0_global,'   ',uB_global/u0_global
       end if
    end if

    ! At regular intervals, dump something we can start from
    if ( modulo(iteration,dump_period_dump)==0 ) then
       attempts = 0
       call DumpDump(Nspecies, num_local, max_per_proc, &
            particles, myid, iteration, attempts, Nretries)
       if (myid == 0) then
          write (*,fmt='(a,i7,a, e12.6)') 'iteration = ', iteration, &
               '  t = ', dble(iteration)*dt
       end if
    end if

end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 4.0 After Main loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! If there's something left in the probe buffers, dump it now!
if (fields_collected >= 1) then
   do ii = 1, Nprobes
      if (probe(ii)%mine) then
         attempts = 0
         call DumpEprobe(Nprobes, probe(ii), fields_collected, &
              attempts, Nretries)
      end if
   end do
   fields_collected = 0
end if

  if (myid==0) then
     write (*,*) 'Picard out'
  end if

  call MPI_finalize(ierr)
  
  stop
end program picard
