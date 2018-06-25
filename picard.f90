
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The program 'picard' was developed by Jesper Lindkvist in 2016-2018 using resources
! provided by the Swedish National Infrastructure for Computing (SNIC) at the 
! High Performance Computing Center North (HPC2N), Umeå University, Sweden.
! Jesper Lindkvist was funded by the Swedish National Space Board (SNSB project 201/15).
! @author    :  Jesper Lindkvist
! Email      :  jesper.lindkvist@umu.se
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


program picard
  use mpi
  implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1.0 Initialization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  logical :: writeParticles, writeFields, xflow, yflow, zflow, updateEfield, updateBfield
  integer :: Nx, Ny, Nz, iteration, iter_end, iter_start, ii, jj, kk, ss, num_local, num_global, &
             writeFields_iteration, writeParticles_iteration, maxit, writeOutput_interval
  real*8  :: Lx, Ly, Lz, Lx_min, Lx_max, Ly_min, Ly_max, Lz_min, Lz_max, & 
             Lx_local, Ly_local, Lz_local, xmin, xmax, ymin, ymax, zmin, zmax, dxyz(3), &
             writeParticles_interval, writeFields_interval, starttime, endtime, time, dt, maxerr, E0(4), B0(4), &
             uE_local, uE_global, uB_local, uB_global, uP_local, uP_global, &
             uEB0_global, uE0_global, uB0_global, uP0_global, u0_global, dV, ve(4)

! Macroparticle mass, charge and eta = charge/mass, temperature
! electrons per cell (ppc), number of particles per macroparticle (n2p)
  logical :: gyration(8), fluid(8), DEBUG
  integer :: ppc(8), ipp(8), hpp(8), epp(8)
  real*8  :: mass(8), charge(8), eta(8), kelvin(8), density(8), n2p(8), tB(4,8), sB(4,8), & 
             v0(4,8), vn(8), Qn(8), nu_i(8), nu_d(8), Ek(8)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1.1 DOMAIN PARAMETERS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ! should be multiple of 2
  integer, parameter :: iprocs = 2
  integer, parameter :: jprocs = 2
  integer, parameter :: kprocs = 2
  ! should be equal to INTEGER*8-2, e.g., 6, 30, 54, 78, 102
  integer, parameter :: Nx_local = 30
  integer, parameter :: Ny_local = 30
  integer, parameter :: Nz_local = 30
  ! number of species included in the simulation
  integer, parameter :: species = 3
  ! should be a multiple of 8
  integer, parameter :: max_per_proc = 2**24


! Physical and mathematical constants in SI units
  real*8, parameter :: tau  = 6.2831853071795864769252867665590d0  ! 2*pi
  real*8, parameter :: ec   = 1.60217646d-19  ! elementary charge [C]
  real*8, parameter :: mp   = 1.6726216d-27  ! proton mass [kg]
  real*8, parameter :: me   = 9.1093819d-31  ! electron mass [kg]
  real*8, parameter :: mm   = 1.8835311d-28  ! muon mass [kg]
  real*8, parameter :: mn   = 1.6749272d-27  ! neutron mass [kg]
  real*8, parameter :: bohr = 5.29177208d-11  ! Bohr radius [m]
  real*8, parameter :: c0   = 2.99792458d8  ! speed of light in vacuum [m/s]
  real*8, parameter :: mu0  = 1.25663706143591729538505735331180d-6  ! Magnetic constant [ SI ]
  real*8, parameter :: eps0 = 8.854187817d-12  ! Electric constant [ SI ]
  real*8, parameter :: kb   = 1.380650d-23  ! Bolzmann's constant [J/K]


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1.2 Domain Definitions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! real_particles
! real_particles: Holds the information of particles in the form
! real_particles(parameter_nr, local_particle_nr).
! parameter_nr:
! 1 - position x
! 2 - position y
! 3 - position z
! 4 - velocity x
! 5 - velocity y
! 6 - velocity z
! 7 - temp1
! 8 - temp2

! int_particles
! int_particles(local_particle_nr).
! hold the species number

  integer :: int_particles(max_per_proc)
  real*8 :: real_particles(8,max_per_proc)

! Inverse distance matrix
  real*8 :: D(2*Nx_local*iprocs+1,2*Ny_local*jprocs+1,2*Nz_local*kprocs+1)

! Electrostatic potential
  real*8 :: U(Nx_local+2,Ny_local+2,Nz_local+2)

! Electric field matrices
  real*8 :: E(3,Nx_local+2,Ny_local+2,Nz_local+2)

! Magnetic field matrices
  real*8 :: B(3,Nx_local+2,Ny_local+2,Nz_local+2)

! Number of particles for each species.
  real*8 :: N(species,Nx_local+2,Ny_local+2,Nz_local+2)

! Particle flux for each species.
  real*8 :: F(species,3,Nx_local+2,Ny_local+2,Nz_local+2)

! Charge in C for all species.
  real*8 :: C(Nx_local+2,Ny_local+2,Nz_local+2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2.0 Initialize MPI
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer :: myid, Nprocs, ierr, rc
  integer, allocatable :: seed(:)

  call MPI_init( ierr )
  call MPI_comm_rank( MPI_comm_world, myid, ierr )
  call MPI_comm_size( MPI_comm_world, Nprocs, ierr )

!  write (*,*) myid, Nprocs, iprocs*jprocs*kprocs

  if (Nprocs .ne. iprocs*jprocs*kprocs) then
    if (myid .eq. 0) then
      write(*,*) 'Wrong Number of Processes!'
!      call MPI_abort(MPI_comm_world, ierr)
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
  DEBUG = .false.  ! should we debug?
  writeParticles = .true.  ! should we write particles to disk?
  writeFields = .true.  ! should we write fields to disk?
  xflow = .true.  ! should we have non-periodic in x?
  yflow = .true.  ! should we have non-periodic in y?
  zflow = .true.  ! should we have non-periodic in z?
  updateEfield = .true.  ! should we update E-field?
  updateBfield = .false.  ! should we update B-field? NOT IMPLEMENTED

! INTEGER
  iter_end = 2**30  ! Number of iterations
  maxit = 2**20  ! maximum number of iterations in Poisson solver

! REAL
  xmin = -600.0d0  ! physical domain
  ymin = -600.0d0  ! physical domain
  zmin = -600.0d0  ! physical domain
  xmax =  600.0d0  ! physical domain
  ymax =  600.0d0  ! physical domain
  zmax =  600.0d0  ! physical domain
  dt   =  4.0d-6   ! time step
  writeParticles_interval = dt*1000  ! interval of particle writing to disk
  writeFields_interval = dt*100  ! interval of field writing to disk
  writeOutput_interval = 10  ! interval of default output writing
  maxerr = 2.0d0**(-16)  ! maximum error in Poisson solver

! ZERO B-FIELD FOR CHOSEN SPECIES. TRUE: F=q(E+vxB). FALSE: F=qE
  gyration(1)       = .true.      ! e-
  gyration(2)       = .true.      ! H+
  gyration(3)       = .true.      ! H2O+
  gyration(4)       = .true.     
  gyration(5)       = .true.     
  gyration(6)       = .true.
  gyration(7)       = .true.
  gyration(8)       = .true.

! FLUID APPROXIMATION. TRUE: v_perp=FxB/(q*B^2), a_para = F dot B... FALSE: Lorentz force
  fluid(1)          = .false.      ! e-
  fluid(2)          = .false.     ! H+
  fluid(3)          = .false.     ! H2O+
  fluid(4)          = .false.    
  fluid(5)          = .false.    
  fluid(6)          = .false.
  fluid(7)          = .false.
  fluid(8)          = .false.

! Micromass [kg].
  mass(1)           = me !*dsqrt(mp/me)      ! e-
  mass(2)           = mp*1.0d0      ! H+
  mass(3)           = mp*18.0d0     ! H2O+
  mass(4)           = mp*4.0d0      ! He++
  mass(5)           = mp*44.0d0     ! CO2+
  mass(6)           = mp*1.0d0      
  mass(7)           = mp*1.0d0    
  mass(8)           = mp*1.0d0      

! Microcharge [C]. Index 1 and 2 are electrons (-ec) and protons (+ec), which are calculated in Section 2.2.
  charge(3)         = ec         ! H2O+
  charge(4)         = ec         
  charge(5)         = ec         
  charge(6)         = ec
  charge(7)         = ec
  charge(8)         = ec

! UPSTREAM number density [m-3]. Index 1 are electrons and are calculated in Section 2.2.
  density(2)        = 1.07d6     ! 7d6 at 1 AU  ! H+
  density(3)        = 0.0d6      ! H2O+
  density(4)        = 0.0d6      
  density(5)        = 0.0d6      
  density(6)        = 0.0d6
  density(7)        = 0.0d6
  density(8)        = 0.0d6

! UPSTREAM number of particles per cell [#]. Index 1 are electrons and are calculated in Section 2.2.
  ppc(2)            = 64     ! H+
  ppc(3)            = 0     ! H2O+
  ppc(4)            = 0    
  ppc(5)            = 0    
  ppc(6)            = 0
  ppc(7)            = 0
  ppc(8)            = 0

! UPSTREAM temperature [K]  ( 1.1605d4 K = 1 eV )
  kelvin(1)         = 8.02d4    ! 15.0d4 at 1 AU.   ! e-
  kelvin(2)         = 4.27d4    ! 8.0d4 at 1 AU.    ! H+
  kelvin(3)         = 4.27d4    ! 8.0d4 at 1 AU.   ! H2O+
  kelvin(4)         = 0.0d4    
  kelvin(5)         = 0.0d4    
  kelvin(6)         = 0.0d4
  kelvin(7)         = 0.0d4
  kelvin(8)         = 0.0d4

! Set directed (bulk) velocity of plasma flow [m/s].
  v0(1,1) = -430.0d3 ! -430.0d0 at 1 AU
  v0(2,1) =  0.0d3
  v0(3,1) =  0.0d3

  v0(:,2) = 1.0d0*v0(:,1)
  v0(:,3) = 1.0d0*v0(:,1)
  v0(:,4) = 1.0d0*v0(:,1)
  v0(:,5) = 1.0d0*v0(:,1)
  v0(:,6) = 1.0d0*v0(:,1)
  v0(:,7) = 1.0d0*v0(:,1)
  v0(:,8) = 1.0d0*v0(:,1)

! Set magnetic field strength in x(1), y(2), z(3) in [T].
  B0(1) =  0.0d-9
  B0(2) =  400.0d0*1.78d-9   ! 6.0d-9 at 1 AU
  B0(3) =  0.0d-9

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2.1.1 HASER MODEL PARAMETERS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Set radial neutral velocity for the Haser model [m/s]
  vn(2) = 0.0d0     ! H+
  vn(3) = 629.0d0     ! H2O+
  vn(4) = 0.0d0     ! He++
  vn(5) = 0.0d0     ! CO2+
  vn(6) = 0.0d0
  vn(7) = 0.0d0
  vn(8) = 0.0d0

! Set production rate of H2O [#/s]
  Qn(2) = 0.0d0     ! H+
  Qn(3) = (400.0d0)**(-3)*1.99d26    !  2.59d28 *R^(-5.18) H2O+
  Qn(4) = 0.0d0     ! He++
  Qn(5) = 0.0d0     ! CO2+
  Qn(6) = 0.0d0
  Qn(7) = 0.0d0
  Qn(8) = 0.0d0

! Excess energy of end products [J]
  Ek(2) = 0.0d0  ! H
  Ek(3) = 1.0d0*1.52d1*ec   ! (1.24d1 for low activity) H2O
  Ek(4) = 0.0d0  ! He+
  Ek(5) = 0.0d0  ! CO2
  Ek(6) = 0.0d0
  Ek(7) = 0.0d0
  Ek(8) = 0.0d0

! Ionization rate of species 3 [1/s]
  nu_i(2) = 0.0d0  ! H
  nu_i(3) = 1.26d-7   ! 8.28d-7   H2O
  nu_i(4) = 0.0d0  ! He+
  nu_i(5) = 0.0d0  ! CO2
  nu_i(6) = 0.0d0
  nu_i(7) = 0.0d0
  nu_i(8) = 0.0d0

! Destruction rate of species 3 [1/s]
  nu_d(2) = 0.0d0  ! H
  nu_d(3) = 1.26d-7    ! 8.28d-7   H2O 
  nu_d(4) = 0.0d0  ! He+
  nu_d(5) = 0.0d0  ! CO2
  nu_d(6) = 0.0d0
  nu_d(7) = 0.0d0
  nu_d(8) = 0.0d0

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

! The electrons quasi-neutralizes

  vn(1) = 0.0d0  ! e-
  Qn(1) = 0.0d0  ! e-
  Ek(1) = 0.0d0  ! e-
  nu_i(1) = 0.0d0  ! e-
  nu_d(1) = 0.0d0  ! e-

  charge(1)         = -ec
  charge(2)         = ec
  density(2:)       = density(2:)*dble(ppc(2:))/dble(max(ppc(2:),1))
  density(1)        = sum( -density(2:)*charge(2:) )*charge(1)**(-1)
  do ss = 1, 8
    if (density(ss) .le. 0.0d0 .or. species .lt. ss) then
      ppc(ss) = 0
      density(ss) = 0.0d0
    end if
  end do
  density(1)        = sum( -density(2:)*charge(2:)/charge(1) )


  epp               = max( idint( 0.5d0 - charge/charge(1) ), 0 )
  hpp               = max( idint( 0.5d0 - charge/charge(2) ), 0 )
  hpp(1)            = 0
  
  ipp(3:)           = ppc(3:)
  ipp(2)            = ppc(2)-sum(ipp(3:)*hpp(3:))
  ipp(1)            = 0

  ppc(1)            = sum(ipp(2:)*epp(2:))

  n2p               = 0.0d0
  do ss = 1, 8
    if (ppc(ss) .gt. 0) then
      n2p(ss)       = dV*density(ss)*dble(ppc(ss))**(-1)
    else
      n2p(ss)       = n2p(1)
    end if
  end do


! q/m [#]
  eta               = charge/mass

! Velocity and B-field
  ve(1) = sum( ipp(:)*(n2p(:)*mass(:)+epp(:)*n2p(1)*mass(1)+hpp(:)*n2p(2)*mass(2))*v0(1,:) ) &
        / sum( ipp(:)*(n2p(:)*mass(:)+epp(:)*n2p(1)*mass(1)+hpp(:)*n2p(2)*mass(2)) )
  ve(2) = sum( ipp(:)*(n2p(:)*mass(:)+epp(:)*n2p(1)*mass(1)+hpp(:)*n2p(2)*mass(2))*v0(2,:) ) &
        / sum( ipp(:)*(n2p(:)*mass(:)+epp(:)*n2p(1)*mass(1)+hpp(:)*n2p(2)*mass(2)) )
  ve(3) = sum( ipp(:)*(n2p(:)*mass(:)+epp(:)*n2p(1)*mass(1)+hpp(:)*n2p(2)*mass(2))*v0(3,:) ) &
        / sum( ipp(:)*(n2p(:)*mass(:)+epp(:)*n2p(1)*mass(1)+hpp(:)*n2p(2)*mass(2)) )

  do ss= 1, 8
    v0(4,ss) = dsqrt( sum( v0(1:3,ss)**2 ) )
  end do
  ve(4)   = dsqrt( sum( ve(1:3)**2 ) )
  B0(4)   = dsqrt( sum( B0(1:3)**2 ) )
  B(1,:,:,:) = B0(1)
  B(2,:,:,:) = B0(2)
  B(3,:,:,:) = B0(3)

  E0(1) = ve(3)*B0(2)-ve(2)*B0(3)
  E0(2) = ve(1)*B0(3)-ve(3)*B0(1)
  E0(3) = ve(2)*B0(1)-ve(1)*B0(2)
  E0(4) = dsqrt( sum( E0(1:3)**2 ) )

! Initialize Distance matrix.
  call CalcDistance(D,Nx,Ny,Nz,dxyz,xflow,yflow,zflow)   


! Initialize fields and potentials to zero.
  U                 = 0.0d0
  N                 = 0.0d0
  F                 = 0.0d0
  E                 = 0.0d0
  real_particles    = 0.0d0
  int_particles     = 0


! Initialize scalars to zero.
  time = 0.0d0
  iter_start = 0
  writeFields_iteration = 0
  writeParticles_iteration = 0
  num_local = 0
  num_global = 0


! For the Boris algorithm from Birdsall and Langdon for constant homogeneous magnetic field
  do ss = 1, 8
    if ( gyration(ss) ) then
      tB(1,ss) = 0.5d0*dt*B0(1)*eta(ss)
      tB(2,ss) = 0.5d0*dt*B0(2)*eta(ss)
      tB(3,ss) = 0.5d0*dt*B0(3)*eta(ss)
      tB(4,ss) = 0.5d0*dt*B0(4)*eta(ss)
    else
      tB(:,ss) = 0.0d0
    end if
  end do

  sB(1,:) = 2.0d0*tB(1,:)/(1.0d0+tB(4,:)**2)
  sB(2,:) = 2.0d0*tB(2,:)/(1.0d0+tB(4,:)**2)
  sB(3,:) = 2.0d0*tB(3,:)/(1.0d0+tB(4,:)**2)
  sB(4,:) = 2.0d0*tB(4,:)/(1.0d0+tB(4,:)**2)


! Diagnostics
  uB_local      = 0.0d0
  uB_global     = 0.0d0
  uE_local      = 0.0d0
  uE_global     = 0.0d0
  uP_local      = 0.0d0
  uP_global     = 0.0d0
  uE0_global    = 0.5d0*( E0(4)**2*eps0 ) * Lx*Ly*Lz
  uB0_global    = 0.5d0*( B0(4)**2/mu0  ) * Lx*Ly*Lz
  uEB0_global   = uE0_global + uB0_global
  uP0_global    = 0.5d0*sum(  density(:)*( mass(:)*v0(4,:)**2 + 3.0d0*kb*kelvin(:) )  ) * Lx*Ly*Lz  ! Addition of both bulk and temperature (not integrated)
  u0_global     = uEB0_global + uP0_global


! Calculate plasma parameters and write warnings
  if (myid .eq. 0) then
     call system('mkdir save/')
     call system('rm -r save/*')
     call system('rmdir save/')
     call system('mkdir save/')
     call WriteWarnings(mass,charge,eta,kelvin,density,v0,ve,dt,n2p,dxyz,B0,Nx,Ny,Nz,xmin,xmax,ymin,ymax,zmin,zmax,&
                          iprocs,jprocs,kprocs,gyration,xflow,yflow,zflow,Qn,vn,Ek,nu_i,nu_d,ppc)
  end if
  call MPI_barrier( MPI_comm_world, ierr)

! Write particle status
  if (myid .eq. 0) then
    write(*,*) '****************************************************************'
    write(*,*) 'PARAMETER   SPECIES_1   SPECIES_2   SPECIES_3   SPECIES_4   SPECIES_5   SPECIES_6   SPECIES_7   SPECIES_8'
    write(*,*) 'n2p', n2p
    write(*,*) 'ppc', ppc
    write(*,*) 'ipp', ipp
    write(*,*) 'epp', epp
    write(*,*) 'hpp', hpp
    write(*,*) '****************************************************************'
  end if


! Initializing the particles in the domain and writing to disk

  call ParticlesInitialize(real_particles,int_particles,Lx_min,Ly_min,Lz_min,dxyz,kelvin,mass,v0,ve,tB,max_per_proc,&
                               num_local,ppc,epp,hpp,ipp,species,Nx_local,Ny_local,Nz_local)
     
  if (writeParticles) then
    if (myid .eq. 0) then
      call system('mkdir save/particles/')
    end if
    call MPI_barrier( MPI_comm_world, ierr)
    if (DEBUG) then
      call DumpParticles( real_particles, int_particles, num_local, max_per_proc, writeParticles_iteration, myid )
    else
      writeParticles_iteration = writeParticles_iteration + 1
    end if
  end if

  if (writeFields) then
    call ParticlesGridifyAll(F,N,C,real_particles,int_particles,Lx_min,Lx_max,Ly_min,Ly_max,Lz_min,Lz_max,dxyz,&
                   xmin,xmax,ymin,ymax,zmin,zmax,species,Nx_local,Ny_local,Nz_local,max_per_proc,num_local,&
                   xflow,yflow,zflow,iprocs,jprocs,kprocs,charge,n2p,myid)
    call BCpotentials(U,C,D,xmin,xmax,ymin,ymax,zmin,zmax,Nx_local,Ny_local,Nz_local,Nx,Ny,Nz,&
                        xflow,yflow,zflow,species,iprocs,jprocs,kprocs,myid)
    call CalcEpotential(U,C,Nx_local,Ny_local,Nz_local,Lx_min,Lx_max,Ly_min,Ly_max,Lz_min,Lz_max,xflow,yflow,zflow,&
                            xmin,xmax,ymin,ymax,zmin,zmax,dxyz,maxerr,maxit,iprocs,jprocs,kprocs,species,myid)
    call CalcEfield(E,U,Nx_local,Ny_local,Nz_local,dxyz,xflow,yflow,zflow,iprocs,jprocs,kprocs,myid)
    call BCfields(E,U,B0,E0,Nx_local,Ny_local,Nz_local,xflow,yflow,zflow,iprocs,jprocs,kprocs,myid)
    if (myid .eq. 0) then
      call system('mkdir save/fields/')
    end if
    call MPI_barrier( MPI_comm_world, ierr)
    if (DEBUG) then
      call DumpFields( U, N, F, E, Nx_local, Ny_local, Nz_local, writeFields_iteration, species, myid )
    else
      writeFields_iteration = writeFields_iteration + 1
    end if
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3.0 Start of Main loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (myid .eq. 0) then
    write(*,*) 'Engaging main loop!'
  end if

  do iteration = iter_start+1, iter_end

    time = time + dt   

    call ParticlesModifyExp( )

    call BCinflow(real_particles,int_particles,Lx_min,Ly_min,Lz_min,Lx_max,Ly_max,Lz_max,dxyz,kelvin,mass,v0,ve,&
                               max_per_proc,num_local,ppc,epp,hpp,ipp,species,Nx_local,Ny_local,Nz_local,&
                               xmin,xmax,ymin,ymax,zmin,zmax,tB,xflow,yflow,zflow)

    call AdvancePosition(real_particles, dt, max_per_proc, num_local) 

    call ParticlesProduce(real_particles,int_particles,Lx_min,Ly_min,Lz_min,dxyz,kelvin,mass,nu_d,nu_i,Ek,Qn,vn,v0,ve,&
                              max_per_proc,num_local,ppc,epp,hpp,ipp,species,Nx_local,Ny_local,Nz_local,dt,n2p) 

    call BCoutflow(real_particles,int_particles,Lx_min,Lx_max,Ly_min,Ly_max,Lz_min,Lz_max,&
                   xmin,xmax,ymin,ymax,zmin,zmax,max_per_proc,num_local,xflow,yflow,zflow)

    call ParticlesTransfer(real_particles,int_particles,Lx_min,Lx_max,Ly_min,Ly_max,Lz_min,Lz_max,&
                           xmin,xmax,ymin,ymax,zmin,zmax,max_per_proc,num_local,iprocs,jprocs,kprocs,myid)

    if (updateEfield) then

      call ParticlesGridify(N,C,real_particles,int_particles,Lx_min,Lx_max,Ly_min,Ly_max,Lz_min,Lz_max,dxyz,&
                   xmin,xmax,ymin,ymax,zmin,zmax,species,Nx_local,Ny_local,Nz_local,max_per_proc,num_local,&
                   xflow,yflow,zflow,iprocs,jprocs,kprocs,charge,n2p,myid)

      call BCpotentials(U,C,D,xmin,xmax,ymin,ymax,zmin,zmax,Nx_local,Ny_local,Nz_local,Nx,Ny,Nz,&
                        xflow,yflow,zflow,species,iprocs,jprocs,kprocs,myid)

      call CalcEpotential(U,C,Nx_local,Ny_local,Nz_local,Lx_min,Lx_max,Ly_min,Ly_max,Lz_min,Lz_max,xflow,yflow,zflow,&
                            xmin,xmax,ymin,ymax,zmin,zmax,dxyz,maxerr,maxit,iprocs,jprocs,kprocs,species,myid)

      call CalcEfield(E,U,Nx_local,Ny_local,Nz_local,dxyz,xflow,yflow,zflow,iprocs,jprocs,kprocs,myid)

      call BCfields(E,U,B0,E0,Nx_local,Ny_local,Nz_local,xflow,yflow,zflow,iprocs,jprocs,kprocs,myid)

    end if

    call ParticlesModifyImp( )



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3.1 Advance velocities and write to disk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    if (writeFields .or. writeParticles) then
      if ( ( dmod(time, writeFields_interval) .lt. dt ) .or. ( dmod(time, writeParticles_interval) .lt. dt ) ) then

        call AdvanceVelocityImp(real_particles,int_particles,E,B0,Lx_min,Ly_min,Lz_min,dxyz,0.5d0*dt,dt,v0,ve,eta,tB,&
                                num_local,max_per_proc,species,Nx_local,Ny_local,Nz_local,fluid,gyration)

        if ( (writeParticles) .and. ( dmod(time, writeParticles_interval) .lt. dt ) ) then
          call DumpParticles( real_particles, int_particles, num_local, max_per_proc, writeParticles_iteration, myid )
        end if

        if ( (writeFields) .and. ( dmod(time, writeFields_interval) .lt. dt ) ) then
          call ParticlesGridifyAll(F,N,C,real_particles,int_particles,Lx_min,Lx_max,Ly_min,Ly_max,Lz_min,Lz_max,dxyz,&
                   xmin,xmax,ymin,ymax,zmin,zmax,species,Nx_local,Ny_local,Nz_local,max_per_proc,num_local,&
                   xflow,yflow,zflow,iprocs,jprocs,kprocs,charge,n2p,myid)
          call DumpFields( U, N, F, E, Nx_local, Ny_local, Nz_local, writeFields_iteration, species, myid )
        end if

        call AdvanceVelocityExp(real_particles,int_particles,E,B0,Lx_min,Ly_min,Lz_min,dxyz,0.5d0*dt,dt,v0,ve,eta,tB,&
                                num_local,max_per_proc,species,Nx_local,Ny_local,Nz_local,fluid,gyration)

      else

        call AdvanceVelocity(real_particles,int_particles,E,B0,Lx_min,Ly_min,Lz_min,dxyz,dt,v0,ve,eta,tB,sB,&
                             num_local,max_per_proc,species,Nx_local,Ny_local,Nz_local,fluid,gyration)
    
      end if
    else

      call AdvanceVelocity(real_particles,int_particles,E,B0,Lx_min,Ly_min,Lz_min,dxyz,dt,v0,ve,eta,tB,sB,&
                           num_local,max_per_proc,species,Nx_local,Ny_local,Nz_local,fluid,gyration)

    end if
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3.2 Write information to default output
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (iteration .eq. 1) then
      if (myid .eq. 0) then
        write (*,*) 'The first iteration was a success!'
        write (*,*) ' '
        write (*,*) '****************************************************************'
        write (*,*) '   ITER           TIME   MPI_RANK  NUM_LOCAL NUM_GLOBAL    PART_ENERGY       E_ENERGY       B_ENERGY'
      end if
    end if


    if (mod(iteration,writeOutput_interval) .eq. 0 ) then

      call MPI_reduce( num_local, num_global, 1, MPI_INT, MPI_SUM, 0, MPI_comm_world, ierr )

      uB_local = sum(     B(:,2:Nx_local+1,2:Ny_local+1,2:Nz_local+1       )**2   )*0.5d0/mu0*dV
      call MPI_reduce( uB_local, uB_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_comm_world, ierr )
      
      uE_local = sum(    (E(1,2:Nx_local+1,2:Ny_local+1,2:Nz_local+1)+E0(1))**2 &
                       + (E(2,2:Nx_local+1,2:Ny_local+1,2:Nz_local+1)+E0(2))**2 &
                       + (E(3,2:Nx_local+1,2:Ny_local+1,2:Nz_local+1)+E0(3))**2   )*0.5d0*eps0*dV
      call MPI_reduce( uE_local, uE_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_comm_world, ierr )

      uP_local = sum(  n2p(int_particles(1:num_local))*mass(int_particles(1:num_local))*&
                       ( real_particles(4,1:num_local)**2 &
                       + real_particles(5,1:num_local)**2 &
                       + real_particles(6,1:num_local)**2 )  )*0.5d0
      call MPI_reduce( uP_local, uP_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_comm_world, ierr )

      if (myid .eq. 0) then
        write(*,fmt='(I8,a,E12.6,a,I8,a,I8,a,I8,a,E12.6,a,E12.6,a,E12.6)') &
          iteration,'   ',time,'   ',myid,'   ',num_local,'   ',num_global,'   ',uP_global/u0_global,'   ',&
          uE_global/u0_global,'   ',uB_global/u0_global
      end if
    end if

  end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 4.0 After Main loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call MPI_finalize(rc)

  stop
end program picard
