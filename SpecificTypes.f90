module SpecificTypes

  type particlearrays
     ! species should be a vector of length max_per_proc, holding the
     ! species number of each particle.
     integer, allocatable :: species(:)
     ! coordinates should be a matrix of dimensions (6,max_per_proc),
     ! containing information about
     ! 1 - position x
     ! 2 - position y
     ! 3 - position z
     ! 4 - velocity x
     ! 5 - velocity y
     ! 6 - velocity z
     real*8,  allocatable :: coordinates(:,:)
     ! When used as a type for message buffers, max_per_proc above
     ! is replaced by the number of particles to be transported. 
  end type particlearrays

  type particlespecies
     ! parameters specific to each species
     ! read from input file
     logical gyration, fluid, cometion
     integer ppc ! particles per cell for this species in the initial plasma
     integer productspecies ! species number for the electron produced
                            ! when ionising this species
     ! v0 is the solar wind velocity
     real*8 mass, charge, upstreamdensity, upstreamkelvin, v0(4)
     ! For cometary ions and their parent neutrals we have the following
     ! Haser and Galand model parameters.
     ! vn: radial neutral velocity for the Haser model [m/s]
     ! Qn: neutral production rate [#/s]
     ! Ek: excess energy of end ionisation products [J]
     ! nu_i: ionization rate of this species [1/s]
     real*8 vn, Qn, Ek, nu_i
     ! calculated
     real*8 n2p, tB(4), sB(4), eta
     ! n2p - number of real particles per macro-particle
     ! eta - charge to mass ratio
     ! tB and sB are used for cross products in the Boris algorithm
     integer ipp, epp, hpp ! Jesper's variables for particle generation handling
  end type particlespecies

  type EandPprobe
   ! Electric field probe
   ! read from input file
     ! desired coordinates of the probe
     real*8 xc, yc, zc
     ! calculated
     ! mine is true is in the space handled by this process
     logical mine
     integer probeno ! probe number
     ! indeces of the probe position in order to get the field
     ! from the matrix within this process
     integer ic(3)
     ! real position of the probe
     real*8 rc(3)
     ! the radius of the probe for particle saving purposes
     real*8 rprobe
     ! The field that we wish to probe
     real*8, allocatable :: E(:,:)
     ! vector of the iterations at which we probe the field
     integer, allocatable :: iterations(:)
  end type EandPprobe

end module SpecificTypes
