# Picard
3D Electrostatic Plasma Solver

This is a particle-in-cell plasma code 'picard' was developed by Jesper Lindkvist with start in 2016 using resources provided
 by the Swedish National Infrastructure for Computing (SNIC) at the
 High Performance Computing Center North (HPC2N), Umeå University, Sweden.
Jesper Lindkvist was funded by the Swedish National Space Board (SNSB project 201/15).
@author    :  Jesper Lindkvist
Email      :  jesper.lindkvist@umu.se

The position and velocity of macroparticles are leap-frogged in time.
Their charge densities are deposited to the grid.
Poisson's equation is solved iteratively on the grid in a central finite difference manner.

The sections of the main program 'picard.f90' is structured as follows:

1.0 Initialization
    Simulation variables are declared here.
1.1 DOMAIN PARAMETERS
    Parameters that define the domain are declared here.
1.2 Domain Definitions
    Field arrays are declared.
2.0 Initialize MPI
    Mpi starts, and initializes a different random seed for each process.
2.1.0 SIMULATION PARAMETERS
    Define your physical world by changing the parameters.
2.1.1 HASER MODEL PARAMETERS
    For outflow from a source with neglible gravity.
2.2 Simulation Definitions
    Preparing the simulation to start. Quasi-neutrality is assumed when calculating new variables.
    Initializing particles initially in the domain, and writing to disk.
    The file 'parameters.dat' contains the simulation parameters used for this setup.
    The file 'warnings.dat' contains the warnings related to the physics of this setup.
3.0 Start of Main loop
    The main loop will be explained further down in the README.
3.1 Advance velocities and write to disk
    Boris algorithm for velocities and writing to disk.
3.2 Write information to default output
    Diagnostics
4.0 After Main loop
    Finalize mpi.

Sections that should be changed for a run:
1.1 DOMAIN PARAMETERS
2.1.0 SIMULATION PARAMETERS
2.1.1 HASER MODEL PARAMETERS
These can be improved by introducing a read from a parameter file instead.

Function breakdown of Section 3. MAIN LOOP

ParticlesModifyExp
    NOT USED. Add any monte carlo process here, and advance their position with a random dt=[0,dt/2].
BCinflow
    Particles at non-periodic boundaries are produced.
AdvancePosition
    The position of all macroparticles is advanced in time by dt.
ParticlesProduce
    Particles are produced inside the domain, and their velocities are advanced with a random dt=[0,dt].
BCoutflow
    Particles outside the domain are removed. So are also particles that would have to move two processes away.
ParticlesTransfer
    Paritcles are transferred to the process that contains their grid position.
ParticlesGridify
    Number densities are deposited to the grid.
BCpotentials
    The potential at all boundaries is calculated in the frame of zero E-field (solar wind frame),
    assuming open boundary conditions (no charge contribution from outside the domain, and zero potential at infinity).
CalcEpotential
    Iterative Poisson solver.
    The convergence criterion can be changed to your need. Might be needed if the solution never converges.
    See the source file for a suggestion!
CalcEfield
    Central difference negative gradient solver.
BCfields
    If any further boundary condition on the fields should be enforced, do so here!
ParticlesModifyImp
    NOT USED. Add any monte carlo process here, and advance their position with a random dt=[0,dt/2].

IF WRITING TO DISK AT THAT TIME STEP:
AdvanceVelocityImp
    Advance the velocity with the Boris algorithm, but stopping midway (half a dt)
DumpParticles
    Write particles to disk in the file-format:
    particles_XXXXX_YYYYY
    where XXXXX is the iteration number (how many times you have written particles to disk)
    and YYYYY is the process number (each process writes their own particles to disk).
    Each row in the file will contain
    POSX POSY POSZ VELX VELY VELZ SPECIES_NUMBER
ParticlesGridifyAll
    Deposits number densities and fluxes to the grid.
DumpFields
    Write fields to disk in the file-format:
    N1_XXXXX_YYYYY
    where XXXXX is the iteration number (how many times you have written fields to disk)
    and YYYYY is the process number (each process writes their own fields to disk).
    Each row contains a value which corresponds to a grid point starting with indeces x: i=1 ,y: j=1, z: k=1,
    which are first looped over 'i', then over 'j', then over 'k'.
    The different fields are:
    UE: electric potential in volts.
    Ex: electric field in x.
    Ey: electric field in y.
    Ez: electric field in z.
    N#: number density of species #.
    X#: flux of species #.
    Y#: flux of species #.
    Z#: flux of species #.
AdvanceVelocityExp
    Advance the velocity with the Boris algorithm, but starting midway (half a dt more)

IF NOT:
AdvanceVelocity
    Advance the velocity with the Boris algorithm


