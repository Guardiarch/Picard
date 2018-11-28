# Picard
3D Electrostatic Plasma Solver


This is a particle-in-cell plasma code 'picard' was developed by
Jesper Lindkvist with start in 2016 using resources provided
by the Swedish National Infrastructure for Computing (SNIC) at the
High Performance Computing Center North (HPC2N), Umeå University, Sweden.
Jesper Lindkvist was funded by the Swedish National Space Board
(SNSB project 201/15).
@author    :  Jesper Lindkvist
Email      :  jesper.lindkvist@umu.se

Some updates were made by Herbert Gunell in 2018.
Email      :  herbert.gunell@physics.org

The position and velocity of macroparticles are leap-frogged in time.
Their charge densities are deposited to the grid.
Poisson's equation is solved iteratively on the grid in a central
finite difference manner.
The potential at all boundaries is calculated in the frame of zero E-field
(solar wind frame), assuming open boundary conditions (no charge contribution
from outside the domain, and zero potential at infinity).
 
