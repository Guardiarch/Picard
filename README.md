# Picard
3D Electrostatic Plasma Solver


This is a particle-in-cell plasma code 'picard' that was developed
by Jesper Lindkvist and Herbert Gunell with start in 2016 using
resources provided by the Swedish National Infrastructure for
Computing (SNIC) at the High Performance Computing Center North
(HPC2N), Umeå University, Sweden. Jesper Lindkvist was funded by
the Swedish National Space Board (SNSB project 201/15) and
Herbert Gunell by the Swedish National Space Agency (SNSA
project 108/18).
@author    :  Jesper Lindkvist
Email      :  jesper.lindkvist.se@gmail.com

@author    :  Herbert Gunell
Email      :  herbert.gunell@physics.org

The position and velocity of macroparticles are leap-frogged in time.
Their charge densities are deposited to the grid.
Poisson's equation is solved iteratively on the grid in a central
finite difference manner. The potential at all boundaries is calculated
in the frame of zero E-field (solar wind frame), assuming open boundary
conditions (no charge contribution from outside the domain, and zero
potential at infinity).
 
dirspicard         creates the necessary subdirectories in your working
		   directory
packthingspicard   moves and packs the output into .tgz files then
		   rebuilds the outp directory
packthingspicard   just packs things - nothing else
m_files		   is a directory containing matlab m-files for
		   postprocessing purposes. The basic procedure
		   is to start by using picardconvert.m to convert the
		   output text files to matlab format and then
		   to use the other.m-files for plotting.


		 
