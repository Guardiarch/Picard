'''
This is a simple tutorial that explains how to plot the particle data from 'picard'.

@author    :  Jesper Lindkvist
Email      :  jesper.lindkvist@umu.se
'''



import sys
import os
import shutil
import numpy as np
import sys
import math
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import gridspec
import gc
import time as timetime
from matplotlib import rc

# Function for opening and reading the "FLASH.PAR" file
def loadPar(parFile, scaName):
    scalar = 0.0;
    try:
        # Open the file to be read
        f = open(parFile, 'r');
    except:
        strErr = '** ERROR ** Could not open file ', parFile;
        raise NameError(strErr);
    try:
        for line in f:              #Find the line containing 'scaName'
            if scaName in line:
                sca_str = line;
                sca_str = sca_str.split('=');
                sca_str = sca_str[1].split();
                for i in range(len(sca_str)):         #test if string contains digits
                    try:
                        scalar = float(sca_str[i]);
                        break;
                    except:
                        if (sca_str[i] == 'T'):
                            scalar = 1.0
                        elif (sca_str[i] == 'F'):
                            scalar = 0.0
                        else:
                            continue
                break
    except:
        strErr = '** ERROR ** Could not access the scalar "', scaName , '"';
        f.close();
        raise NameError(strErr);
    # Close file
    f.close();
    return scalar;

# Function for opening and reading the "LIST" files
def loadTable(tableFile):
    i = 0;   
    column_numbers = 0;
    try:
        # Open the file to be read
        f = open(tableFile, 'r');
    except:
        strErr = '** ERROR ** Could not open file ', tableFile;
        raise NameError(strErr);
        
    # Getting the number of lines in table
    try:
        for i, line in enumerate(f):
            pass;
        temp = line.split();
        column_numbers = 8
        line_numbers=i+1;
    except:
        strErr = '** ERROR ** Could not count lines in table';
        f.close();
        raise NameError(strErr);
    table = np.zeros((line_numbers,column_numbers));
    f.close()
    
    # Setting values to the table
    i=0;
    try:
        # Open the file to be read
        f = open(tableFile, 'r');
    except:
        strErr = '** ERROR ** Could not open file ', tableFile;
        raise NameError(strErr);
        
    try:
        for i, line in enumerate(f):
            for j in range(7):
                start = j*(len(line)/7)
                stop = (j+1)*(len(line)/7)
                table[i,j] = float(line[start:stop])
#            table[i,6] = float(line[stop:])
    except:
        strErr = '** ERROR ** Could not read table';
        f.close();
        raise NameError(strErr);
        
    # Close file
    f.close();
    return table;


def getTheParticlesData(p, r, L, species):

    # Selected particle indexes.
    idx = np.where( (p[:, 0]-r[0])**2+(p[:, 1]-r[1])**2+(p[:, 2]-r[2])**2 < L**2*(1.0-np.abs(p[:, 6]-species)) )[0]
    # Number of selected particles
    N = idx.shape[0];
    # Allocate memory for the selected particles
    particle = np.zeros( (N, 8) );
    # If the number of particles are larger than zero
    if(N > 0):
        particle[:,0] = p[ idx, 0 ];      # rx
        particle[:,1] = p[ idx, 1 ];      # ry
        particle[:,2] = p[ idx, 2 ];      # rz
        
        particle[:,3] = p[ idx, 3 ];      # vx
        particle[:,4] = p[ idx, 4 ];      # vy
        particle[:,5] = p[ idx, 5 ];      # vz
        
        particle[:,6] = p[ idx, 6 ];      # species
        particle[:,7] = p[ idx, 7 ];      # dummy



    return (N, particle );
  
if __name__ == '__main__':

    t0 = timetime.time()
    matplotlib.rcParams['pdf.fonttype'] = 42
    gc.enable()

    try:
        moon = str(sys.argv[1])
    except:
        moon = 'comet'
    try:
        particles_iteration = str(sys.argv[2]);
    except:
        particles_iteration = ''
    try:
        version_number = str(sys.argv[3]);
    except:
        version_number = ''
    try:
        procname = str(sys.argv[4]);
    except:
        procname = ''
    
    root                = '/media/jesper/BANK_0001/picard/' + moon + '/' + moon + '_' + version_number
#    root                = os.path.dirname(os.path.abspath(__file__))

    save                = root + '/save'
    parFilename         = save + '/parameters.dat'
    resFilename         = save + '/warnings.dat'
    fields              = save + '/fields'
    particles           = save + '/particles'
    save_location       = save + '/plots_particles_' + particles_iteration   

    # Physical constants
    
    el_ch         = 1.60217657e-19;  # electrical charge
    mu_0          = 4.0e-7*np.pi;  # V*s/(A*m)
    kb            = 1.3806488e-23; # J/K
    amu           = 1.66053892e-27; # kg
    mass_electron = 0.000548579909;  # mass of electron [amu]
    eps_0         = 8.854188E-12;   # farads/meter
    c             = 299792458.0;      # Speed of light m/s
    
    
    # Load the parameters from FLASH.PAR
    
    n2p         = np.zeros(8)
    mass        = np.zeros(8)
    charge      = np.zeros(8)
    density     = np.zeros(8)
    kelvin      = np.zeros(8)
    Qn          = np.zeros(8)
    vn          = np.zeros(8)
    nu_i        = np.zeros(8)
    nu_d        = np.zeros(8)
    ppc         = np.zeros(8)
    gyration    = np.zeros(8)
    v0          = np.zeros(4)
    B0          = np.zeros(4)
    dxyz        = np.zeros(3)
    xmin        = 0.0
    xmax        = 0.0
    ymin        = 0.0
    ymax        = 0.0
    zmin        = 0.0
    zmax        = 0.0
    dt          = 0.0
    Nx          = 0.0
    Ny          = 0.0
    Nz          = 0.0
    iprocs      = 0.0
    jprocs      = 0.0
    kprocs      = 0.0
    
    
    for i in range(8):
        n2p[i]      = loadPar(parFilename, 'n2p_' + str(i+1))
        mass[i]     = loadPar(parFilename, 'mass_' + str(i+1))
        charge[i]   = loadPar(parFilename, 'charge_' + str(i+1))
        density[i]  = loadPar(parFilename, 'density_' + str(i+1))
        kelvin[i]   = loadPar(parFilename, 'kelvin_' + str(i+1))
        Qn[i]       = loadPar(parFilename, 'Qn_' + str(i+1))
        vn[i]       = loadPar(parFilename, 'vn_' + str(i+1))
        nu_i[i]     = loadPar(parFilename, 'nu_i_' + str(i+1))
        nu_d[i]     = loadPar(parFilename, 'nu_d_' + str(i+1))
        ppc[i]      = loadPar(parFilename, 'ppc_' + str(i+1))
        gyration[i] = loadPar(parFilename, 'gyration_' + str(i+1))
        if i<3:
            v0[i]   = loadPar(parFilename, 'v0_' + str(i+1))
            B0[i]   = loadPar(parFilename, 'B0_' + str(i+1))
            dxyz[i] = loadPar(parFilename, 'dxyz_' + str(i+1))
        if i<1:
            xmin    = loadPar(parFilename, 'xmin_' + str(i+1))
            xmax    = loadPar(parFilename, 'xmax_' + str(i+1))
            ymin    = loadPar(parFilename, 'ymin_' + str(i+1))
            ymax    = loadPar(parFilename, 'ymax_' + str(i+1))
            zmin    = loadPar(parFilename, 'zmin_' + str(i+1))
            zmax    = loadPar(parFilename, 'zmax_' + str(i+1))
            dt      = loadPar(parFilename, 'dt_' + str(i+1))
            Nx      = loadPar(parFilename, 'Nx_' + str(i+1))
            Ny      = loadPar(parFilename, 'Ny_' + str(i+1))
            Nz      = loadPar(parFilename, 'Nz_' + str(i+1))
            iprocs  = loadPar(parFilename, 'iprocs_' + str(i+1))
            jprocs  = loadPar(parFilename, 'jprocs_' + str(i+1))
            kprocs  = loadPar(parFilename, 'kprocs_' + str(i+1))

    B0[3] = np.sqrt(B0[0]**2+B0[1]**2+B0[2]**2)
    v0[3] = np.sqrt(v0[0]**2+v0[1]**2+v0[2]**2)

    density0    = np.sum(density[:]*np.sign(ppc[:]))
    mass0       = np.sum(mass[:]*density[:]*np.sign(ppc[:]))/density0
    charge0     = np.sum(charge[:]*density[:]*np.sign(ppc[:]))/density0
    r0          = np.array([0.0e3, 0.0e3, 0.0e3])
    
    print 'Mean density [m^-3]: ', density0, density
    print 'Mean mass [amu]: ', mass0/amu, mass/amu
    print 'Mean charge [e]: ', charge0/el_ch, charge/el_ch

    save_particles = save + '/particles_' + particles_iteration + '_' + str(int(round(r0[0]/1000.0))) + '_' + str(int(round(r0[1]/1000.0))) + '_' + str(int(round(r0[2]/1000.0))) + '.npz'

    try: 
        vars_npz            = np.load( save_particles )     
        N                   = vars_npz['N']  
        N0                  = vars_npz['N0']  
        N1                  = vars_npz['N1']  
        N2                  = vars_npz['N2']  
        p0                  = vars_npz['p0']  
        p1                  = vars_npz['p1']  
        p2                  = vars_npz['p2']
        vars_npz.close()
    except:
        # Checking what lists you have in the directory you ran the script.    
        filename = []
        for files in os.listdir(particles):
            if files.endswith(procname) and files.startswith("particles_" + particles_iteration):
                filename.append(particles + '/' + files);
        filename.sort()

        N = 0
        N0 = 0
        N1 = 0
        N2 = 0

        for i in range(len(filename[:])):    
          
            ptable = loadTable(filename[i])

            N_temp           = ptable.shape[0]
            N0_temp, p0_temp = getTheParticlesData(ptable, r0, dxyz[0], 1)
            N1_temp, p1_temp = getTheParticlesData(ptable, r0, dxyz[0], 2)
            N2_temp, p2_temp = getTheParticlesData(ptable, r0, dxyz[0], 3)

            N   = N + N_temp
            N0  = N0 + N0_temp
            N1  = N1 + N1_temp
            N2  = N2 + N2_temp
            
            if i == 0:
                p0 = 1.0*p0_temp
                p1 = 1.0*p1_temp
                p2 = 1.0*p2_temp
            else:
                if N0 > 0:
                    p0_temptemp = np.zeros((N0,8))
                    p0_temptemp[0:N0_temp,:] = 1.0*p0_temp
                    p0_temptemp[N0_temp:,:] = 1.0*p0
                    p0 = 1.0*p0_temptemp
                if N1 > 0:
                    p1_temptemp = np.zeros((N1,8))
                    p1_temptemp[0:N1_temp,:] = 1.0*p1_temp
                    p1_temptemp[N1_temp:,:] = 1.0*p1
                    p1 = 1.0*p1_temptemp
                if N2 > 0:
                    p2_temptemp = np.zeros((N2,8))
                    p2_temptemp[0:N2_temp,:] = 1.0*p2_temp
                    p2_temptemp[N2_temp:,:] = 1.0*p2
                    p2 = 1.0*p2_temptemp

        np.savez(save_particles
            ,   N                   = N
            ,   N0                  = N0
            ,   N1                  = N1
            ,   N2                  = N2
            ,   p0                  = p0
            ,   p1                  = p1
            ,   p2                  = p2
            )


    print N, N0, N1, N2

    dV          = dxyz[0]**3*4.0*np.pi/3.0
    deg         = 5

    lat = np.linspace(-90,90,num=np.round(180.0/(1.0*deg))+1,endpoint=True)
    lon = np.linspace(-180,180,num=np.round(360.0/(1.0*deg))+1,endpoint=True)

    dcosrad     = np.cos(np.deg2rad(deg/np.sqrt(2)))
    dOmega      = 2.0*np.pi*(1.0-dcosrad)
    u0          = 1.0*v0[3]
    vth         = np.sqrt( 2.0*kb*kelvin[1]/mass[1] )

    vx0         = 0.5*u0*(1.0+math.erf(u0/vth)) + 0.5*vth/np.sqrt(np.pi)*np.exp(-(u0/vth)**2)
    vy0         = 0.5*vth/np.sqrt(np.pi)
    vxm         = u0 + vth/np.sqrt(np.pi)*np.exp(-(u0/vth)**2)*(1.0+math.erf(u0/vth))**(-1)
    vym         = vth/np.sqrt(np.pi)
    dcosrad_0   = vx0/(vx0+vy0)
    dOmega_0    = 2.0*np.pi*(1.0-dcosrad_0)


    x               = np.zeros( (lon.shape[0],lat.shape[0]) )
    y               = 0.0*x
    z               = 0.0*x
    flux0_dOmega    = 0.0*x
    Eflux0_dOmega   = 0.0*x
    flux1_dOmega    = 0.0*x
    Eflux1_dOmega   = 0.0*x
    flux2_dOmega    = 0.0*x
    Eflux2_dOmega   = 0.0*x

    momentum        = np.zeros((8,4))
    momentum[0,0]   = np.sum( p0[:,3] )/N0
    momentum[0,1]   = np.sum( p0[:,4] )/N0
    momentum[0,2]   = np.sum( p0[:,5] )/N0
    momentum[1,0]   = np.sum( p1[:,3] )/N1
    momentum[1,1]   = np.sum( p1[:,4] )/N1
    momentum[1,2]   = np.sum( p1[:,5] )/N1
    momentum[2,0]   = np.sum( p2[:,3] )/N2
    momentum[2,1]   = np.sum( p2[:,4] )/N2
    momentum[2,2]   = np.sum( p2[:,5] )/N2
    momentum[:,3]   = np.sqrt( momentum[:,0]**2+momentum[:,1]**2+momentum[:,2]**2 )

    print 'Momentum electrons', momentum[0,:]
    print 'Momentum protons', momentum[1,:]
    print 'Momentum H2O+', momentum[2,:]

    for i in range(lon.shape[0]):
        for j in range(lat.shape[0]):
            if (lat[j]>=90.0 or lat[j]<=-90.0):
                x[i,j]  =  0.0
                y[i,j]  =  1.0*np.sign(lat[j])
                z[i,j]  =  0.0
            elif lat[j]==0.0:
                x[i,j]  =  1.0*np.cos(np.deg2rad(lon[i]))
                y[i,j]  =  0.0
                z[i,j]  =  1.0*np.sin(np.deg2rad(lon[i]))
            else:
                x[i,j]  =  1.0*np.cos(np.deg2rad(lat[j]))*np.cos(np.deg2rad(lon[i]))
                y[i,j]  =  1.0*np.sin(np.deg2rad(lat[j]))
                z[i,j]  =  1.0*np.cos(np.deg2rad(lat[j]))*np.sin(np.deg2rad(lon[i]))
            
    
    for i in range(lon.shape[0]):
        for j in range(lat.shape[0]):
            for p in range(p0.shape[0]):
                if (-(p0[p,3]*x[i,j]+p0[p,4]*y[i,j]+p0[p,5]*z[i,j]) >= np.sqrt(p0[p,3]**2+p0[p,4]**2+p0[p,5]**2)*dcosrad):
                    flux0_dOmega[i,j]   = flux0_dOmega[i,j]  + np.sqrt(p0[p,3]**2+p0[p,4]**2+p0[p,5]**2)
                    Eflux0_dOmega[i,j]  = Eflux0_dOmega[i,j] + np.sqrt(p0[p,3]**2+p0[p,4]**2+p0[p,5]**2)**3
            for p in range(p1.shape[0]):
                if (-(p1[p,3]*x[i,j]+p1[p,4]*y[i,j]+p1[p,5]*z[i,j]) >= np.sqrt(p1[p,3]**2+p1[p,4]**2+p1[p,5]**2)*dcosrad):
                    flux1_dOmega[i,j]   = flux1_dOmega[i,j]  + np.sqrt(p1[p,3]**2+p1[p,4]**2+p1[p,5]**2)
                    Eflux1_dOmega[i,j]  = Eflux1_dOmega[i,j] + np.sqrt(p1[p,3]**2+p1[p,4]**2+p1[p,5]**2)**3
            for p in range(p2.shape[0]):
                if (-(p2[p,3]*x[i,j]+p2[p,4]*y[i,j]+p2[p,5]*z[i,j]) >= np.sqrt(p2[p,3]**2+p2[p,4]**2+p2[p,5]**2)*dcosrad):
                    flux2_dOmega[i,j]   = flux2_dOmega[i,j]  + np.sqrt(p2[p,3]**2+p2[p,4]**2+p2[p,5]**2)
                    Eflux2_dOmega[i,j]  = Eflux2_dOmega[i,j] + np.sqrt(p2[p,3]**2+p2[p,4]**2+p2[p,5]**2)**3

    fluxmin_dOmega  = density[1]*vy0/np.pi
    fluxmax_dOmega  = density[1]*vx0/min(dOmega, np.pi)
    Efluxmin_dOmega = density[1]*vy0/np.pi*vx0**2*mass[1]*0.5
    Efluxmax_dOmega = density[1]*vx0/min(dOmega, np.pi)*vx0**2*mass[1]*0.5

    minfluxmin_dOmega  = fluxmin_dOmega*fluxmin_dOmega/fluxmax_dOmega
    minfluxmax_dOmega  = fluxmax_dOmega*fluxmin_dOmega/fluxmax_dOmega
    minEfluxmin_dOmega = Efluxmin_dOmega*Efluxmin_dOmega/Efluxmax_dOmega
    minEfluxmax_dOmega = Efluxmax_dOmega*Efluxmin_dOmega/Efluxmax_dOmega

    flux0_dOmega    = flux0_dOmega*n2p[0]/dV/dOmega
    flux1_dOmega    = flux1_dOmega*n2p[1]/dV/dOmega
    flux2_dOmega    = flux2_dOmega*n2p[2]/dV/dOmega
    Eflux0_dOmega   = Eflux0_dOmega*n2p[0]/dV/dOmega*mass[0]*0.5
    Eflux1_dOmega   = Eflux1_dOmega*n2p[1]/dV/dOmega*mass[1]*0.5
    Eflux2_dOmega   = Eflux2_dOmega*n2p[2]/dV/dOmega*mass[2]*0.5

    print deg, dOmega_0*180.0/np.pi
    print fluxmin_dOmega, fluxmax_dOmega, fluxmax_dOmega/fluxmin_dOmega
    print Efluxmin_dOmega, Efluxmax_dOmega, Efluxmax_dOmega/Efluxmin_dOmega
    print np.max(flux0_dOmega), np.max(flux1_dOmega), np.max(flux2_dOmega)

    thecmap = cm.jet
    theextent = np.array([-180-0.5*deg, 180+0.5*deg, -90-0.5*deg, 90+0.5*deg])
    my_xticks = np.linspace(-180,180,num=13,endpoint=True)
    my_yticks = np.linspace(-90,90,num=7,endpoint=True)
    my_cticks = np.linspace(0,10**13,num=11,endpoint=True)


    flux0_dOmega[flux0_dOmega < fluxmin_dOmega] = fluxmin_dOmega
    flux1_dOmega[flux1_dOmega < fluxmin_dOmega] = fluxmin_dOmega
    flux2_dOmega[flux2_dOmega < minfluxmin_dOmega] = minfluxmin_dOmega
    Eflux0_dOmega[Eflux0_dOmega < Efluxmin_dOmega] = Efluxmin_dOmega
    Eflux1_dOmega[Eflux1_dOmega < Efluxmin_dOmega] = Efluxmin_dOmega
    Eflux2_dOmega[Eflux2_dOmega < minEfluxmin_dOmega] = minEfluxmin_dOmega

    try:
        os.stat(save_location);
        shutil.rmtree(save_location);
    except:
        pass
    os.mkdir(save_location)


    cmap_title = r' '

    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 5.0))

    axes = plt.subplot(gs[1,:])
    pc = imshow( flux0_dOmega.T, aspect='equal', interpolation="nearest", norm=matplotlib.colors.LogNorm(vmin=fluxmin_dOmega, vmax=fluxmax_dOmega), origin='lower', extent=theextent, cmap=thecmap )
    ax = plt.gca()
    xlabel(r'$\mathrm{Azimuth} \; [ \mathrm{deg} ]$', fontsize=18)
    ylabel(r'$\mathrm{Elevation} \; [ \mathrm{deg} ]$', fontsize=18)
    plt.xticks(my_xticks, fontsize=12)
    plt.yticks(my_yticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{e^- differential \, flux} \; [ \mathrm{m^{-2} \, s^{-1} \, sr^{-1}} ]$', fontsize=18)

    axes = plt.subplot(gs[0,:])
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
#    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/flux0_' + particles_iteration + '.'   
    saveloc_png = saveloc + 'png'
    saveloc_eps = saveloc + 'eps'
    if saveloc != 'none':
        try:
            savefig(saveloc_png);
        except:
            print saveloc_png
            strErr = '** ERROR ** Could not save file';
            raise NameError(strErr);
        try:
            savefig(saveloc_eps);
        except:
            print saveloc_eps
            strErr = '** ERROR ** Could not save file';
            raise NameError(strErr);
    else:
        strErr = '** ERROR ** Save location is not specified'
        raise NameError(strErr);

    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 5.0))

    axes = plt.subplot(gs[1,:])
    pc = imshow( flux1_dOmega.T, aspect='equal', interpolation="nearest", norm=matplotlib.colors.LogNorm(vmin=fluxmin_dOmega, vmax=fluxmax_dOmega), origin='lower', extent=theextent, cmap=thecmap )
    ax = plt.gca()
    xlabel(r'$\mathrm{Azimuth} \; [ \mathrm{deg} ]$', fontsize=18)
    ylabel(r'$\mathrm{Elevation} \; [ \mathrm{deg} ]$', fontsize=18)
    plt.xticks(my_xticks, fontsize=12)
    plt.yticks(my_yticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{H^+ differential \, flux} \; [ \mathrm{m^{-2} \, s^{-1} \, sr^{-1}} ]$', fontsize=18)

    axes = plt.subplot(gs[0,:])
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
#    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/flux1_' + particles_iteration + '.'   
    saveloc_png = saveloc + 'png'
    saveloc_eps = saveloc + 'eps'
    if saveloc != 'none':
        try:
            savefig(saveloc_png);
        except:
            print saveloc_png
            strErr = '** ERROR ** Could not save file';
            raise NameError(strErr);
        try:
            savefig(saveloc_eps);
        except:
            print saveloc_eps
            strErr = '** ERROR ** Could not save file';
            raise NameError(strErr);
    else:
        strErr = '** ERROR ** Save location is not specified'
        raise NameError(strErr);

    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 5.0))

    axes = plt.subplot(gs[1,:])
    pc = imshow( flux2_dOmega.T, aspect='equal', interpolation="nearest", norm=matplotlib.colors.LogNorm(vmin=minfluxmin_dOmega, vmax=minfluxmax_dOmega), origin='lower', extent=theextent, cmap=thecmap )
    ax = plt.gca()
    xlabel(r'$\mathrm{Azimuth} \; [ \mathrm{deg} ]$', fontsize=18)
    ylabel(r'$\mathrm{Elevation} \; [ \mathrm{deg} ]$', fontsize=18)
    plt.xticks(my_xticks, fontsize=12)
    plt.yticks(my_yticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{H_2O^+ differential \, flux} \; [ \mathrm{m^{-2} \, s^{-1} \, sr^{-1}} ]$', fontsize=18)

    axes = plt.subplot(gs[0,:])
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
#    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/flux2_' + particles_iteration + '.'    
    saveloc_png = saveloc + 'png'
    saveloc_eps = saveloc + 'eps'
    if saveloc != 'none':
        try:
            savefig(saveloc_png);
        except:
            print saveloc_png
            strErr = '** ERROR ** Could not save file';
            raise NameError(strErr);
        try:
            savefig(saveloc_eps);
        except:
            print saveloc_eps
            strErr = '** ERROR ** Could not save file';
            raise NameError(strErr);
    else:
        strErr = '** ERROR ** Save location is not specified'
        raise NameError(strErr);


    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 5.0))

    axes = plt.subplot(gs[1,:])
    pc = imshow( Eflux0_dOmega.T/el_ch, aspect='equal', interpolation="nearest", norm=matplotlib.colors.LogNorm(vmin=Efluxmin_dOmega/el_ch, vmax=Efluxmax_dOmega/el_ch), origin='lower', extent=theextent, cmap=thecmap )
    ax = plt.gca()
    xlabel(r'$\mathrm{Azimuth} \; [ \mathrm{deg} ]$', fontsize=18)
    ylabel(r'$\mathrm{Elevation} \; [ \mathrm{deg} ]$', fontsize=18)
    plt.xticks(my_xticks, fontsize=12)
    plt.yticks(my_yticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{e^- differential \, Eflux} \; [ \mathrm{eV \, m^{-2} \, s^{-1} \, sr^{-1}} ]$', fontsize=18)

    axes = plt.subplot(gs[0,:])
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
#    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/Eflux0_' + particles_iteration + '.'   
    saveloc_png = saveloc + 'png'
    saveloc_eps = saveloc + 'eps'
    if saveloc != 'none':
        try:
            savefig(saveloc_png);
        except:
            print saveloc_png
            strErr = '** ERROR ** Could not save file';
            raise NameError(strErr);
        try:
            savefig(saveloc_eps);
        except:
            print saveloc_eps
            strErr = '** ERROR ** Could not save file';
            raise NameError(strErr);
    else:
        strErr = '** ERROR ** Save location is not specified'
        raise NameError(strErr);

    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 5.0))

    axes = plt.subplot(gs[1,:])
    pc = imshow( Eflux1_dOmega.T/el_ch, aspect='equal', interpolation="nearest", norm=matplotlib.colors.LogNorm(vmin=Efluxmin_dOmega/el_ch, vmax=Efluxmax_dOmega/el_ch), origin='lower', extent=theextent, cmap=thecmap )
    ax = plt.gca()
    xlabel(r'$\mathrm{Azimuth} \; [ \mathrm{deg} ]$', fontsize=18)
    ylabel(r'$\mathrm{Elevation} \; [ \mathrm{deg} ]$', fontsize=18)
    plt.xticks(my_xticks, fontsize=12)
    plt.yticks(my_yticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{H^+ differential \, Eflux} \; [ \mathrm{eV \, m^{-2} \, s^{-1} \, sr^{-1}} ]$', fontsize=18)

    axes = plt.subplot(gs[0,:])
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
#    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/Eflux1_' + particles_iteration + '.'    
    saveloc_png = saveloc + 'png'
    saveloc_eps = saveloc + 'eps'
    if saveloc != 'none':
        try:
            savefig(saveloc_png);
        except:
            print saveloc_png
            strErr = '** ERROR ** Could not save file';
            raise NameError(strErr);
        try:
            savefig(saveloc_eps);
        except:
            print saveloc_eps
            strErr = '** ERROR ** Could not save file';
            raise NameError(strErr);
    else:
        strErr = '** ERROR ** Save location is not specified'
        raise NameError(strErr);

    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 5.0))

    axes = plt.subplot(gs[1,:])
    pc = imshow( Eflux2_dOmega.T/el_ch, aspect='equal', interpolation="nearest", norm=matplotlib.colors.LogNorm(vmin=minEfluxmin_dOmega/el_ch, vmax=minEfluxmax_dOmega/el_ch), origin='lower', extent=theextent, cmap=thecmap )
    ax = plt.gca()
    xlabel(r'$\mathrm{Azimuth} \; [ \mathrm{deg} ]$', fontsize=18)
    ylabel(r'$\mathrm{Elevation} \; [ \mathrm{deg} ]$', fontsize=18)
    plt.xticks(my_xticks, fontsize=12)
    plt.yticks(my_yticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{H_2O^+ differential \, Eflux} \; [ \mathrm{eV \, m^{-2} \, s^{-1} \, sr^{-1}} ]$', fontsize=18)

    axes = plt.subplot(gs[0,:])
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
#    cbar.ax.minorticks_on()



    rc('text', usetex=False)
    saveloc   = save_location + '/Eflux2_' + particles_iteration + '.'  
    saveloc_png = saveloc + 'png'
    saveloc_eps = saveloc + 'eps'
    if saveloc != 'none':
        try:
            savefig(saveloc_png);
        except:
            print saveloc_png
            strErr = '** ERROR ** Could not save file';
            raise NameError(strErr);
        try:
            savefig(saveloc_eps);
        except:
            print saveloc_eps
            strErr = '** ERROR ** Could not save file';
            raise NameError(strErr);
    else:
        strErr = '** ERROR ** Save location is not specified'
        raise NameError(strErr);

    
#    show()
    cla()
    clf()
    close('all')
    rc('text', usetex=False)
    print 'END OF PROGRAM';
