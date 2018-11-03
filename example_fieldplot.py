'''
This is a simple tutorial that explains how to plot the field data from 'picard'.

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
        column_numbers = 1
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
            table[i,0] = float(line)
    except:
        strErr = '** ERROR ** Could not read table';
        f.close();
        raise NameError(strErr);
        
    # Close file
    f.close();
    return table;


def getMyParticlesData(p, r, L, species):
    # Selected particle indexes.
    idx = np.where( (p[:, 0]-r[0])**2+(p[:, 1]-r[1])**2+(p[:, 2]-r[2])**2 <= L**2 )[0]
                 
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
        fields_iteration = str(sys.argv[2])
    except:
        fields_iteration = ''
    try:
        version_number = str(sys.argv[3]);
    except:
        version_number = ''
    try:
        procname = str(sys.argv[4]);
    except:
        procname = ''
    
    root                = '/home/herbert/picard/a1'
#    root                = os.path.dirname(os.path.abspath(__file__))

    save                = root + '/save'
    parFilename         = save + '/parameters.dat'
    resFilename         = save + '/warnings.dat'
    fields              = save + '/fields'
    particles           = save + '/particles'
    save_location       = save + '/plots_fields_' + fields_iteration  
    save_fields         = save + '/fields_' + fields_iteration + '.npz'
    
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
    Nx          = 0
    Ny          = 0
    Nz          = 0
    iprocs      = 0
    jprocs      = 0
    kprocs      = 0
    
    
    for i in range(8):
        n2p[i]      = loadPar(parFilename, 'n2p_' + str(i+1))
        mass[i]     = loadPar(parFilename, 'mass_' + str(i+1))
        charge[i]   = loadPar(parFilename, 'charge_' + str(i+1))
        density[i]  = loadPar(parFilename, 'density_' + str(i+1))
        kelvin[i]   = loadPar(parFilename, 'kelvin_' + str(i+1))
        Qn[i]       = loadPar(parFilename, 'Qn_' + str(i+1))
        vn[i]       = loadPar(parFilename, 'vn_' + str(i+1))
        nu_i[i]     = loadPar(parFilename, 'nui_' + str(i+1))
        nu_d[i]     = loadPar(parFilename, 'nud_' + str(i+1))
        ppc[i]      = np.int_(loadPar(parFilename, 'ppc_' + str(i+1)))
        gyration[i] = np.int_(loadPar(parFilename, 'gyration_' + str(i+1)))
        if i<3:
            v0[i]   = loadPar(parFilename, 'v0_' + str(i+1))
            B0[i]   = loadPar(parFilename, 'B0_' + str(i+1))
            dxyz[i] = loadPar(parFilename, 'dxyz_' + str(i+1))
        if i<1:
            xmin    = loadPar(parFilename, 'xmin')
            xmax    = loadPar(parFilename, 'xmax')
            ymin    = loadPar(parFilename, 'ymin')
            ymax    = loadPar(parFilename, 'ymax')
            zmin    = loadPar(parFilename, 'zmin')
            zmax    = loadPar(parFilename, 'zmax')
            dt      = loadPar(parFilename, 'dt')
            Nx      = np.int_(loadPar(parFilename, 'Nx'))
            Ny      = np.int_(loadPar(parFilename, 'Ny'))
            Nz      = np.int_(loadPar(parFilename, 'Nz'))
            iprocs  = np.int_(loadPar(parFilename, 'iprocs'))
            jprocs  = np.int_(loadPar(parFilename, 'jprocs'))
            kprocs  = np.int_(loadPar(parFilename, 'kprocs'))
            xflow   = np.int_(loadPar(parFilename, 'xflow'))
            yflow   = np.int_(loadPar(parFilename, 'yflow'))
            zflow   = np.int_(loadPar(parFilename, 'zflow'))

    B0[3] = np.sqrt(B0[0]**2+B0[1]**2+B0[2]**2)
    v0[3] = np.sqrt(v0[0]**2+v0[1]**2+v0[2]**2)

    print Nx, Ny, Nz
    density0    = np.sum(density[:]*np.sign(ppc[:]))
    mass0       = np.sum(mass[:]*density[:]*np.sign(ppc[:]))/density0
    charge0     = np.sum(charge[:]*density[:]*np.sign(ppc[:]))/density0
    E0          = np.zeros(4)
    E0[0]       = -(v0[1]*B0[2]-v0[2]*B0[1])
    E0[1]       = -(v0[2]*B0[0]-v0[0]*B0[2])
    E0[2]       = -(v0[0]*B0[1]-v0[1]*B0[0])
    E0[3]       = np.sqrt(E0[0]**2+E0[1]**2+E0[2]**2)
    Emax        = np.sqrt(E0[3]**2 + (c*B0[3])**2)
    J0          = np.zeros(4)
    J0[0:3]     = charge[0]*B0[3]**2/(mu_0*mass[0]*v0[3]) * v0[0:3]/v0[3]  * (c/v0[3])**(-1)
    J0[3]       = np.sqrt(J0[0]**2+J0[1]**2+J0[2]**2)
    UE0         = np.abs(mass[0]*v0[3]**2/charge[0]) * (c/v0[3])
    r0          = np.array([0.0e3, 0.0e3, 0.0e3])


    
    print 'Running for', moon, 'version', version_number, 'iteration', fields_iteration
    print 'Mean density [m^-3]: ', density0, density
    print 'Mean mass [amu]: ', mass0/amu, mass/amu
    print 'Mean charge [e]: ', charge0/el_ch, charge/el_ch

    x   = np.zeros((Nx+2, Ny+2, Nz+2))
    y   = 0.0*x
    z   = 0.0*x
    UE  = 0.0*x
    Ex  = 0.0*x
    Ey  = 0.0*x
    Ez  = 0.0*x
    N0  = 0.0*x   # Neutrals
    N1  = 0.0*x
    N2  = 0.0*x
    N3  = 0.0*x
    X1  = 0.0*x
    X2  = 0.0*x
    X3  = 0.0*x
    Y1  = 0.0*x
    Y2  = 0.0*x
    Y3  = 0.0*x
    Z1  = 0.0*x
    Z2  = 0.0*x
    Z3  = 0.0*x





    try: 
        vars_npz            = np.load( save_fields ) 
        UE  = vars_npz['UE']  
        Ex  = vars_npz['Ex']  
        Ey  = vars_npz['Ey']  
        Ez  = vars_npz['Ez']  
        N1  = vars_npz['N1']  
        N2  = vars_npz['N2']  
        N3  = vars_npz['N3']  
        X1  = vars_npz['X1']  
        X2  = vars_npz['X2']  
        X3  = vars_npz['X3']  
        Y1  = vars_npz['Y1']  
        Y2  = vars_npz['Y2']  
        Y3  = vars_npz['Y3']  
        Z1  = vars_npz['Z1']  
        Z2  = vars_npz['Z2']  
        Z3  = vars_npz['Z3']  
        vars_npz.close()
        print 'Loading Data SMALL'
    except:    
        print 'Loading Data LARGE'
        # U
        filename_UE = []
        for files in os.listdir(fields):
            if files.endswith(procname) and files.startswith("UE_" + fields_iteration):
                filename_UE.append(fields + '/' + files);
        filename_UE.sort()

        # Ex
        filename_Ex = []
        for files in os.listdir(fields):
            if files.endswith(procname) and files.startswith("Ex_" + fields_iteration):
                filename_Ex.append(fields + '/' + files);
        filename_Ex.sort()

        # Ey
        filename_Ey = []
        for files in os.listdir(fields):
            if files.endswith(procname) and files.startswith("Ey_" + fields_iteration):
                filename_Ey.append(fields + '/' + files);
        filename_Ey.sort()

        # Ez
        filename_Ez = []
        for files in os.listdir(fields):
            if files.endswith(procname) and files.startswith("Ez_" + fields_iteration):
                filename_Ez.append(fields + '/' + files);
        filename_Ez.sort()

        # N1
        filename_N1 = []
        for files in os.listdir(fields):
            if files.endswith(procname) and files.startswith("N1_" + fields_iteration):
                filename_N1.append(fields + '/' + files);  
        filename_N1.sort()

        # X1
        filename_X1 = []
        for files in os.listdir(fields):
            if files.endswith(procname) and files.startswith("X1_" + fields_iteration):
                filename_X1.append(fields + '/' + files);  
        filename_X1.sort()

        # Y1
        filename_Y1 = []
        for files in os.listdir(fields):
            if files.endswith(procname) and files.startswith("Y1_" + fields_iteration):
                filename_Y1.append(fields + '/' + files);  
        filename_Y1.sort()

        # Z1
        filename_Z1 = []
        for files in os.listdir(fields):
            if files.endswith(procname) and files.startswith("Z1_" + fields_iteration):
                filename_Z1.append(fields + '/' + files);  
        filename_Z1.sort()

        # N2
        filename_N2 = []
        for files in os.listdir(fields):
            if files.endswith(procname) and files.startswith("N2_" + fields_iteration):
                filename_N2.append(fields + '/' + files);
        filename_N2.sort()

        # X2
        filename_X2 = []
        for files in os.listdir(fields):
            if files.endswith(procname) and files.startswith("X2_" + fields_iteration):
                filename_X2.append(fields + '/' + files);
        filename_X2.sort()

        # Y2
        filename_Y2 = []
        for files in os.listdir(fields):
            if files.endswith(procname) and files.startswith("Y2_" + fields_iteration):
                filename_Y2.append(fields + '/' + files);
        filename_Y2.sort()

        # Z2
        filename_Z2 = []
        for files in os.listdir(fields):
            if files.endswith(procname) and files.startswith("Z2_" + fields_iteration):
                filename_Z2.append(fields + '/' + files);
        filename_Z2.sort()

        # N3
        filename_N3 = []
        for files in os.listdir(fields):
            if files.endswith(procname) and files.startswith("N3_" + fields_iteration):
                filename_N3.append(fields + '/' + files);
        filename_N3.sort()

        # X3
        filename_X3 = []
        for files in os.listdir(fields):
            if files.endswith(procname) and files.startswith("X3_" + fields_iteration):
                filename_X3.append(fields + '/' + files);
        filename_X3.sort()

        # Y3
        filename_Y3 = []
        for files in os.listdir(fields):
            if files.endswith(procname) and files.startswith("Y3_" + fields_iteration):
                filename_Y3.append(fields + '/' + files);
        filename_Y3.sort()

        # Z3
        filename_Z3 = []
        for files in os.listdir(fields):
            if files.endswith(procname) and files.startswith("Z3_" + fields_iteration):
                filename_Z3.append(fields + '/' + files);
        filename_Z3.sort()

        for kk in range(kprocs):
            for jj in range(jprocs):
                for ii in range(iprocs):
                    UE_table = loadTable(filename_UE[ii+iprocs*jj+iprocs*jprocs*kk])
                    Ex_table = loadTable(filename_Ex[ii+iprocs*jj+iprocs*jprocs*kk])
                    Ey_table = loadTable(filename_Ey[ii+iprocs*jj+iprocs*jprocs*kk])
                    Ez_table = loadTable(filename_Ez[ii+iprocs*jj+iprocs*jprocs*kk])
                    N1_table = loadTable(filename_N1[ii+iprocs*jj+iprocs*jprocs*kk])
                    X1_table = loadTable(filename_X1[ii+iprocs*jj+iprocs*jprocs*kk])
                    Y1_table = loadTable(filename_Y1[ii+iprocs*jj+iprocs*jprocs*kk])
                    Z1_table = loadTable(filename_Z1[ii+iprocs*jj+iprocs*jprocs*kk])
                    N2_table = loadTable(filename_N2[ii+iprocs*jj+iprocs*jprocs*kk])
                    X2_table = loadTable(filename_X2[ii+iprocs*jj+iprocs*jprocs*kk])
                    Y2_table = loadTable(filename_Y2[ii+iprocs*jj+iprocs*jprocs*kk])
                    Z2_table = loadTable(filename_Z2[ii+iprocs*jj+iprocs*jprocs*kk])
                    N3_table = loadTable(filename_N3[ii+iprocs*jj+iprocs*jprocs*kk])
                    X3_table = loadTable(filename_X3[ii+iprocs*jj+iprocs*jprocs*kk])
                    Y3_table = loadTable(filename_Y3[ii+iprocs*jj+iprocs*jprocs*kk])
                    Z3_table = loadTable(filename_Z3[ii+iprocs*jj+iprocs*jprocs*kk])

                    l=0
                    for k in range(Nz/kprocs+2):
                        for j in range(Ny/jprocs+2):        
                            for i in range(Nx/iprocs+2):
                                UE[i+ii*Nx/iprocs,j+jj*Ny/jprocs,k+kk*Nz/kprocs] = UE_table[l]
                                Ex[i+ii*Nx/iprocs,j+jj*Ny/jprocs,k+kk*Nz/kprocs] = Ex_table[l]
                                Ey[i+ii*Nx/iprocs,j+jj*Ny/jprocs,k+kk*Nz/kprocs] = Ey_table[l]
                                Ez[i+ii*Nx/iprocs,j+jj*Ny/jprocs,k+kk*Nz/kprocs] = Ez_table[l]
                                N1[i+ii*Nx/iprocs,j+jj*Ny/jprocs,k+kk*Nz/kprocs] = N1_table[l]
                                X1[i+ii*Nx/iprocs,j+jj*Ny/jprocs,k+kk*Nz/kprocs] = X1_table[l]
                                Y1[i+ii*Nx/iprocs,j+jj*Ny/jprocs,k+kk*Nz/kprocs] = Y1_table[l]
                                Z1[i+ii*Nx/iprocs,j+jj*Ny/jprocs,k+kk*Nz/kprocs] = Z1_table[l]
                                N2[i+ii*Nx/iprocs,j+jj*Ny/jprocs,k+kk*Nz/kprocs] = N2_table[l]
                                X2[i+ii*Nx/iprocs,j+jj*Ny/jprocs,k+kk*Nz/kprocs] = X2_table[l]
                                Y2[i+ii*Nx/iprocs,j+jj*Ny/jprocs,k+kk*Nz/kprocs] = Y2_table[l]
                                Z2[i+ii*Nx/iprocs,j+jj*Ny/jprocs,k+kk*Nz/kprocs] = Z2_table[l]
                                N3[i+ii*Nx/iprocs,j+jj*Ny/jprocs,k+kk*Nz/kprocs] = N3_table[l]
                                X3[i+ii*Nx/iprocs,j+jj*Ny/jprocs,k+kk*Nz/kprocs] = X3_table[l]
                                Y3[i+ii*Nx/iprocs,j+jj*Ny/jprocs,k+kk*Nz/kprocs] = Y3_table[l]
                                Z3[i+ii*Nx/iprocs,j+jj*Ny/jprocs,k+kk*Nz/kprocs] = Z3_table[l]
                                l=l+1
        np.savez(save_fields
            ,   UE                  = UE
            ,   Ex                  = Ex
            ,   Ey                  = Ey
            ,   Ez                  = Ez
            ,   N1                  = N1
            ,   N2                  = N2
            ,   N3                  = N3
            ,   X1                  = X1
            ,   X2                  = X2
            ,   X3                  = X3
            ,   Y1                  = Y1
            ,   Y2                  = Y2
            ,   Y3                  = Y3
            ,   Z1                  = Z1
            ,   Z2                  = Z2
            ,   Z3                  = Z3
            )

    Xval        = 0.0*x[:,0,0]
    Yval        = 0.0*x[0,:,0]
    Zval        = 0.0*x[0,0,:]

    dV          = dxyz[0]*dxyz[1]*dxyz[2]

    Jmag        = 0.0*x
    Jx          = 0.0*x
    Jy          = 0.0*x
    Jz          = 0.0*x
    Jr          = 0.0*x
    Er          = 0.0*x
    Emag        = 0.0*x
    Distance    = 0.0*x

    Jx          = X1*charge[0]+X2*charge[1]+X3*charge[2]
    Jy          = Y1*charge[0]+Y2*charge[1]+Y3*charge[2]
    Jz          = Z1*charge[0]+Z2*charge[1]+Z3*charge[2]

    Emag        = np.sqrt((Ex+E0[0])**2 + (Ey+E0[1])**2 + (Ez+E0[2])**2)
    Jmag        = np.sqrt(Jx**2 + Jy**2 + Jz**2)
    for i in range(Xval.shape[0]):
        for j in range(Yval.shape[0]):    
            for k in range(Zval.shape[0]):
                Xval[i]         = xmin+(i-0.5)*dxyz[0]
                Yval[j]         = ymin+(j-0.5)*dxyz[1]
                Zval[k]         = zmin+(k-0.5)*dxyz[2]
                Distance[i,j,k] = np.sqrt( (Xval[i])**2 + (Yval[j])**2 + (Zval[k])**2 )
                Er[i,j,k]       = (Ex[i,j,k]*Xval[i]+Ey[i,j,k]*Yval[j]+Ez[i,j,k]*Zval[k])*Distance[i,j,k]**(-1)
                Jr[i,j,k]       = (Jx[i,j,k]*Xval[i]+Jy[i,j,k]*Yval[j]+Jz[i,j,k]*Zval[k])*Distance[i,j,k]**(-1)
                if (Distance[i,j,k] > 0.0 and Qn[2] > 0.0 and vn[2] > 0.0):
                    N0[i,j,k]   = Qn[2]/( 4.0*np.pi*vn[2]*Distance[i,j,k]**2 )*np.exp( -Distance[i,j,k]*nu_d[2]/vn[2] )*dV
                

    UE_origo    = 0.125*( UE[-(Nx/2+2),Ny/2+1,Nz/2+1]+UE[-(Nx/2+2),Ny/2+1,-(Nz/2+2)]+UE[-(Nx/2+2),-(Ny/2+2),Nz/2+1]+UE[-(Nx/2+2),-(Ny/2+2),-(Nz/2+2)] 
                        + UE[  Nx/2+1 ,Ny/2+1,Nz/2+1]+UE[  Nx/2+1 ,Ny/2+1,-(Nz/2+2)]+UE[  Nx/2+1 ,-(Ny/2+2),Nz/2+1]+UE[  Nx/2+1 ,-(Ny/2+2),-(Nz/2+2)] )
    UE_X1       = 0.25*( UE[-1,Ny/2+1,Nz/2+1]+UE[-1,Ny/2+1,-(Nz/2+2)]+UE[-1,-(Ny/2+2),Nz/2+1]+UE[-1,-(Ny/2+2),-(Nz/2+2)] )
    UE_X0       = 0.25*( UE[+0,Ny/2+1,Nz/2+1]+UE[+0,Ny/2+1,-(Nz/2+2)]+UE[+0,-(Ny/2+2),Nz/2+1]+UE[+0,-(Ny/2+2),-(Nz/2+2)] )
    UE_Y1       = 0.25*( UE[Nx/2+1,-1,Nz/2+1]+UE[Nx/2+1,-1,-(Nz/2+2)]+UE[-(Nx/2+2),-1,Nz/2+1]+UE[-(Nx/2+2),-1,-(Nz/2+2)] )
    UE_Y0       = 0.25*( UE[Nx/2+1,+0,Nz/2+1]+UE[Nx/2+1,+0,-(Nz/2+2)]+UE[-(Nx/2+2),+0,Nz/2+1]+UE[-(Nx/2+2),+0,-(Nz/2+2)] )
    UE_Z1       = 0.25*( UE[Nx/2+1,Ny/2+1,-1]+UE[Nx/2+1,-(Ny/2+2),-1]+UE[-(Nx/2+2),Ny/2+1,-1]+UE[-(Nx/2+2),-(Ny/2+2),-1] )
    UE_Z0       = 0.25*( UE[Nx/2+1,Ny/2+1,+0]+UE[Nx/2+1,-(Ny/2+2),+0]+UE[-(Nx/2+2),Ny/2+1,+0]+UE[-(Nx/2+2),-(Ny/2+2),+0] )

    u0          = 1.0*v0[3]
    vth         = np.sqrt( 2.0*kb*kelvin[1]/mass[1] )

    vx0         = 0.5*u0*(1.0+math.erf(u0/vth)) + 0.5*vth/np.sqrt(np.pi)*np.exp(-(u0/vth)**2)
    vy0         = 0.5*vth/np.sqrt(np.pi)
    vxm         = u0 + vth/np.sqrt(np.pi)*np.exp(-(u0/vth)**2)*(1.0+math.erf(u0/vth))**(-1)
    vym         = vth/np.sqrt(np.pi)
    dcosrad_0   = vx0/(vx0+vy0)
    dOmega_0    = 2.0*np.pi*(1.0-dcosrad_0)

    nrden_max       = 16.0
    nrden_min       = nrden_max**(-1)
    rmin            = min(np.abs(xmin),np.abs(xmax),np.abs(ymin),np.abs(ymax),np.abs(xmin),np.abs(ymax))
    if (Qn[2] > 0.0 and vn[2] > 0.0):
        mrden_max       = ( Qn[2]/( 4.0*np.pi*vn[2]*(0.25*dxyz[0])**2 )*np.exp( -(0.25*dxyz[0])*nu_d[2]/vn[2] ) )/( density[0] )
        mrden_min       = ( Qn[2]/( 4.0*np.pi*vn[2]*rmin**2 )*np.exp( -rmin*nu_d[2]/vn[2] ) )/( density[0] )
    else:
        mrden_max       = nrden_max
        mrden_min       = nrden_min


# CHANGE OF XMIN, XMAX, etc...
    scale_length = np.abs( mass[0]*v0[3]/(charge[0]*B0[3]) )

    xmin        = xmin/scale_length
    xmax        = xmax/scale_length
    ymin        = ymin/scale_length
    ymax        = ymax/scale_length
    zmin        = zmin/scale_length
    zmax        = zmax/scale_length
    rmin        = rmin/scale_length
    dxyz[:]     = dxyz[:]/scale_length

    Rcyc        = mass[2]/mass[0]
    Rele        = 1.0

#    param       = np.linspace(0.0, min(0.5*np.pi, np.arccos( 1.0 - min(1.0, zmax/Rcyc) )), num=1024, endpoint= False)
    param       = np.linspace(0.0, 2.0*np.pi, num=1024, endpoint= True)

    Xcyc        = -Rcyc*(param-np.sin(param))
    Ycyc        = 0.0*param
    Zcyc        = Rcyc*(1.0-np.cos(param))

    param       = np.linspace(0.0, np.around(mass[2]/mass[0])*2.0*np.pi, num=np.around(mass[2]/mass[0])*1024, endpoint= True)

    Xele        = -Rele*(param-np.sin(param))
    Yele        = 0.0*param
    Zele        = -Rele*(1.0-np.cos(param))

    remove      = np.zeros(3)

    remove[0]   = 1.0+int(round( Nx*(1.5-np.sqrt(2)) ))
    remove[1]   = 1.0+int(round( Ny*(1.5-np.sqrt(2)) ))
    remove[2]   = 1.0+int(round( Nz*(1.5-np.sqrt(2)) ))

#    remove[0]   = Nx/2
#    remove[1]   = Ny/2
#    remove[2]   = Nz/2

#    remove[:]   = 1

#    remove[:]   = 0


    indeces         = [0]*6
    indeces[0]      = int(round(0.0        + remove[0]))
    indeces[1]      = int(round(x.shape[0] - remove[0]))
    indeces[2]      = int(round(0.0        + remove[1]))
    indeces[3]      = int(round(x.shape[1] - remove[1]))
    indeces[4]      = int(round(0.0        + remove[2]))
    indeces[5]      = int(round(x.shape[2] - remove[2]))


    theextent       = np.zeros((3,4))
    theextent[0]    = np.around(np.array([1.0*ymin-1.0*dxyz[1]+remove[1]*dxyz[1], 1.0*ymax+1.0*dxyz[1]-remove[1]*dxyz[1], \
                                          1.0*zmin-1.0*dxyz[2]+remove[2]*dxyz[2], 1.0*zmax+1.0*dxyz[2]-remove[2]*dxyz[2]])*1000.0)/1000.0
    theextent[1]    = np.around(np.array([1.0*xmin-1.0*dxyz[0]+remove[0]*dxyz[0], 1.0*xmax+1.0*dxyz[0]-remove[0]*dxyz[0], \
                                          1.0*zmin-1.0*dxyz[2]+remove[2]*dxyz[2], 1.0*zmax+1.0*dxyz[2]-remove[2]*dxyz[2]])*1000.0)/1000.0
    theextent[2]    = np.around(np.array([1.0*xmin-1.0*dxyz[0]+remove[0]*dxyz[0], 1.0*xmax+1.0*dxyz[0]-remove[0]*dxyz[0], \
                                          1.0*ymin-1.0*dxyz[1]+remove[1]*dxyz[1], 1.0*ymax+1.0*dxyz[1]-remove[1]*dxyz[1]])*1000.0)/1000.0
    my_xticks       = np.around(np.linspace(theextent[1,0]/1000.0,theextent[1,1]/1000.0,num=11,endpoint=True))
    my_yticks       = np.around(np.linspace(theextent[0,0]/1000.0,theextent[0,1]/1000.0,num=11,endpoint=True))
    my_zticks       = np.around(np.linspace(theextent[0,2]/1000.0,theextent[0,3]/1000.0,num=11,endpoint=True))
    my_cticks       = np.around(1000.0*np.linspace(2.0*nrden_min,2.0*nrden_max,num=21,endpoint=True))/1000.0
    my_cticks       = np.around(1000.0*np.linspace(2.0*mrden_min,2.0*mrden_max,num=21,endpoint=True))/1000.0
    ef_min          = -2.0
    ef_max          = 2.0
    my_cticks_ef    = np.around(1000.0*np.linspace(2.0*ef_min,2.0*ef_max,num=21,endpoint=True))/1000.0
    J_min           = -4.0
    J_max           = 4.0
    my_cticks_J     = np.around(1000.0*np.linspace(2.0*J_min,2.0*J_max,num=21,endpoint=True))/1000.0
    volt_min        = -1.0
    volt_max        = 1.0
    my_cticks_volt  = np.around(1000.0*np.linspace(2.0*volt_min,2.0*volt_max,num=21,endpoint=True))/1000.0

    print remove
    print Nx, Ny, Nz
    print indeces
    print theextent
    
    try:
        os.stat(save_location);
        shutil.rmtree(save_location);
    except:
        pass
    os.mkdir(save_location)


    cmap_title = r' '

    thecmap         = cm.jet

# MASS 0 (neutrals)
    N0_min = mrden_min
    N0_max = mrden_max
    N0[N0 < N0_min*dV*density[0]] = N0_min*dV*density[0]

    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    pc = imshow( 0.5*(N0[Nx/2+1,indeces[2]:indeces[3],indeces[4]:indeces[5]].T+N0[-(Nx/2+2),indeces[2]:indeces[3],indeces[4]:indeces[5]].T)/dV/density[0], aspect='equal', interpolation="nearest", origin='lower', vmin=N0_min, vmax=N0_max, extent=theextent[0], cmap=thecmap, norm=matplotlib.colors.LogNorm() )
    ax = plt.gca()
    xlabel(r'$y \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$z \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_yticks, fontsize=12)
#    plt.yticks(my_zticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{H_2O \, number \, density} \; [ \mathrm{n_0} ]$', fontsize=18)
#    title(r'$\mathrm{e^- \, number \, density} \; [ \mathrm{m^{-3}} ]$', fontsize=18)
    text(1.1, 0.65, r'$\uparrow \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\odot \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\rightarrow \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
#    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/N0_YZ_' + fields_iteration + '.'   
    saveloc_png = saveloc + 'png'
    saveloc_eps = saveloc + 'eps'
    if saveloc != 'none':
#        try:
#            savefig(saveloc_png);
#        except:
#            print saveloc_png
#            strErr = '** ERROR ** Could not save file';
#            raise NameError(strErr);
        savefig(saveloc_png);
        try:
            savefig(saveloc_eps);
        except:
            print saveloc_eps
            strErr = '** ERROR ** Could not save file';
            raise NameError(strErr);
    else:
        strErr = '** ERROR ** Save location is not specified'
        raise NameError(strErr);
    cla()
    clf()
    close('all')



    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    plot(Xcyc, Zcyc, color='k', linestyle='dotted') 
    plot(Xele, Zele, color='k', linestyle='dotted')  
    pc = imshow( 0.5*(N0[indeces[0]:indeces[1],Ny/2+1,indeces[4]:indeces[5]].T+N0[indeces[0]:indeces[1],-(Ny/2+2),indeces[4]:indeces[5]].T)/dV/density[0], aspect='equal', interpolation="nearest", origin='lower', vmin=N0_min, vmax=N0_max, extent=theextent[1], cmap=thecmap, norm=matplotlib.colors.LogNorm() )
    ax = plt.gca()
    xlabel(r'$x \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$z \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_xticks, fontsize=12)
#    plt.yticks(my_zticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{H_2O \, number \, density} \; [ \mathrm{n_0} ]$', fontsize=18)
#    title(r'$\mathrm{e^- \, number \, density} \; [ \mathrm{m^{-3}} ]$', fontsize=18)
    text(1.1, 0.65, r'$\uparrow \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\rightarrow \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\otimes \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
#    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/N0_XZ_' + fields_iteration + '.'   
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
    cla()
    clf()
    close('all')



    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    pc = imshow( 0.5*(N0[indeces[0]:indeces[1],indeces[2]:indeces[3],Nz/2+1].T+N0[indeces[0]:indeces[1],indeces[2]:indeces[3],-(Nz/2+2)].T)/dV/density[0], aspect='equal', interpolation="nearest", origin='lower', vmin=N0_min, vmax=N0_max, extent=theextent[2], cmap=thecmap, norm=matplotlib.colors.LogNorm() )
    ax = plt.gca()
    xlabel(r'$x \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$y \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_xticks, fontsize=12)
#    plt.yticks(my_yticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{H_2O \, number \, density} \; [ \mathrm{n_0} ]$', fontsize=18)
#    title(r'$\mathrm{e^- \, number \, density} \; [ \mathrm{m^{-3}} ]$', fontsize=18)
    text(1.1, 0.65, r'$\odot \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\rightarrow \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\uparrow \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
#    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/N0_XY_' + fields_iteration + '.'   
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
    cla()
    clf()
    close('all')


# MASS 1
    N1_min = nrden_min
    N1_max = nrden_max
    N1[N1 < N1_min*dV*density[0]] = N1_min*dV*density[0]

    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    pc = imshow( 0.5*(N1[Nx/2+1,indeces[2]:indeces[3],indeces[4]:indeces[5]].T+N1[-(Nx/2+2),indeces[2]:indeces[3],indeces[4]:indeces[5]].T)/dV/density[0], aspect='equal', interpolation="nearest", origin='lower', vmin=N1_min, vmax=N1_max, extent=theextent[0], cmap=thecmap, norm=matplotlib.colors.LogNorm() )
    ax = plt.gca()
    xlabel(r'$y \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$z \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_yticks, fontsize=12)
#    plt.yticks(my_zticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{e^- \, number \, density} \; [ \mathrm{n_0} ]$', fontsize=18)
#    title(r'$\mathrm{e^- \, number \, density} \; [ \mathrm{m^{-3}} ]$', fontsize=18)
    text(1.1, 0.65, r'$\uparrow \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\odot \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\rightarrow \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
#    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/N1_YZ_' + fields_iteration + '.'   
    saveloc_png = saveloc + 'png'
    saveloc_eps = saveloc + 'eps'
    if saveloc != 'none':
#        try:
#            savefig(saveloc_png);
#        except:
#            print saveloc_png
#            strErr = '** ERROR ** Could not save file';
#            raise NameError(strErr);
        savefig(saveloc_png);
        try:
            savefig(saveloc_eps);
        except:
            print saveloc_eps
            strErr = '** ERROR ** Could not save file';
            raise NameError(strErr);
    else:
        strErr = '** ERROR ** Save location is not specified'
        raise NameError(strErr);
    cla()
    clf()
    close('all')



    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    plot(Xcyc, Zcyc, color='k', linestyle='dotted') 
    plot(Xele, Zele, color='k', linestyle='dotted')  
    pc = imshow( 0.5*(N1[indeces[0]:indeces[1],Ny/2+1,indeces[4]:indeces[5]].T+N1[indeces[0]:indeces[1],-(Ny/2+2),indeces[4]:indeces[5]].T)/dV/density[0], aspect='equal', interpolation="nearest", origin='lower', vmin=N1_min, vmax=N1_max, extent=theextent[1], cmap=thecmap, norm=matplotlib.colors.LogNorm() )
    ax = plt.gca()
    xlabel(r'$x \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$z \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_xticks, fontsize=12)
#    plt.yticks(my_zticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{e^- \, number \, density} \; [ \mathrm{n_0} ]$', fontsize=18)
#    title(r'$\mathrm{e^- \, number \, density} \; [ \mathrm{m^{-3}} ]$', fontsize=18)
    text(1.1, 0.65, r'$\uparrow \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\rightarrow \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\otimes \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
#    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/N1_XZ_' + fields_iteration + '.'   
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
    cla()
    clf()
    close('all')



    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    pc = imshow( 0.5*(N1[indeces[0]:indeces[1],indeces[2]:indeces[3],Nz/2+1].T+N1[indeces[0]:indeces[1],indeces[2]:indeces[3],-(Nz/2+2)].T)/dV/density[0], aspect='equal', interpolation="nearest", origin='lower', vmin=N1_min, vmax=N1_max, extent=theextent[2], cmap=thecmap, norm=matplotlib.colors.LogNorm() )
    ax = plt.gca()
    xlabel(r'$x \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$y \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_xticks, fontsize=12)
#    plt.yticks(my_yticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{e^- \, number \, density} \; [ \mathrm{n_0} ]$', fontsize=18)
#    title(r'$\mathrm{e^- \, number \, density} \; [ \mathrm{m^{-3}} ]$', fontsize=18)
    text(1.1, 0.65, r'$\odot \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\rightarrow \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\uparrow \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
#    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/N1_XY_' + fields_iteration + '.'   
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
    cla()
    clf()
    close('all')


# MASS 2
    N2_min = nrden_min
    N2_max = nrden_max
    N2[N2 < N2_min*dV*density[0]] = N2_min*dV*density[0]


    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    pc = imshow( 0.5*(N2[Nx/2+1,indeces[2]:indeces[3],indeces[4]:indeces[5]].T+N2[-(Nx/2+2),indeces[2]:indeces[3],indeces[4]:indeces[5]].T)/dV/density[0], aspect='equal', interpolation="nearest", origin='lower', vmin=N2_min, vmax=N2_max, extent=theextent[0], cmap=thecmap, norm=matplotlib.colors.LogNorm() )
    ax = plt.gca()
    xlabel(r'$y \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$z \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_yticks, fontsize=12)
#    plt.yticks(my_zticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{H^+ \, number \, density} \; [ \mathrm{n_0} ]$', fontsize=18)
#    title(r'$\mathrm{H^+ \, number \, density} \; [ \mathrm{m^{-3}} ]$', fontsize=18)
    text(1.1, 0.65, r'$\uparrow \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\odot \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\rightarrow \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
#    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/N2_YZ_' + fields_iteration + '.'   
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
    cla()
    clf()
    close('all')

    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    plot(Xcyc, Zcyc, color='k', linestyle='dotted') 
    plot(Xele, Zele, color='k', linestyle='dotted')  
    pc = imshow( 0.5*(N2[indeces[0]:indeces[1],Ny/2+1,indeces[4]:indeces[5]].T+N2[indeces[0]:indeces[1],-(Ny/2+2),indeces[4]:indeces[5]].T)/dV/density[0], aspect='equal', interpolation="nearest", origin='lower', vmin=N2_min, vmax=N2_max, extent=theextent[1], cmap=thecmap, norm=matplotlib.colors.LogNorm() )
    ax = plt.gca()
    xlabel(r'$x \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$z \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_xticks, fontsize=12)
#    plt.yticks(my_zticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{H^+ \, number \, density} \; [ \mathrm{n_0} ]$', fontsize=18)
#    title(r'$\mathrm{H^+ \, number \, density} \; [ \mathrm{m^{-3}} ]$', fontsize=18)
    text(1.1, 0.65, r'$\uparrow \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\rightarrow \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\otimes \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
#    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/N2_XZ_' + fields_iteration + '.'   
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
    cla()
    clf()
    close('all')

    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    pc = imshow( 0.5*(N2[indeces[0]:indeces[1],indeces[2]:indeces[3],Nz/2+1].T+N2[indeces[0]:indeces[1],indeces[2]:indeces[3],-(Nz/2+2)].T)/dV/density[0], aspect='equal', interpolation="nearest", origin='lower', vmin=N2_min, vmax=N2_max, extent=theextent[2], cmap=thecmap, norm=matplotlib.colors.LogNorm() )
    ax = plt.gca()
    xlabel(r'$x \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$y \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_xticks, fontsize=12)
#    plt.yticks(my_yticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{H^+ \, number \, density} \; [ \mathrm{n_0} ]$', fontsize=18)
#    title(r'$\mathrm{H^+ \, number \, density} \; [ \mathrm{m^{-3}} ]$', fontsize=18)
    text(1.1, 0.65, r'$\odot \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\rightarrow \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\uparrow \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
#    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/N2_XY_' + fields_iteration + '.'   
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
    cla()
    clf()
    close('all')


# MASS 3

    N3_min = nrden_min
    N3_max = nrden_max
    N3[N3 < N3_min] = N3_min


    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    pc = imshow( 0.5*(N3[Nx/2+1,indeces[2]:indeces[3],indeces[4]:indeces[5]].T+N3[-(Nx/2+2),indeces[2]:indeces[3],indeces[4]:indeces[5]].T)/dV/density[0], aspect='equal', interpolation="nearest", origin='lower', vmin=N3_min, vmax=N3_max, extent=theextent[0], cmap=thecmap, norm=matplotlib.colors.LogNorm() )
    ax = plt.gca()
    xlabel(r'$y \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$z \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_yticks, fontsize=12)
#    plt.yticks(my_zticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{H_2O^+ \, number \, density} \; [ \mathrm{n_0} ]$', fontsize=18)
#    title(r'$\mathrm{H_2O^+ \, number \, density} \; [ \mathrm{m^{-3}} ]$', fontsize=18)
    text(1.1, 0.65, r'$\uparrow \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\odot \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\rightarrow \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm())
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
#    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/N3_YZ_' + fields_iteration + '.'  
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
    cla()
    clf()
    close('all')



    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    plot(Xcyc, Zcyc, color='w', linestyle='dotted') 
    plot(Xele, Zele, color='w', linestyle='dotted')  
    pc = imshow( 0.5*(N3[indeces[0]:indeces[1],Ny/2+1,indeces[4]:indeces[5]].T+N3[indeces[0]:indeces[1],-(Ny/2+2),indeces[4]:indeces[5]].T)/dV/density[0], aspect='equal', interpolation="nearest", origin='lower', vmin=N3_min, vmax=N3_max, extent=theextent[1], cmap=thecmap, norm=matplotlib.colors.LogNorm() )
    ax = plt.gca()
    xlabel(r'$x \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$z \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_xticks, fontsize=12)
#    plt.yticks(my_zticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{H_2O^+ \, number \, density} \; [ \mathrm{n_0} ]$', fontsize=18)
#    title(r'$\mathrm{H_2O^+ \, number \, density} \; [ \mathrm{m^{-3}} ]$', fontsize=18)
    text(1.1, 0.65, r'$\uparrow \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\rightarrow \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\otimes \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm() )
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
#    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/N3_XZ_' + fields_iteration + '.'  
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
    cla()
    clf()
    close('all')



    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    pc = imshow( 0.5*(N3[indeces[0]:indeces[1],indeces[2]:indeces[3],Nz/2+1].T+N3[indeces[0]:indeces[1],indeces[2]:indeces[3],-(Nz/2+2)].T)/dV/density[0], aspect='equal', interpolation="nearest", origin='lower', vmin=N3_min, vmax=N3_max, extent=theextent[2], cmap=thecmap, norm=matplotlib.colors.LogNorm() )
    ax = plt.gca()

    xlabel(r'$x \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$y \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_xticks, fontsize=12)
#    plt.yticks(my_yticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{H_2O^+ \, number \, density} \; [ \mathrm{n_0} ]$', fontsize=18)
#    title(r'$\mathrm{H_2O^+ \, number \, density} \; [ \mathrm{m^{-3}} ]$', fontsize=18)
    text(1.1, 0.65, r'$\odot \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\rightarrow \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\uparrow \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm())
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
#    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/N3_XY_' + fields_iteration + '.'  
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
    cla()
    clf()
    close('all')


    thecmap         = cm.bwr

# ELECTRIC POTENTIAL, U
    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    pc = imshow( 0.5*(UE[Nx/2+1,indeces[2]:indeces[3],indeces[4]:indeces[5]].T+UE[-(Nx/2+2),indeces[2]:indeces[3],indeces[4]:indeces[5]].T)/UE0, aspect='equal', interpolation="nearest", origin='lower', vmin=volt_min, vmax=volt_max, extent=theextent[0], cmap=thecmap )
    ax = plt.gca()
    xlabel(r'$y \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$z \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_yticks, fontsize=12)
#    plt.yticks(my_zticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{U} \; [ \mathrm{U}_0 ]$', fontsize=18)
    text(1.1, 0.65, r'$\uparrow \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\odot \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\rightarrow \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks_volt)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/UE_YZ_' + fields_iteration + '.'   
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
    cla()
    clf()
    close('all')

    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    plot(Xcyc, Zcyc, color='k', linestyle='dotted') 
    plot(Xele, Zele, color='k', linestyle='dotted')  
    plot(Xele, Zele, color='k', linestyle='dotted') 
    pc = imshow( 0.5*(UE[indeces[0]:indeces[1],Ny/2+1,indeces[4]:indeces[5]].T+UE[indeces[0]:indeces[1],-(Ny/2+2),indeces[4]:indeces[5]].T)/UE0, aspect='equal', interpolation="nearest", origin='lower', vmin=volt_min, vmax=volt_max, extent=theextent[1], cmap=thecmap )
    ax = plt.gca()
    xlabel(r'$x \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$z \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_xticks, fontsize=12)
#    plt.yticks(my_zticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{U} \; [ \mathrm{U}_0 ]$', fontsize=18)
    text(1.1, 0.65, r'$\uparrow \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\rightarrow \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\otimes \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks_volt)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/UE_XZ_' + fields_iteration + '.'   
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
    cla()
    clf()
    close('all')

    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
#    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
#    vlines(-3.0-1.5*dxyz[0], -3.0, 3.0, color='k', linestyle='dashed')
#    hlines(-3.0, 3.0-1.5*dxyz[0], -3.0-1.5*dxyz[0], color='k', linestyle='dashed')
#    vlines(3.0-1.5*dxyz[0], 3.0, -3.0, color='k', linestyle='dashed')
#    hlines(3.0, -3.0-1.5*dxyz[0], 3.0-1.5*dxyz[0], color='k', linestyle='dashed')
    pc = imshow( 0.5*(UE[indeces[0]:indeces[1],indeces[2]:indeces[3],Nz/2+1].T+UE[indeces[0]:indeces[1],indeces[2]:indeces[3],-(Nz/2+2)].T)/UE0, aspect='equal', interpolation="nearest", origin='lower', vmin=volt_min, vmax=volt_max, extent=theextent[2], cmap=thecmap )
    ax = plt.gca()
    xlabel(r'$x \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$y \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_xticks, fontsize=12)
#    plt.yticks(my_yticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{U} \; [ \mathrm{U}_0 ]$', fontsize=18)
    text(1.1, 0.65, r'$\odot \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\rightarrow \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\uparrow \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks_volt)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/UE_XY_' + fields_iteration + '.'   
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
    cla()
    clf()
    close('all')


# ELECTRIC POTENTIAL, U, scaled
    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    pc = imshow( (0.5*(UE[Nx/2+1,indeces[2]:indeces[3],indeces[4]:indeces[5]].T+UE[-(Nx/2+2),indeces[2]:indeces[3],indeces[4]:indeces[5]].T)-UE_X1)/UE0, aspect='equal', interpolation="nearest", origin='lower', vmin=volt_min, vmax=volt_max, extent=theextent[0], cmap=thecmap )
    ax = plt.gca()
    xlabel(r'$y \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$z \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_yticks, fontsize=12)
#    plt.yticks(my_zticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{U} \; [ \mathrm{U}_0 ]$', fontsize=18)
    text(1.1, 0.65, r'$\uparrow \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\odot \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\rightarrow \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks_volt)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/UD_YZ_' + fields_iteration + '.'   
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
    cla()
    clf()
    close('all')

    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    plot(Xcyc, Zcyc, color='k', linestyle='dotted') 
    plot(Xele, Zele, color='k', linestyle='dotted')  
    plot(Xele, Zele, color='k', linestyle='dotted') 
    pc = imshow( (0.5*(UE[indeces[0]:indeces[1],Ny/2+1,indeces[4]:indeces[5]].T+UE[indeces[0]:indeces[1],-(Ny/2+2),indeces[4]:indeces[5]].T)-UE_X1)/UE0, aspect='equal', interpolation="nearest", origin='lower', vmin=volt_min, vmax=volt_max, extent=theextent[1], cmap=thecmap )
    ax = plt.gca()
    xlabel(r'$x \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$z \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_xticks, fontsize=12)
#    plt.yticks(my_zticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{U} \; [ \mathrm{U}_0 ]$', fontsize=18)
    text(1.1, 0.65, r'$\uparrow \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\rightarrow \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\otimes \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks_volt)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/UD_XZ_' + fields_iteration + '.'   
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
    cla()
    clf()
    close('all')

    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
#    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
#    vlines(-3.0-1.5*dxyz[0], -3.0, 3.0, color='k', linestyle='dashed')
#    hlines(-3.0, 3.0-1.5*dxyz[0], -3.0-1.5*dxyz[0], color='k', linestyle='dashed')
#    vlines(3.0-1.5*dxyz[0], 3.0, -3.0, color='k', linestyle='dashed')
#    hlines(3.0, -3.0-1.5*dxyz[0], 3.0-1.5*dxyz[0], color='k', linestyle='dashed')
    pc = imshow( (0.5*(UE[indeces[0]:indeces[1],indeces[2]:indeces[3],Nz/2+1].T+UE[indeces[0]:indeces[1],indeces[2]:indeces[3],-(Nz/2+2)].T)-UE_X1)/UE0, aspect='equal', interpolation="nearest", origin='lower', vmin=volt_min, vmax=volt_max, extent=theextent[2], cmap=thecmap )
    ax = plt.gca()
    xlabel(r'$x \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$y \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_xticks, fontsize=12)
#    plt.yticks(my_yticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{U} \; [ \mathrm{U}_0 ]$', fontsize=18)
    text(1.1, 0.65, r'$\odot \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\rightarrow \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\uparrow \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks_volt)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/UD_XY_' + fields_iteration + '.'   
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
    cla()
    clf()
    close('all')


# ELECTRIC POTENTIAL, U, center
    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    pc = imshow( (0.5*(UE[Nx/2+1,indeces[2]:indeces[3],indeces[4]:indeces[5]].T+UE[-(Nx/2+2),indeces[2]:indeces[3],indeces[4]:indeces[5]].T)-UE_origo)/UE0, aspect='equal', interpolation="nearest", origin='lower', vmin=volt_min, vmax=volt_max, extent=theextent[0], cmap=thecmap )
    ax = plt.gca()
    xlabel(r'$y \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$z \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_yticks, fontsize=12)
#    plt.yticks(my_zticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{U} \; [ \mathrm{U}_0 ]$', fontsize=18)
    text(1.1, 0.65, r'$\uparrow \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\odot \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\rightarrow \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks_volt)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/UC_YZ_' + fields_iteration + '.'   
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
    cla()
    clf()
    close('all')

    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    plot(Xcyc, Zcyc, color='k', linestyle='dotted') 
    plot(Xele, Zele, color='k', linestyle='dotted')  
    plot(Xele, Zele, color='k', linestyle='dotted') 
    pc = imshow( (0.5*(UE[indeces[0]:indeces[1],Ny/2+1,indeces[4]:indeces[5]].T+UE[indeces[0]:indeces[1],-(Ny/2+2),indeces[4]:indeces[5]].T)-UE_origo)/UE0, aspect='equal', interpolation="nearest", origin='lower', vmin=volt_min, vmax=volt_max, extent=theextent[1], cmap=thecmap )
    ax = plt.gca()
    xlabel(r'$x \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$z \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_xticks, fontsize=12)
#    plt.yticks(my_zticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{U} \; [ \mathrm{U}_0 ]$', fontsize=18)
    text(1.1, 0.65, r'$\uparrow \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\rightarrow \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\otimes \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks_volt)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/UC_XZ_' + fields_iteration + '.'   
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
    cla()
    clf()
    close('all')

    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
#    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
#    vlines(-3.0-1.5*dxyz[0], -3.0, 3.0, color='k', linestyle='dashed')
#    hlines(-3.0, 3.0-1.5*dxyz[0], -3.0-1.5*dxyz[0], color='k', linestyle='dashed')
#    vlines(3.0-1.5*dxyz[0], 3.0, -3.0, color='k', linestyle='dashed')
#    hlines(3.0, -3.0-1.5*dxyz[0], 3.0-1.5*dxyz[0], color='k', linestyle='dashed')
    pc = imshow( (0.5*(UE[indeces[0]:indeces[1],indeces[2]:indeces[3],Nz/2+1].T+UE[indeces[0]:indeces[1],indeces[2]:indeces[3],-(Nz/2+2)].T)-UE_origo)/UE0, aspect='equal', interpolation="nearest", origin='lower', vmin=volt_min, vmax=volt_max, extent=theextent[2], cmap=thecmap )
    ax = plt.gca()
    xlabel(r'$x \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$y \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_xticks, fontsize=12)
#    plt.yticks(my_yticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{U} \; [ \mathrm{U}_0 ]$', fontsize=18)
    text(1.1, 0.65, r'$\odot \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\rightarrow \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\uparrow \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks_volt)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/UC_XY_' + fields_iteration + '.'   
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
    cla()
    clf()
    close('all')


    thecmap         = cm.bwr

# ELECTRIC FIELD X
    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    pc = imshow( 0.5*(Ex[Nx/2+1,indeces[2]:indeces[3],indeces[4]:indeces[5]].T+Ex[-(Nx/2+2),indeces[2]:indeces[3],indeces[4]:indeces[5]].T)/E0[3], aspect='equal', interpolation="nearest", origin='lower', vmin=ef_min, vmax=ef_max, extent=theextent[0], cmap=thecmap )
    ax = plt.gca()
    xlabel(r'$y \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$z \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_yticks, fontsize=12)
#    plt.yticks(my_zticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{E_x} \; [ \mathrm{ E_{0} } ]$', fontsize=18)
#    title(r'$\mathrm{E_x} \; [ \mathrm{V \, m^{-1}} ]$', fontsize=18)
    text(1.1, 0.65, r'$\uparrow \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\odot \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\rightarrow \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks_ef)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/Ex_YZ_' + fields_iteration + '.'   
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
    cla()
    clf()
    close('all')


    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    plot(Xcyc, Zcyc, color='k', linestyle='dotted') 
    plot(Xele, Zele, color='k', linestyle='dotted')  
    pc = imshow( 0.5*(Ex[indeces[0]:indeces[1],Ny/2+1,indeces[4]:indeces[5]].T+Ex[indeces[0]:indeces[1],-(Ny/2+2),indeces[4]:indeces[5]].T)/E0[3], aspect='equal', interpolation="nearest", origin='lower', vmin=ef_min, vmax=ef_max, extent=theextent[1], cmap=thecmap )
    ax = plt.gca()
    xlabel(r'$x \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$z \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_xticks, fontsize=12)
#    plt.yticks(my_zticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{E_x} \; [ \mathrm{ E_{0} } ]$', fontsize=18)

#    title(r'$\mathrm{E_x} \; [ \mathrm{V \, m^{-1}} ]$', fontsize=18)
    text(1.1, 0.65, r'$\uparrow \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\rightarrow \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\otimes \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks_ef)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/Ex_XZ_' + fields_iteration + '.'   
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
    cla()
    clf()
    close('all')


    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
#    vlines(-3.0-1.5*dxyz[0], -3.0, 3.0, color='k', linestyle='dashed')
#    hlines(-3.0, 3.0-1.5*dxyz[0], -3.0-1.5*dxyz[0], color='k', linestyle='dashed')
#    vlines(3.0-1.5*dxyz[0], 3.0, -3.0, color='k', linestyle='dashed')
#    hlines(3.0, -3.0-1.5*dxyz[0], 3.0-1.5*dxyz[0], color='k', linestyle='dashed')
    pc = imshow( 0.5*(Ex[indeces[0]:indeces[1],indeces[2]:indeces[3],Nz/2+1].T+Ex[indeces[0]:indeces[1],indeces[2]:indeces[3],-(Nz/2+2)].T)/E0[3], aspect='equal', interpolation="nearest", origin='lower', vmin=ef_min, vmax=ef_max, extent=theextent[2], cmap=thecmap )
    ax = plt.gca()
    xlabel(r'$x \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$y \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_xticks, fontsize=12)
#    plt.yticks(my_yticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{E_x} \; [ \mathrm{ E_{0} } ]$', fontsize=18)
#    title(r'$\mathrm{E_x} \; [ \mathrm{V \, m^{-1}} ]$', fontsize=18)
    text(1.1, 0.65, r'$\odot \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\rightarrow \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\uparrow \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks_ef)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/Ex_XY_' + fields_iteration + '.'   
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
    cla()
    clf()
    close('all')


# ELECTRIC FIELD Y
    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    pc = imshow( 0.5*(Ey[Nx/2+1,indeces[2]:indeces[3],indeces[4]:indeces[5]].T+Ey[-(Nx/2+2),indeces[2]:indeces[3],indeces[4]:indeces[5]].T)/E0[3], aspect='equal', interpolation="nearest", origin='lower', vmin=ef_min, vmax=ef_max, extent=theextent[0], cmap=thecmap )
    ax = plt.gca()
    xlabel(r'$y \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$z \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_yticks, fontsize=12)
#    plt.yticks(my_zticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{E_y} \; [ \mathrm{ E_{0} } ]$', fontsize=18)
#    title(r'$\mathrm{E_y} \; [ \mathrm{V \, m^{-1}} ]$', fontsize=18)
    text(1.1, 0.65, r'$\uparrow \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\odot \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\rightarrow \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks_ef)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')

    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/Ey_YZ_' + fields_iteration + '.'   
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
    cla()
    clf()
    close('all')

    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    plot(Xcyc, Zcyc, color='k', linestyle='dotted') 
    plot(Xele, Zele, color='k', linestyle='dotted')  
    pc = imshow( 0.5*(Ey[indeces[0]:indeces[1],Ny/2+1,indeces[4]:indeces[5]].T+Ey[indeces[0]:indeces[1],-(Ny/2+2),indeces[4]:indeces[5]].T)/E0[3], aspect='equal', interpolation="nearest", origin='lower', vmin=ef_min, vmax=ef_max, extent=theextent[1], cmap=thecmap )

    ax = plt.gca()
    xlabel(r'$x \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$z \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_xticks, fontsize=12)
#    plt.yticks(my_zticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{E_y} \; [ \mathrm{ E_{0} } ]$', fontsize=18)
#    title(r'$\mathrm{E_y} \; [ \mathrm{V \, m^{-1}} ]$', fontsize=18)
    text(1.1, 0.65, r'$\uparrow \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\rightarrow \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\otimes \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks_ef)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/Ey_XZ_' + fields_iteration + '.'   
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
    cla()
    clf()
    close('all')

    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
#    vlines(-3.0-1.5*dxyz[0], -3.0, 3.0, color='k', linestyle='dashed')
#    hlines(-3.0, 3.0-1.5*dxyz[0], -3.0-1.5*dxyz[0], color='k', linestyle='dashed')
#    vlines(3.0-1.5*dxyz[0], 3.0, -3.0, color='k', linestyle='dashed')
#    hlines(3.0, -3.0-1.5*dxyz[0], 3.0-1.5*dxyz[0], color='k', linestyle='dashed')
    pc = imshow( 0.5*(Ey[indeces[0]:indeces[1],indeces[2]:indeces[3],Nz/2+1].T+Ey[indeces[0]:indeces[1],indeces[2]:indeces[3],-(Nz/2+2)].T)/E0[3], aspect='equal', interpolation="nearest", origin='lower', vmin=ef_min, vmax=ef_max, extent=theextent[2], cmap=thecmap )
    ax = plt.gca()
    xlabel(r'$x \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$y \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_xticks, fontsize=12)
#    plt.yticks(my_yticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{E_y} \; [ \mathrm{ E_{0} } ]$', fontsize=18)
#    title(r'$\mathrm{E_y} \; [ \mathrm{V \, m^{-1}} ]$', fontsize=18)
    text(1.1, 0.65, r'$\odot \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\rightarrow \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\uparrow \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks_ef)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/Ey_XY_' + fields_iteration + '.'   
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
    cla()
    clf()
    close('all')


# ELECTRIC FIELD Z
    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    pc = imshow( 0.5*(Ez[Nx/2+1,indeces[2]:indeces[3],indeces[4]:indeces[5]].T+Ez[-(Nx/2+2),indeces[2]:indeces[3],indeces[4]:indeces[5]].T)/E0[3], aspect='equal', interpolation="nearest", origin='lower', vmin=ef_min, vmax=ef_max, extent=theextent[0], cmap=thecmap )
    ax = plt.gca()
    xlabel(r'$y \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$z \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_yticks, fontsize=12)
#    plt.yticks(my_zticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{E_z} \; [ \mathrm{ E_{0} } ]$', fontsize=18)
#   title(r'$\mathrm{E_z} \; [ \mathrm{V \, m^{-1}} ]$', fontsize=18)
    text(1.1, 0.65, r'$\uparrow \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\odot \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\rightarrow \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks_ef)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/Ez_YZ_' + fields_iteration + '.'   
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
    cla()
    clf()
    close('all')

    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    plot(Xcyc, Zcyc, color='k', linestyle='dotted') 
    plot(Xele, Zele, color='k', linestyle='dotted')  
    pc = imshow( 0.5*(Ez[indeces[0]:indeces[1],Ny/2+1,indeces[4]:indeces[5]].T+Ez[indeces[0]:indeces[1],-(Ny/2+2),indeces[4]:indeces[5]].T)/E0[3], aspect='equal', interpolation="nearest", origin='lower', vmin=ef_min, vmax=ef_max, extent=theextent[1], cmap=thecmap )
    ax = plt.gca()
    xlabel(r'$x \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$z \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_xticks, fontsize=12)
#    plt.yticks(my_zticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{E_z} \; [ \mathrm{ E_{0} } ]$', fontsize=18)
#   title(r'$\mathrm{E_z} \; [ \mathrm{V \, m^{-1}} ]$', fontsize=18)
    text(1.1, 0.65, r'$\uparrow \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\rightarrow \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\otimes \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks_ef)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/Ez_XZ_' + fields_iteration + '.'   
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
    cla()
    clf()
    close('all')

    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    pc = imshow( 0.5*(Ez[indeces[0]:indeces[1],indeces[2]:indeces[3],Nz/2+1].T+Ez[indeces[0]:indeces[1],indeces[2]:indeces[3],-(Nz/2+2)].T)/E0[3], aspect='equal', interpolation="nearest", origin='lower', vmin=ef_min, vmax=ef_max, extent=theextent[2], cmap=thecmap )
    ax = plt.gca()
    xlabel(r'$x \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$y \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_xticks, fontsize=12)
#    plt.yticks(my_yticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{E_z} \; [ \mathrm{ E_{0} } ]$', fontsize=18)
#   title(r'$\mathrm{E_z} \; [ \mathrm{V \, m^{-1}} ]$', fontsize=18)
    text(1.1, 0.65, r'$\odot \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\rightarrow \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\uparrow \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks_ef)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/Ez_XY_' + fields_iteration + '.'   
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
    cla()
    clf()
    close('all')

# ELECTRIC FIELD RADIAL
    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    pc = imshow( 0.5*(Er[Nx/2+1,indeces[2]:indeces[3],indeces[4]:indeces[5]].T+Er[-(Nx/2+2),indeces[2]:indeces[3],indeces[4]:indeces[5]].T)/E0[3], aspect='equal', interpolation="nearest", origin='lower', vmin=ef_min, vmax=ef_max, extent=theextent[0], cmap=thecmap )

    ax = plt.gca()
    xlabel(r'$y \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$z \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_yticks, fontsize=12)
#    plt.yticks(my_zticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{E_r} \; [ \mathrm{ E_{0} } ]$', fontsize=18)
#   title(r'$\mathrm{E_r} \; [ \mathrm{V \, m^{-1}} ]$', fontsize=18)
    text(1.1, 0.65, r'$\uparrow \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\odot \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\rightarrow \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks_ef)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/Er_YZ_' + fields_iteration + '.'   
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
    cla()
    clf()
    close('all')

    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    plot(Xcyc, Zcyc, color='k', linestyle='dotted') 
    plot(Xele, Zele, color='k', linestyle='dotted')  
    pc = imshow( 0.5*(Er[indeces[0]:indeces[1],Ny/2+1,indeces[4]:indeces[5]].T+Er[indeces[0]:indeces[1],-(Ny/2+2),indeces[4]:indeces[5]].T)/E0[3], aspect='equal', interpolation="nearest", origin='lower', vmin=ef_min, vmax=ef_max, extent=theextent[1], cmap=thecmap )
    ax = plt.gca()
    xlabel(r'$x \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$z \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_xticks, fontsize=12)
#    plt.yticks(my_zticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{E_r} \; [ \mathrm{ E_{0} } ]$', fontsize=18)
#   title(r'$\mathrm{E_r} \; [ \mathrm{V \, m^{-1}} ]$', fontsize=18)
    text(1.1, 0.65, r'$\uparrow \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\rightarrow \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\otimes \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks_ef)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/Er_XZ_' + fields_iteration + '.'   
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
    cla()
    clf()
    close('all')

    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    pc = imshow( 0.5*(Er[indeces[0]:indeces[1],indeces[2]:indeces[3],Nz/2+1].T+Er[indeces[0]:indeces[1],indeces[2]:indeces[3],-(Nz/2+2)].T)/E0[3], aspect='equal', interpolation="nearest", origin='lower', vmin=ef_min, vmax=ef_max, extent=theextent[2], cmap=thecmap )
    ax = plt.gca()
    xlabel(r'$x \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$y \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_xticks, fontsize=12)
#    plt.yticks(my_yticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{E_r} \; [ \mathrm{ E_{0} } ]$', fontsize=18)
#   title(r'$\mathrm{E_r} \; [ \mathrm{V \, m^{-1}} ]$', fontsize=18)
    text(1.1, 0.65, r'$\odot \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\rightarrow \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\uparrow \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])

#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks_ef)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/Er_XY_' + fields_iteration + '.'   
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
    cla()
    clf()
    close('all')

    thecmap         = cm.jet

# ELECTRIC FIELD MAG
    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    pc = imshow( 0.5*(Emag[Nx/2+1,indeces[2]:indeces[3],indeces[4]:indeces[5]].T+Emag[-(Nx/2+2),indeces[2]:indeces[3],indeces[4]:indeces[5]].T)/E0[3], aspect='equal', interpolation="nearest", origin='lower', vmin=1.0+ef_min, vmax=1.0+ef_max, extent=theextent[0], cmap=thecmap )
    ax = plt.gca()
    xlabel(r'$y \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$z \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_yticks, fontsize=12)
#    plt.yticks(my_zticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{E} \; [ \mathrm{ E_{0} } ]$', fontsize=18)
#   title(r'$\mathrm{E} \; [ \mathrm{V \, m^{-1}} ]$', fontsize=18)
    text(1.1, 0.65, r'$\uparrow \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\odot \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\rightarrow \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks_ef)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/Emag_YZ_' + fields_iteration + '.'   
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
    cla()
    clf()
    close('all')

    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    plot(Xcyc, Zcyc, color='k', linestyle='dotted') 
    plot(Xele, Zele, color='k', linestyle='dotted')  
    pc = imshow( 0.5*(Emag[indeces[0]:indeces[1],Ny/2+1,indeces[4]:indeces[5]].T+Emag[indeces[0]:indeces[1],-(Ny/2+2),indeces[4]:indeces[5]].T)/E0[3], aspect='equal', interpolation="nearest", origin='lower', vmin=1.0+ef_min, vmax=1.0+ef_max, extent=theextent[1], cmap=thecmap )
    ax = plt.gca()
    xlabel(r'$x \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$z \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_xticks, fontsize=12)
#    plt.yticks(my_zticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{E} \; [ \mathrm{ E_{0} } ]$', fontsize=18)
#   title(r'$\mathrm{E} \; [ \mathrm{V \, m^{-1}} ]$', fontsize=18)
    text(1.1, 0.65, r'$\uparrow \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\rightarrow \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\otimes \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks_ef)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/Emag_XZ_' + fields_iteration + '.'   
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
    cla()
    clf()
    close('all')

    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    pc = imshow( 0.5*(Emag[indeces[0]:indeces[1],indeces[2]:indeces[3],Nz/2+1].T+Emag[indeces[0]:indeces[1],indeces[2]:indeces[3],-(Nz/2+2)].T)/E0[3], aspect='equal', interpolation="nearest", origin='lower', vmin=1.0+ef_min, vmax=1.0+ef_max, extent=theextent[2], cmap=thecmap )
    ax = plt.gca()
    xlabel(r'$x \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$y \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_xticks, fontsize=12)
#    plt.yticks(my_yticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{E} \; [ \mathrm{ E_{0} } ]$', fontsize=18)
#   title(r'$\mathrm{E} \; [ \mathrm{V \, m^{-1}} ]$', fontsize=18)
    text(1.1, 0.65, r'$\odot \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\rightarrow \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\uparrow \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks_ef)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/Emag_XY_' + fields_iteration + '.'   
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
    cla()
    clf()
    close('all')

    thecmap         = cm.bwr

# CURRENT X
    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    pc = imshow( 0.5*(Jx[Nx/2+1,indeces[2]:indeces[3],indeces[4]:indeces[5]].T+Jx[-(Nx/2+2),indeces[2]:indeces[3],indeces[4]:indeces[5]].T)/J0[3], aspect='equal', interpolation="nearest", origin='lower', vmin=J_min, vmax=J_max, extent=theextent[0], cmap=thecmap )
    ax = plt.gca()
    xlabel(r'$y \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$z \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_yticks, fontsize=12)
#    plt.yticks(my_zticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{J_x} \; [ \mathrm{ J_{0} } ]$', fontsize=18)
#    title(r'$\mathrm{E_x} \; [ \mathrm{V \, m^{-1}} ]$', fontsize=18)
    text(1.1, 0.65, r'$\uparrow \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\odot \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\rightarrow \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks_J)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/Jx_YZ_' + fields_iteration + '.'   
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
    cla()
    clf()
    close('all')


    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    plot(Xcyc, Zcyc, color='k', linestyle='dotted') 
    plot(Xele, Zele, color='k', linestyle='dotted')  
    pc = imshow( 0.5*(Jx[indeces[0]:indeces[1],Ny/2+1,indeces[4]:indeces[5]].T+Jx[indeces[0]:indeces[1],-(Ny/2+2),indeces[4]:indeces[5]].T)/J0[3], aspect='equal', interpolation="nearest", origin='lower', vmin=J_min, vmax=J_max, extent=theextent[1], cmap=thecmap )
    ax = plt.gca()
    xlabel(r'$x \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$z \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_xticks, fontsize=12)
#    plt.yticks(my_zticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{J_x} \; [ \mathrm{ J_{0} } ]$', fontsize=18)
#    title(r'$\mathrm{E_x} \; [ \mathrm{V \, m^{-1}} ]$', fontsize=18)
    text(1.1, 0.65, r'$\uparrow \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\rightarrow \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\otimes \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks_J)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/Jx_XZ_' + fields_iteration + '.'   
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
    cla()
    clf()
    close('all')


    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    pc = imshow( 0.5*(Jx[indeces[0]:indeces[1],indeces[2]:indeces[3],Nz/2+1].T+Jx[indeces[0]:indeces[1],indeces[2]:indeces[3],-(Nz/2+2)].T)/J0[3], aspect='equal', interpolation="nearest", origin='lower', vmin=J_min, vmax=J_max, extent=theextent[2], cmap=thecmap )
    ax = plt.gca()
    xlabel(r'$x \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$y \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_xticks, fontsize=12)
#    plt.yticks(my_yticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{J_x} \; [ \mathrm{ J_{0} } ]$', fontsize=18)
#    title(r'$\mathrm{E_x} \; [ \mathrm{V \, m^{-1}} ]$', fontsize=18)
    text(1.1, 0.65, r'$\odot \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\rightarrow \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\uparrow \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks_J)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/Jx_XY_' + fields_iteration + '.'   
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
    cla()
    clf()
    close('all')


# CURRENT Y
    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    pc = imshow( 0.5*(Jy[Nx/2+1,indeces[2]:indeces[3],indeces[4]:indeces[5]].T+Jy[-(Nx/2+2),indeces[2]:indeces[3],indeces[4]:indeces[5]].T)/J0[3], aspect='equal', interpolation="nearest", origin='lower', vmin=J_min, vmax=J_max, extent=theextent[0], cmap=thecmap )
    ax = plt.gca()
    xlabel(r'$y \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$z \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_yticks, fontsize=12)
#    plt.yticks(my_zticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{J_y} \; [ \mathrm{ J_{0} } ]$', fontsize=18)
#    title(r'$\mathrm{E_y} \; [ \mathrm{V \, m^{-1}} ]$', fontsize=18)
    text(1.1, 0.65, r'$\uparrow \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\odot \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\rightarrow \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks_J)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/Jy_YZ_' + fields_iteration + '.'   
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
    cla()
    clf()
    close('all')

    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    plot(Xcyc, Zcyc, color='k', linestyle='dotted') 
    plot(Xele, Zele, color='k', linestyle='dotted')  
    pc = imshow( 0.5*(Jy[indeces[0]:indeces[1],Ny/2+1,indeces[4]:indeces[5]].T+Jy[indeces[0]:indeces[1],-(Ny/2+2),indeces[4]:indeces[5]].T)/J0[3], aspect='equal', interpolation="nearest", origin='lower', vmin=J_min, vmax=J_max, extent=theextent[1], cmap=thecmap )
    ax = plt.gca()
    xlabel(r'$x \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$z \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_xticks, fontsize=12)
#    plt.yticks(my_zticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{J_y} \; [ \mathrm{ J_{0} } ]$', fontsize=18)
#    title(r'$\mathrm{E_y} \; [ \mathrm{V \, m^{-1}} ]$', fontsize=18)
    text(1.1, 0.65, r'$\uparrow \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\rightarrow \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\otimes \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks_J)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/Jy_XZ_' + fields_iteration + '.'   
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
    cla()
    clf()
    close('all')

    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    pc = imshow( 0.5*(Jy[indeces[0]:indeces[1],indeces[2]:indeces[3],Nz/2+1].T+Jy[indeces[0]:indeces[1],indeces[2]:indeces[3],-(Nz/2+2)].T)/J0[3], aspect='equal', interpolation="nearest", origin='lower', vmin=J_min, vmax=J_max, extent=theextent[2], cmap=thecmap )
    ax = plt.gca()
    xlabel(r'$x \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$y \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_xticks, fontsize=12)
#    plt.yticks(my_yticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{J_y} \; [ \mathrm{ J_{0} } ]$', fontsize=18)
#    title(r'$\mathrm{E_y} \; [ \mathrm{V \, m^{-1}} ]$', fontsize=18)
    text(1.1, 0.65, r'$\odot \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\rightarrow \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\uparrow \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks_J)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/Jy_XY_' + fields_iteration + '.'   
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
    cla()
    clf()
    close('all')


# CURRENT Z
    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    pc = imshow( 0.5*(Jz[Nx/2+1,indeces[2]:indeces[3],indeces[4]:indeces[5]].T+Jz[-(Nx/2+2),indeces[2]:indeces[3],indeces[4]:indeces[5]].T)/J0[3], aspect='equal', interpolation="nearest", origin='lower', vmin=J_min, vmax=J_max, extent=theextent[0], cmap=thecmap )
    ax = plt.gca()
    xlabel(r'$y \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$z \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_yticks, fontsize=12)
#    plt.yticks(my_zticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{J_z} \; [ \mathrm{ J_{0} } ]$', fontsize=18)
#   title(r'$\mathrm{E_z} \; [ \mathrm{V \, m^{-1}} ]$', fontsize=18)
    text(1.1, 0.65, r'$\uparrow \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\odot \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\rightarrow \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks_J)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/Jz_YZ_' + fields_iteration + '.'   
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
    cla()
    clf()
    close('all')

    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    plot(Xcyc, Zcyc, color='k', linestyle='dotted') 
    plot(Xele, Zele, color='k', linestyle='dotted')  
    pc = imshow( 0.5*(Jz[indeces[0]:indeces[1],Ny/2+1,indeces[4]:indeces[5]].T+Jz[indeces[0]:indeces[1],-(Ny/2+2),indeces[4]:indeces[5]].T)/J0[3], aspect='equal', interpolation="nearest", origin='lower', vmin=J_min, vmax=J_max, extent=theextent[1], cmap=thecmap )
    ax = plt.gca()
    xlabel(r'$x \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$z \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_xticks, fontsize=12)
#    plt.yticks(my_zticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{J_z} \; [ \mathrm{ J_{0} } ]$', fontsize=18)
#   title(r'$\mathrm{E_z} \; [ \mathrm{V \, m^{-1}} ]$', fontsize=18)
    text(1.1, 0.65, r'$\uparrow \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\rightarrow \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\otimes \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks_J)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/Jz_XZ_' + fields_iteration + '.'   
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
    cla()
    clf()
    close('all')

    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    pc = imshow( 0.5*(Jz[indeces[0]:indeces[1],indeces[2]:indeces[3],Nz/2+1].T+Jz[indeces[0]:indeces[1],indeces[2]:indeces[3],-(Nz/2+2)].T)/J0[3], aspect='equal', interpolation="nearest", origin='lower', vmin=J_min, vmax=J_max, extent=theextent[2], cmap=thecmap )
    ax = plt.gca()
    xlabel(r'$x \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$y \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_xticks, fontsize=12)
#    plt.yticks(my_yticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{J_z} \; [ \mathrm{ J_{0} } ]$', fontsize=18)
#   title(r'$\mathrm{E_z} \; [ \mathrm{V \, m^{-1}} ]$', fontsize=18)
    text(1.1, 0.65, r'$\odot \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\rightarrow \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\uparrow \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks_J)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')

    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/Jz_XY_' + fields_iteration + '.'   
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
    cla()
    clf()
    close('all')

# CURRENT RADIAL
    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    pc = imshow( 0.5*(Jr[Nx/2+1,indeces[2]:indeces[3],indeces[4]:indeces[5]].T+Jr[-(Nx/2+2),indeces[2]:indeces[3],indeces[4]:indeces[5]].T)/J0[3], aspect='equal', interpolation="nearest", origin='lower', vmin=J_min, vmax=J_max, extent=theextent[0], cmap=thecmap )
    ax = plt.gca()
    xlabel(r'$y \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$z \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_yticks, fontsize=12)
#    plt.yticks(my_zticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{J_r} \; [ \mathrm{ J_{0} } ]$', fontsize=18)
#   title(r'$\mathrm{E_r} \; [ \mathrm{V \, m^{-1}} ]$', fontsize=18)
    text(1.1, 0.65, r'$\uparrow \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\odot \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\rightarrow \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks_J)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/Jr_YZ_' + fields_iteration + '.'   
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
    cla()
    clf()
    close('all')

    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    plot(Xcyc, Zcyc, color='k', linestyle='dotted') 
    plot(Xele, Zele, color='k', linestyle='dotted')  
    pc = imshow( 0.5*(Jr[indeces[0]:indeces[1],Ny/2+1,indeces[4]:indeces[5]].T+Jr[indeces[0]:indeces[1],-(Ny/2+2),indeces[4]:indeces[5]].T)/J0[3], aspect='equal', interpolation="nearest", origin='lower', vmin=J_min, vmax=J_max, extent=theextent[1], cmap=thecmap )
    ax = plt.gca()
    xlabel(r'$x \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$z \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_xticks, fontsize=12)
#    plt.yticks(my_zticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{J_r} \; [ \mathrm{ J_{0} } ]$', fontsize=18)
#   title(r'$\mathrm{E_r} \; [ \mathrm{V \, m^{-1}} ]$', fontsize=18)
    text(1.1, 0.65, r'$\uparrow \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\rightarrow \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\otimes \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks_J)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/Jr_XZ_' + fields_iteration + '.'   
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
    cla()
    clf()
    close('all')

    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    pc = imshow( 0.5*(Jr[indeces[0]:indeces[1],indeces[2]:indeces[3],Nz/2+1].T+Jr[indeces[0]:indeces[1],indeces[2]:indeces[3],-(Nz/2+2)].T)/J0[3], aspect='equal', interpolation="nearest", origin='lower', vmin=J_min, vmax=J_max, extent=theextent[2], cmap=thecmap )
    ax = plt.gca()
    xlabel(r'$x \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$y \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_xticks, fontsize=12)
#    plt.yticks(my_yticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{J_r} \; [ \mathrm{ J_{0} } ]$', fontsize=18)
#   title(r'$\mathrm{E_r} \; [ \mathrm{V \, m^{-1}} ]$', fontsize=18)
    text(1.1, 0.65, r'$\odot \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\rightarrow \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\uparrow \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks_J)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')

    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/Jr_XY_' + fields_iteration + '.'   
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
    cla()
    clf()
    close('all')

    thecmap         = cm.jet

# CURRENT MAG
    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    pc = imshow( 0.5*(Jmag[Nx/2+1,indeces[2]:indeces[3],indeces[4]:indeces[5]].T+Jmag[-(Nx/2+2),indeces[2]:indeces[3],indeces[4]:indeces[5]].T)/J0[3], aspect='equal', interpolation="nearest", origin='lower', vmin=0.0, vmax=J_max-J_min, extent=theextent[0], cmap=thecmap )
    ax = plt.gca()
    xlabel(r'$y \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$z \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_yticks, fontsize=12)
#    plt.yticks(my_zticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{J} \; [ \mathrm{ J_{0} } ]$', fontsize=18)
#   title(r'$\mathrm{E} \; [ \mathrm{V \, m^{-1}} ]$', fontsize=18)
    text(1.1, 0.65, r'$\uparrow \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\odot \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\rightarrow \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks_J)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/Jmag_YZ_' + fields_iteration + '.'   
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
    cla()
    clf()
    close('all')

    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    plot(Xcyc, Zcyc, color='k', linestyle='dotted') 
    plot(Xele, Zele, color='k', linestyle='dotted')  
    pc = imshow( 0.5*(Jmag[indeces[0]:indeces[1],Ny/2+1,indeces[4]:indeces[5]].T+Jmag[indeces[0]:indeces[1],-(Ny/2+2),indeces[4]:indeces[5]].T)/J0[3], aspect='equal', interpolation="nearest", origin='lower', vmin=0.0, vmax=J_max-J_min, extent=theextent[1], cmap=thecmap )
    ax = plt.gca()
    xlabel(r'$x \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$z \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_xticks, fontsize=12)
#    plt.yticks(my_zticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{J} \; [ \mathrm{ J_{0} } ]$', fontsize=18)
#   title(r'$\mathrm{E} \; [ \mathrm{V \, m^{-1}} ]$', fontsize=18)
    text(1.1, 0.65, r'$\uparrow \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\rightarrow \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\otimes \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks_J)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/Jmag_XZ_' + fields_iteration + '.'   
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
    cla()
    clf()
    close('all')

    rc('text', usetex=True)
    gs = gridspec.GridSpec(2, 2, height_ratios=(1,8), width_ratios=(1, 1) )
    plt.figure(figsize=(7.0, 7.0))

    axes = plt.subplot(gs[1,:])
    vlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    hlines(0, -0.5*dxyz[0], 0.5*dxyz[0], color='k')
    pc = imshow( 0.5*(Jmag[indeces[0]:indeces[1],indeces[2]:indeces[3],Nz/2+1].T+Jmag[indeces[0]:indeces[1],indeces[2]:indeces[3],-(Nz/2+2)].T)/J0[3], aspect='equal', interpolation="nearest", origin='lower', vmin=0.0, vmax=J_max-J_min, extent=theextent[2], cmap=thecmap )
    ax = plt.gca()
    xlabel(r'$x \; [ \mathrm{r_{e}} ]$', fontsize=18)
    ylabel(r'$y \; [ \mathrm{r_{e}} ]$', fontsize=18)
#    plt.xticks(my_xticks, fontsize=12)
#    plt.yticks(my_yticks, fontsize=12)
    plt.minorticks_on()
    ax.tick_params(axis='both',which='both',color='white')
    title(r'$\mathrm{J} \; [ \mathrm{ J_{0} } ]$', fontsize=18)
#   title(r'$\mathrm{E} \; [ \mathrm{V \, m^{-1}} ]$', fontsize=18)
    text(1.1, 0.65, r'$\odot \, \mathbf{E_0}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.5, r'$\rightarrow \, \mathrm{Sun}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')
    text(1.1, 0.35, r'$\uparrow \, \mathbf{B}$', horizontalalignment=u'center', verticalalignment=u'center', transform = ax.transAxes, fontsize=18, color='k')

    axes = plt.subplot(gs[0,:])
#    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=LogLocator(subs=range(10)), norm=matplotlib.colors.LogNorm());
    cbar = plt.colorbar(pc, cax=axes, orientation='horizontal', ticks=my_cticks_J)
    cbar.set_label(cmap_title, rotation=0, fontsize = 18, verticalalignment=u'top')

    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.minorticks_on()


    rc('text', usetex=False)
    saveloc   = save_location + '/Jmag_XY_' + fields_iteration + '.'   
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
    cla()
    clf()
    close('all')


# ENDING
    cla()
    clf()
    close('all')
    rc('text', usetex=False)
    print 'END OF PROGRAM';
