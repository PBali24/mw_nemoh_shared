"""
nemoh_utilities
Created on Thu Feb  5 09:15:11 2015
This module is a collection of several functions
used to run the Nemoh software package in Windows
@author: tverbrug
@modified: gverao&pbalitsky on 23.10.17
"""
import os
import shutil as sh
import platform as pt
import numpy as np
import subprocess
from scipy import interpolate #necesary module for interpolating
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
# Used Functions

def calcM(rho=1025.0,sel=2):
    pathFile = './mesh/KH.dat'
    with open(pathFile,'r') as f:
        khRaw = f.readlines()
        goodLine = khRaw[sel]
        khAll = goodLine.split()
        c = float(khAll[sel])

    if sel < 3:
        pathFile = './mesh/Hydrostatics.dat'
        with open(pathFile,'r') as f:
            mRaw = f.readlines()
            goodLine = mRaw[3]
            mAll = goodLine.split()
            mass = rho*float(mAll[2])
    else:
        pathFile = './mesh/Inertia_hull.dat'
        with open(pathFile,'r') as f:
            mRaw = f.readlines()
            goodLine = mRaw[sel-3]
            mAll = goodLine.split()
            mass = float(mAll[sel-3])
        
    return (mass,c)

"""def getAB(dof,nbody=1,sel=10):
    # Selected DOF
    if sel<6:
        selDof = sel
    else:
        selDof = dof.index(1)
    # Open Radiation Coefficients    
    pathFile = './results/RadiationCoefficients.tec'
    with open(pathFile,'r') as f:
        irfRaw = f.readlines()
    strPat = 'body    1 in DoF   ' + str(sum(dof[0:selDof+1]))
    irfInd = 1 + sum(dof[0:selDof])*2
    for iL in range(0,len(irfRaw)):
        irfInfo = irfRaw[iL]
        test = irfInfo.find(strPat)
        if test!=-1:
            indStart = iL
            irfInfo = irfInfo.split()
            irfNrlong = irfInfo[8]
            indComma = irfNrlong.index(',')
            irfNr = int(irfNrlong[0:indComma])
    # Arrange data in RadiationCoefficients file
    A = np.zeros(irfNr)
    B = np.zeros(irfNr)
    omega = np.zeros(irfNr)
    for iL in range(indStart+1,irfNr+indStart+1):
        irfLine = irfRaw[iL].split()
        A[iL-indStart-1] = float(irfLine[irfInd])
        B[iL-indStart-1] = float(irfLine[irfInd+1])
        omega[iL-indStart-1] = float(irfLine[0])
    # Return values
    return(A,B,omega) """

"""def getFe(dof,sel=10):
    # Selected DOF
    if sel<6:
        selDof = sel
    else:
        selDof = dof.index(1)
    # Open Radiation Coefficients    
    pathFile = './results/ExcitationForce.tec'
    with open(pathFile,'r') as f:
        irfRaw = f.readlines()
    strPat = 'I='
    irfInd = 1 + sum(dof[0:selDof])*2
    for iL in range(0,len(irfRaw)):
        irfInfo = irfRaw[iL]
        test = irfInfo.find(strPat)
        if test!=-1:
            indStart = iL
            irfInfo = irfInfo.split()
            irfNrlong = irfInfo[8]
            indComma = irfNrlong.index(',')
            irfNr = int(irfNrlong[0:indComma])
    # Arrange data in RadiationCoefficients file
    FeMod = np.zeros(irfNr)
    FeAng = np.zeros(irfNr)
    for iL in range(indStart+1,irfNr+indStart+1):
        irfLine = irfRaw[iL].split()
        FeMod[iL-indStart-1] = float(irfLine[irfInd])
        FeAng[iL-indStart-1] = abs(float(irfLine[irfInd+1]))
    # Return values
    return(FeMod,FeAng)"""
    
def wavenum(omega,depth):
    k0 = omega**2/9.81
    k_prev = k0    
    test = 10.0
    while(test > 0.0001):
        k_new = k0/np.tanh(k_prev*depth)       
        test = np.abs(k_new-k_prev)
        k_prev = k_new
    return k_new             

def createMeshOpt(zG,nPanels,nsym,rho=1025.0,g=9.81,nbody=1,xG=0.0,yG=0.0):
    if nbody==1:
        fid = open('./Mesh.cal','w')
        fid.write('axisym\n')
        fid.write('{:d}\n'.format(nsym))
        fid.write('0. 0.\n')
        value = '0. 0. {0:f} \n'.format(zG)
        fid.write(value)
        fid.write(str(nPanels) + '\n')
        fid.write('2\n0.\n1.\n')
        fid.write('{0:f}\n'.format(rho))
        fid.write('{0:f}\n'.format(9.81))
        fid.close()
        os.system('Mesh.exe')
    else:
        for iB in range(nbody):
            fid = open('./Mesh.cal','w')
            fid.write('axisym{:d}\n'.format(iB+1))
            fid.write('{:d}\n'.format(nsym))
            fid.write('0. 0.\n')
            value = '{0:f} {1:f} {2:f} \n'.format(xG[iB],yG[iB],zG)
            fid.write(value)
            fid.write(str(nPanels) + '\n')
            fid.write('2\n0.\n1.\n')
            fid.write('{0:f}\n'.format(rho))
            fid.write('{0:f}\n'.format(9.81))
            fid.close()
            os.system('Mesh.exe')


def writeCalFile(rhoW,depW,omega,zG,dof,aO={},nbody=1,xG=[0.0],yG=[0.0]):
    # Read info on the mesh
    nrNode = [0]*nbody
    nrPanel = [0]*nbody
    for iB in range(nbody):
        if nbody == 1:
            f1 = open('./mesh/axisym_info.dat','r')
        else:
            f1 = open('./mesh/axisym{:d}_info.dat'.format(iB+1),'r')
        lineInfo = f1.readline()
        lineInfo = lineInfo.split()
        nrNode[iB] = int(lineInfo[0])
        nrPanel[iB] = int(lineInfo[1])
        f1.close()
    # Read advanced options if there are any
    dirCheck = aO['dirCheck']
    irfCheck = aO['irfCheck']
    kochCheck = aO['kochCheck']
    fsCheck = aO['fsCheck']

    # Create the Nemoh calibration file
    fid = open('./Nemoh.cal','w')
    fid.write('--- Environment ---\n')
    fid.write(str(rhoW) + '				! RHO 			! KG/M**3 	! Fluid specific volume\n')
    fid.write('9.81				! G			! M/S**2	! Gravity\n')
    fid.write(str(depW) + '				! DEPTH			! M		! Water depth\n')
    fid.write('0.	0.			! XEFF YEFF		! M		! Wave measurement point\n')
    fid.write('--- Description of floating bodies ---\n')
    fid.write('{:d}				! Number of bodies\n'.format(nbody))
    for iB in range(nbody):      
        fid.write('--- Body {:d} ---\n'.format(iB+1))
        if nbody == 1:
            fid.write('axisym.dat			! Name of mesh file\n')
        else:
            fid.write('axisym{:d}.dat			! Name of mesh file\n'.format(iB+1))
        fid.write(str(nrNode[iB]) + '\t' + str(nrPanel[iB]) + '			! Number of points and number of panels\n')
        fid.write('{:d}				! Number of degrees of freedom\n'.format(sum(dof)))
        for iDof in range(len(dof)):
            if (iDof == 0 and dof[iDof] == 1):
                fid.write('1 1. 0.	0. 0. 0. 0.		! Surge\n')
            elif (iDof == 1 and dof[iDof] == 1):
                fid.write('1 0. 1.	0. 0. 0. 0.		! Sway\n')
            elif (iDof == 2 and dof[iDof] == 1):
                fid.write('1 0. 0. 1. 0. 0. 0.		! Heave\n')
            elif (iDof == 3 and dof[iDof] == 1):
                fid.write('2 1. 0. 0. {0:f} {1:f} {2:f}		! Roll about CdG\n'.format(xG[iB],yG[iB],zG))
            elif (iDof == 4 and dof[iDof] == 1):
                fid.write('2 0. 1. 0. {0:f} {1:f} {2:f}		! Pitch about CdG\n'.format(xG[iB],yG[iB],zG))
            elif (iDof == 5 and dof[iDof] == 1):
                fid.write('2 0. 0. 1. {0:f} {1:f} {2:f}		! Yaw about CdG\n'.format(xG[iB],yG[iB],zG))
        fid.write('{:d}				! Number of resulting generalised forces\n'.format(sum(dof)))
        for iDof in range(len(dof)):
            if (iDof == 0 and dof[iDof] == 1):
                fid.write('1 1. 0.	0. 0. 0. 0.		! Force in X direction\n')
            elif (iDof == 1 and dof[iDof] == 1):
                fid.write('1 0. 1.	0. 0. 0. 0.		! Force in Y direction\n')
            elif (iDof == 2 and dof[iDof] == 1):
                fid.write('1 0. 0. 1. 0. 0. 0.		! Force in Z direction\n')
            elif (iDof == 3 and dof[iDof] == 1):
                fid.write('2 1. 0. 0. {0:f} 0. {1:f}		! Roll Moment about CdG\n'.format(xG[iB],zG))
            elif (iDof == 4 and dof[iDof] == 1):
                fid.write('2 0. 1. 0. {0:f} 0. {1:f}		! Pitch Moment about CdG\n'.format(xG[iB],zG))
            elif (iDof == 5 and dof[iDof] == 1):
                fid.write('2 0. 0. 1. {0:f} 0. {1:f}		! Yaw Moment about CdG\n'.format(xG[iB],zG))
        fid.write('0				! Number of lines of additional information\n')
    
    fid.write('--- Load cases to be solved ---\n')
    fid.write('{0:d}\t{1:.12f}\t{2:.12f}		! Number of wave frequencies, Min, and Max (rad/s)\n'.format(omega[0],omega[1],omega[2]))
    if dirCheck:
        fid.write(str(aO['dirStep']) + '\t' + str(aO['dirStart']) + '\t' + str(aO['dirStop']) + '		! Number of wave directions, Min and Max (degrees)\n')
    else:
        fid.write('1	0.	0.		! Number of wave directions, Min and Max (degrees)\n')
    fid.write('--- Post processing ---\n')
    if irfCheck:
        fid.write('1' + '\t' + str(aO['irfStep']) + '\t' + str(aO['irfDur']) + '\t\t! IRF 				! IRF calculation (0 for no calculation), time step and duration\n')
    else:
        fid.write('0' + '\t0.01\t20.\t\t! IRF 				! IRF calculation (0 for no calculation), time step and duration\n')
    fid.write('0				! Show pressure\n')
    if kochCheck:
        fid.write(str(aO['kochStep']) + '\t' + str(aO['kochStart']) + '\t' + str(aO['kochStop']) + '		! Kochin function 		! Number of directions of calculation (0 for no calculations), Min and Max (degrees)\n')
    else:
        fid.write('0	0.	180.		! Kochin function 		! Number of directions of calculation (0 for no calculations), Min and Max (degrees)\n')
    if fsCheck:
        fid.write(str(aO['fsNx']) + '\t' + str(aO['fsNy']) + '\t' + str(aO['fsLengthX']) + '\t' + str(aO['fsLengthY']) + '	! Free surface elevation 	! Number of points in x direction (0 for no calcutions) and y direction and dimensions of domain in x and y direction	\n')
    else:
        fid.write('0	2	1000.	2.	! Free surface elevation 	! Number of points in x direction (0 for no calcutions) and y direction and dimensions of domain in x and y direction	\n')
    
    fid.close()

def readNemohResults(dof,dirname):
    # Open index file
    with open(dirname+'/index.dat','r') as fid:
        allData = fid.readlines()
        
    Nw,Nbeta,Nradiation,Nintegration,Ntheta = [int(a) for a in allData[0].split()]
    w = np.zeros(Nw)
    beta = np.zeros(Nbeta)
    beta[:] = [float(a) for a in allData[Nintegration+Nradiation+3].split()]
    w[:] = [float(a) for a in allData[Nintegration+Nradiation+4].split()]
    
    with open(dirname+'/ExcitationForce.tec','r') as fid:
        allData = fid.readlines()

    Famp = np.zeros((Nw,Nbeta,Nintegration))
    Fphi = np.zeros((Nw,Nbeta,Nintegration))

    pp=0
    for ll in range(Nbeta):
        for kk in range(Nw): #it works fine for just one frequency and one DOF, consider moving to more frequencies saving results as an estructure array
                #data is saved in the order of the DOF followed in all the
                #scripts
            line = [float(a) for a in allData[pp+2+Nradiation].split()]#is nintegration or nradiation?
            c=1            
            for iB in range(Nintegration):
                    Famp[kk,ll,iB]=line[c]  #Read the Amplitude of the Fe, in each line of ExcitationForce.tec the even element
                    Fphi[kk,ll,iB]=line[c+1]#Read the Phase of the Fe, in each line of ExcitationForce.tec the odd element (excluding the 1)
                    c += 2
            pp += 1
        pp += 1

    ExcitationForce=Famp*(np.cos(Fphi)+1j*np.sin(Fphi))
     
    with open(dirname+'/RadiationCoefficients.tec','r') as fid:
            allDataAB = fid.readlines()
            
    ARaw=np.zeros((Nw,Nradiation,Nintegration))#this will work only for one frequency (check more than one frequency)
    BRaw=np.zeros((Nw,Nradiation,Nintegration))#this will work only for one frequency (check more than one frequency)
    
    pp=0
    for iB in range(Nintegration):
        for iw in range(Nw):
            line = [float(a) for a in allDataAB[pp+2+Nradiation].split()]#is nintegration or nradiation?
            c = 1
            for id in range(Nradiation):
                ARaw[iw,id,iB]=line[c]  #Read the Amplitude of the Fe, in each line of ExcitationForce.tec the even element
                BRaw[iw,id,iB]=line[c+1]#Read the Phase of the Fe, in each line of ExcitationForce.tec the odd element (excluding the 1)
                c += 2
            pp += 1
        pp += 1           
                            
        fid.close()
                
        # Reshape added mass and hydrodynamic damping
        ndof = np.sum(dof)
        Ma = np.zeros((Nw,ndof,Nintegration,Nintegration))
        Bhyd = np.zeros((Nw,ndof,Nintegration,Nintegration))
        for iW in range(Nw):
            for iDof in range(ndof):
                Ma[iW,iDof,:,:] = ARaw[iW,iDof:Nradiation:ndof,:]
                Bhyd[iW,iDof,:,:] = BRaw[iW,iDof:Nradiation:ndof,:]
            
    return w, beta, Ma, Bhyd, ExcitationForce

def runNemoh(show_console,nbody=1):
    if nbody == 1:
        sh.copyfile('./mesh/axisym.dat','axisym.dat')
    else:
        for iB in range(nbody):
            sh.copyfile('./mesh/axisym{:d}.dat'.format(iB+1),'axisym{:d}.dat'.format(iB+1))
    if pt.system()=='Linux':
        subprocess.call('./preProc')
        subprocess.call('./solver')
        os.system('./postProc')
    elif show_console:
        os.system('preProcessor.exe')
        os.system('Solver.exe')
        os.system('postProcessor.exe')
    else: 
        subprocess.call('preProcessor.exe')
        subprocess.check_call('Solver.exe')
        subprocess.check_call('postProcessor.exe')

def wavlen(H,T,d):

    g = 9.81
    dpi = 2*np.pi
    
    L0 = g*T**2/dpi
    L1 = g*T**2/dpi*np.tanh(dpi*d/L0)
    
    while (np.abs(L1-L0) > 0.001):
        L0 = L1
        L1 = g*T**2/dpi*np.tanh(dpi*d/L0)
        
    return L1

def nemoh_interp(NEM_ini,NEM_GRID,rad_wave,diff_wave,pert_wave,inc_wave,tot_wave,X,Y):
    #get rid of number of frequncies only use nf = 1 
#-----------------------------------------------------------------------------
#                              0. LOAD NEMOH DATA
#-----------------------------------------------------------------------------
    #-----------------------------------------------------------------------------
    #                              0. CREATE TARGET GRID
    #-----------------------------------------------------------------------------

    dxN=NEM_GRID['dxN']#Target delta for the interpolation. Introduce desired grid dx
    dyN=NEM_GRID['dyN']#Target delta for the interpolation. Introduce desired grid dy
    xwidth=NEM_GRID['Lg'] #Value of the lenght of the NEMOH domain
    ylength=NEM_GRID['Wg']#Value of the width fo the NEMOH
    Nx = int(np.round(NEM_GRID['Lg']/dxN) + 1) #redifine number of gridpoints as Nxg = np.round(Lg/dx) + 1
    Ny = int(np.round(NEM_GRID['Wg']/dyN) + 1) #redifine number of gridpoints as Nyg = np.round(Wg/dy) + 1
    xvectN = np.arange(-xwidth/2,xwidth/2+1,dxN)
    xvectN = xvectN[:Nx]
    yvectN = np.arange(-ylength/2,ylength/2+1,dyN)
    yvectN = yvectN[:Ny]
    [xvectN,yvectN]=np.meshgrid(xvectN,yvectN)
    #-----------------------------------------------------------------------------
    #                              0. INTERPOLATE 
    #-----------------------------------------------------------------------------
    #-----------------------------------------------------------------------------
    #                              0. RAD_WAVE 
    #-----------------------------------------------------------------------------
    #AFTER INTERPOLATING AT LEAST FOR PLOTTING THE MATRIX NEEDS TO BE TRANSPOSED ???/!!!!!
    #the module and the angle of the complex wave field has to be interpolated independently

    diff_waveN = np.zeros((Ny,Nx),dtype=np.complex_)
    rad_waveN = np.zeros((Ny,Nx),dtype=np.complex_)
    tot_waveN = np.zeros((Ny,Nx),dtype=np.complex_)
    pert_waveN = np.zeros((Ny,Nx),dtype=np.complex_)
    inc_waveN= np.zeros((Ny,Nx),dtype=np.complex_) 
 
    rad_mod=np.array(np.abs(rad_wave[:]))
    rad_pha=np.array(np.angle(rad_wave[:]))
    #Alternative loop version for interpolating GVF
    fin = interpolate.RegularGridInterpolator((Y[:,0],X[0,:]),rad_mod[:,:])
    
    '''rad_modN = fin(yvect_mw[:, 0], xvect_mw[0, :])
    rad_modN = np.zeros((np.shape(yvect_mw)[0],np.shape(xvect_mw)[0]))
    for jj in range(np.shape(yvect_mw)[0]):
        for ii in range(np.shape(xvect_mw)[0]):
            yaux=yvect_mw[jj,0]
            xaux=xvect_mw[0,ii]
            rad_modN[jj,ii] = fin(np.array([yaux,xaux]))'''
    
    #Interpolate module of the complex number
    pointlist = [[a,b] for a in yvectN[:,0] for b in xvectN[0,:]]
    rad_modN = np.reshape(fin(np.array(pointlist)),(Ny,Nx))
    mask = np.isnan(rad_modN)
    rad_modN[mask] = 0.0
    #rad_modN = np.transpose(rad_modN)
    
    #Interpolate angle of the complex number
    fin = interpolate.RegularGridInterpolator((Y[:,0],X[0,:]),rad_pha[:,:])
    pointlist = [[a,b] for a in yvectN[:,0] for b in xvectN[0,:]]
    rad_phaN = np.reshape(fin(np.array(pointlist)),(Ny,Nx))
    mask = np.isnan(rad_phaN)
    rad_phaN[mask] = 0.0
    #rad_phaN = np.transpose(rad_phaN)
    
    #Obtain the complex wave field again
    
    #rad_waveN=np.zeros(np.shape(yvect_mw)[0],np.shape(xvect_mw[1]))
    #rad_waveN=rad_modN+1j*rad_phaN #Alternatively the real and imaginary part could be interpoalted and use this equation
    rad_waveN_aux = rad_modN * np.exp( 1j * rad_phaN)
    rad_waveN[:,:] = rad_waveN_aux
    #-----------------------------------------------------------------------------
    #                              0. DIFF_WAVE 
    #-----------------------------------------------------------------------------
    #AFTER INTERPOLATING AT LEAST FOR PLOTTING THE MATRIX NEEDS TO BE TRANSPOSED ???/!!!!!
    #the module and the angle of the complex wave field has to be interpolated independently

    diff_mod=np.array(np.abs(diff_wave))
    diff_pha=np.array(np.angle(diff_wave))
    #Alternative loop version for interpolating GVF
    fin = interpolate.RegularGridInterpolator((Y[:,0],X[0,:]),diff_mod[:,:])
    
    '''rad_modN = fin(yvect_mw[:, 0], xvect_mw[0, :])
    rad_modN = np.zeros((np.shape(yvect_mw)[0],np.shape(xvect_mw)[0]))
    for jj in range(np.shape(yvect_mw)[0]):
        for ii in range(np.shape(xvect_mw)[0]):
            yaux=yvect_mw[jj,0]
            xaux=xvect_mw[0,ii]
            rad_modN[jj,ii] = fin(np.array([yaux,xaux]))'''
    
    #Interpolate module of the complex number
    pointlist = [[a,b] for a in yvectN[:,0] for b in xvectN[0,:]]
    diff_modN = np.reshape(fin(np.array(pointlist)),(Ny,Nx))
    mask = np.isnan(diff_modN)
    diff_modN[mask] = 0.0
    #diff_modN = np.transpose(diff_modN)
    
    #Interpolate angle of the complex number
    fin = interpolate.RegularGridInterpolator((Y[:,0],X[0,:]),diff_pha[:,:])
    pointlist = [[a,b] for a in yvectN[:,0] for b in xvectN[0,:]]
    diff_phaN = np.reshape(fin(np.array(pointlist)),(Ny,Nx))
    mask = np.isnan(diff_phaN)
    diff_phaN[mask] = 0.0
    #diff_phaN = np.transpose(diff_phaN)
    
    #Obtain the complex wave field again
    
    #diff_waveN=np.zeros(np.shape(yvect_mw)[0],np.shape(xvect_mw[1]))
    #diff_waveN=diff_modN+1j*diff_phaN #Alternatively the real and imaginary part could be interpoalted and use this equation
    diff_waveN_aux = diff_modN * np.exp( 1j * diff_phaN)
    diff_waveN[:,:] = diff_waveN_aux
    #-----------------------------------------------------------------------------
    #                              0. INC_WAVE 
    #-----------------------------------------------------------------------------
    #AFTER INTERPOLATING AT LEAST FOR PLOTTING THE MATRIX NEEDS TO BE TRANSPOSED ???/!!!!!
    #the module and the angle of the complex wave field has to be interpolated independently
    
    inc_mod=np.array(np.abs(inc_wave))
    inc_pha=np.array(np.angle(inc_wave))
    #Alternative loop version for interpolating GVF
    fin = interpolate.RegularGridInterpolator((Y[:,0],X[0,:]),inc_mod[:,:])
    
    '''rad_modN = fin(yvect_mw[:, 0], xvect_mw[0, :])
    rad_modN = np.zeros((np.shape(yvect_mw)[0],np.shape(xvect_mw)[0]))
    for jj in range(np.shape(yvect_mw)[0]):
        for ii in range(np.shape(xvect_mw)[0]):
            yaux=yvect_mw[jj,0]
            xaux=xvect_mw[0,ii]
            rad_modN[jj,ii] = fin(np.array([yaux,xaux]))'''
    
    #Interpolate module of the complex number
    pointlist = [[a,b] for a in yvectN[:,0] for b in xvectN[0,:]]
    inc_modN = np.reshape(fin(np.array(pointlist)),(Ny,Nx))
    mask = np.isnan(inc_modN)
    inc_modN[mask] = 0.0
    #inc_modN = np.transpose(inc_modN)
    
    #Interpolate angle of the complex number
    fin = interpolate.RegularGridInterpolator((Y[:,0],X[0,:]),inc_pha[:,:])
    pointlist = [[a,b] for a in yvectN[:,0] for b in xvectN[0,:]]
    inc_phaN = np.reshape(fin(np.array(pointlist)),(Ny,Nx))
    mask = np.isnan(inc_phaN)
    inc_phaN[mask] = 0.0
    #inc_phaN = np.transpose(inc_phaN)
    
    #Obtain the complex wave field again
    
    #rad_waveN=np.zeros(np.shape(yvect_mw)[0],np.shape(xvect_mw[1]))
    #rad_waveN=rad_modN+1j*rad_phaN #Alternatively the real and imaginary part could be interpoalted and use this equation
    inc_waveN_aux = inc_modN * np.exp( 1j * inc_phaN)
    inc_waveN[:,:] = inc_waveN_aux
    #-----------------------------------------------------------------------------
    #                              0. PERT_WAVE 
    #-----------------------------------------------------------------------------
    #pert_waveN = rad_waveN + diff_waveN
    pert_mod=np.array(np.abs(pert_wave))
    pert_pha=np.array(np.angle(pert_wave))
    #Alternative loop version for interpolating GVF
    fin = interpolate.RegularGridInterpolator((Y[:,0],X[0,:]),pert_mod[:,:])
    
    '''rad_modN = fin(yvect_mw[:, 0], xvect_mw[0, :])
    rad_modN = np.zeros((np.shape(yvect_mw)[0],np.shape(xvect_mw)[0]))
    for jj in range(np.shape(yvect_mw)[0]):
        for ii in range(np.shape(xvect_mw)[0]):
            yaux=yvect_mw[jj,0]
            xaux=xvect_mw[0,ii]
            rad_modN[jj,ii] = fin(np.array([yaux,xaux]))'''
    
    #Interpolate module of the complex number
    pointlist = [[a,b] for a in yvectN[:,0] for b in xvectN[0,:]]
    pert_modN = np.reshape(fin(np.array(pointlist)),(Ny,Nx))
    mask = np.isnan(pert_modN)
    pert_modN[mask] = 0.0
    #pert_modN = np.transpose(pert_modN)
    
    #Interpolate angle of the complex number
    fin = interpolate.RegularGridInterpolator((Y[:,0],X[0,:]),pert_pha[:,:])
    pointlist = [[a,b] for a in yvectN[:,0] for b in xvectN[0,:]]
    pert_phaN = np.reshape(fin(np.array(pointlist)),(Ny,Nx))
    mask = np.isnan(pert_phaN)
    pert_phaN[mask] = 0.0
    #pert_phaN = np.transpose(pert_phaN)
    
    #Obtain the complex wave field again
    
    #rad_waveN=np.zeros(np.shape(yvect_mw)[0],np.shape(xvect_mw[1]))
    #rad_waveN=rad_modN+1j*rad_phaN #Alternatively the real and imaginary part could be interpoalted and use this equation
    pert_waveN_aux = pert_modN * np.exp( 1j * pert_phaN)
    pert_waveN[:,:] = pert_waveN_aux
    #-----------------------------------------------------------------------------
    #                              0. TOT_WAVE 
    #-----------------------------------------------------------------------------
    #tot_waveN = inc_waveN + pert_waveN
    tot_mod=np.array(np.abs(tot_wave))
    tot_pha=np.array(np.angle(tot_wave))
    #Alternative loop version for interpolating GVF
    fin = interpolate.RegularGridInterpolator((Y[:,0],X[0,:]),tot_mod[:,:])
    
    '''rad_modN = fin(yvect_mw[:, 0], xvect_mw[0, :])
    rad_modN = np.zeros((np.shape(yvect_mw)[0],np.shape(xvect_mw)[0]))
    for jj in range(np.shape(yvect_mw)[0]):
        for ii in range(np.shape(xvect_mw)[0]):
            yaux=yvect_mw[jj,0]
            xaux=xvect_mw[0,ii]
            rad_modN[jj,ii] = fin(np.array([yaux,xaux]))'''
    
    #Interpolate module of the complex number
    pointlist = [[a,b] for a in yvectN[:,0] for b in xvectN[0,:]]
    tot_modN = np.reshape(fin(np.array(pointlist)),(Ny,Nx))
    mask = np.isnan(tot_modN)
    tot_modN[mask] = 0.0
    #tot_modN = np.transpose(tot_modN)
    
    #Interpolate angle of the complex number
    fin = interpolate.RegularGridInterpolator((Y[:,0],X[0,:]),tot_pha[:,:])
    pointlist = [[a,b] for a in yvectN[:,0] for b in xvectN[0,:]]
    tot_phaN = np.reshape(fin(np.array(pointlist)),(Ny,Nx))
    mask = np.isnan(tot_phaN)
    tot_phaN[mask] = 0.0
    #tot_phaN = np.transpose(tot_phaN)
    
    #Obtain the complex wave field again
    
    #rad_waveN=np.zeros(np.shape(yvect_mw)[0],np.shape(xvect_mw[1]))
    #rad_waveN=rad_modN+1j*rad_phaN #Alternatively the real and imaginary part could be interpoalted and use this equation
    tot_waveN_aux = tot_modN * np.exp( 1j * tot_phaN)
    tot_waveN[:,:] = tot_waveN_aux

    return rad_waveN, diff_waveN, pert_waveN, inc_waveN, tot_waveN, xvectN, yvectN, Nx, Ny

def plotmesh(nbodies,dir_mesh,xlim,xBody,yBody):
    if nbodies == 1:
        ylim = xlim
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.set_aspect('equal')
        with open(dir_mesh+'/Calculation/mesh/axisym_info.dat','r') as fid:
            C = fid.readlines()
            npoints = int(C[0].split()[0])
        #dir_mesh = r"E:/WEC_WAKES_IRREGULAR_VALIDATION/NEMOH/Calculation"
        for ii in range(nbodies):
            with open(dir_mesh+'/Calculation/mesh/axisym.dat','r') as fid:
                B = fid.readlines()
            points = B[2:npoints-1]
            
            points = [point.split() for point in points]
            [point.pop(0) for point in points]
     
            x = [float(point[0]) for point in points]
            x = [xloc - xBody[0] for xloc in x]
            y = [float(point[1]) for point in points]
            y = [yloc - yBody[0] for yloc in y]
            z = [float(point[2]) for point in points]
        
            ax.scatter(x, y, z,c='b',s = 0.5)
        
            #ax.plot_trisurf(x, y, z, linewidth=0.2, antialiased=True)
            
        ax.set_xlabel('X Label')
        ax.set_ylabel('Y Label')
        ax.set_zlabel('Z Label')
        ax.set_xlim(-70,70)
        ax.set_ylim(-70,70)
        ax.set_zlim(-70,70)
        ax.set_xticks(np.arange(-70, 70+1, 50.0))
        ax.set_yticks(np.arange(-70, 70+1, 50.0))
        ax.set_zticks(np.arange(-70, 70+1, 50.0))
        #ax.view_init(azim=0, elev=90)
        plt.show()
        plt.waitforbuttonpress()
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.set_aspect('equal')
        with open(dir_mesh+'/Calculation/mesh/axisym_info.dat','r') as fid:
            C = fid.readlines()
            npoints = int(C[0].split()[0])
        #dir_mesh = r"E:/WEC_WAKES_IRREGULAR_VALIDATION/NEMOH/Calculation"
        for ii in range(nbodies):
            with open(dir_mesh+'/Calculation/mesh/axisym.dat','r') as fid:
                B = fid.readlines()
            points = B[2:npoints-1]
            points = [point.split() for point in points]
            [point.pop(0) for point in points]
     
            x = [float(point[0]) for point in points]
            y = [float(point[1]) for point in points]
            z = [float(point[2]) for point in points]
        
            ax.scatter(x, y, z,c='b',s = 0.01)
        
            #ax.plot_trisurf(x, y, z, linewidth=0.2, antialiased=True)
            
        ax.set_xlabel('X Label')
        ax.set_ylabel('Y Label')
        ax.set_zlabel('Z Label')
        ax.set_xlim(-xlim,xlim)
        ax.set_ylim(-ylim,ylim)
        ax.set_zlim(-ylim,ylim)
        ax.set_xticks(np.arange(-xlim, xlim+1, 50.0))
        ax.set_yticks(np.arange(-ylim, ylim+1, 50.0))
        ax.set_zticks(np.arange(-ylim, ylim+1, 100.0))
        ax.view_init(azim=0, elev=90)
        plt.show()
        plt.waitforbuttonpress()
    else:
        ylim = xlim
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.set_aspect('equal')
        with open(dir_mesh+'/Calculation/mesh/axisym1_info.dat','r') as fid:
            C = fid.readlines()
            npoints = int(C[0].split()[0])
        #dir_mesh = r"E:/WEC_WAKES_IRREGULAR_VALIDATION/NEMOH/Calculation"
        for ii in range(nbodies):
            with open(dir_mesh+'/Calculation/mesh/axisym1.dat','r') as fid:
                B = fid.readlines()
            points = B[2:npoints-1]
            
            points = [point.split() for point in points]
            [point.pop(0) for point in points]
     
            x = [float(point[0]) for point in points]
            x = [xloc - xBody[0] for xloc in x]
            y = [float(point[1]) for point in points]
            y = [yloc - yBody[0] for yloc in y]
            z = [float(point[2]) for point in points]
        
            ax.scatter(x, y, z,c='b',s = 0.5)
        
            #ax.plot_trisurf(x, y, z, linewidth=0.2, antialiased=True)
            
        ax.set_xlabel('X Label')
        ax.set_ylabel('Y Label')
        ax.set_zlabel('Z Label')
        ax.set_xlim(-70,70)
        ax.set_ylim(-70,70)
        ax.set_zlim(-70,70)
        ax.set_xticks(np.arange(-70, 70+1, 50.0))
        ax.set_yticks(np.arange(-70, 70+1, 50.0))
        ax.set_zticks(np.arange(-70, 70+1, 50.0))
        #ax.view_init(azim=0, elev=90)
        plt.show()
        plt.waitforbuttonpress()
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.set_aspect('equal')
        with open(dir_mesh+'/Calculation/mesh/axisym1_info.dat','r') as fid:
            C = fid.readlines()
            npoints = int(C[0].split()[0])
        #dir_mesh = r"E:/WEC_WAKES_IRREGULAR_VALIDATION/NEMOH/Calculation"
        for ii in range(nbodies):
            with open(dir_mesh+'/Calculation/mesh/axisym'+str(ii+1)+'.dat','r') as fid:
                B = fid.readlines()
            points = B[2:npoints-1]
            
            points = [point.split() for point in points]
            [point.pop(0) for point in points]
     
            x = [float(point[0]) for point in points]
            y = [float(point[1]) for point in points]
            z = [float(point[2]) for point in points]
        
            ax.scatter(x, y, z,c='b',s = 0.01)
        
            #ax.plot_trisurf(x, y, z, linewidth=0.2, antialiased=True)
            
        ax.set_xlabel('X Label')
        ax.set_ylabel('Y Label')
        ax.set_zlabel('Z Label')
        ax.set_xlim(-xlim,xlim)
        ax.set_ylim(-ylim,ylim)
        ax.set_zlim(-ylim,ylim)
        ax.set_xticks(np.arange(-xlim, xlim+1, 50.0))
        ax.set_yticks(np.arange(-ylim, ylim+1, 50.0))
        ax.set_zticks(np.arange(-ylim, ylim+1, 100.0))
        ax.view_init(azim=0, elev=90)
        plt.show()
        plt.waitforbuttonpress()
    return
