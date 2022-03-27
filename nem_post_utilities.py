# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 15:59:10 2018
Original script created by GVF 19/10/2018
Python script containing the functions to obtain the A,B and Fe coeffecients, RAO and surface elevation of a nemoh run
It generates the NEMOH postprocessed files.
"""
import numpy as np
import os
import sys
#put the relative path to your NM_MW_shared folder
sys.path.append(os.path.abspath(os.path.join(os.path.dirname("..//MW_NEM_shared//"), '..')))
from MW_NEM_shared import nem_utilities as ne

#from scipy import signal

#-----------------------------------------------------------------------------
#                              1.4. NEMOH SIMULATION PARAMETERS (ADVANCED OPTIONS)
#-----------------------------------------------------------------------------
def ini(NEM_ini):#This funcion initializes NEM_OUT data structure
    NEM_OUT = {}#Structure containing Ouput Results for the NEMOH simulation.
    NEM_OUT['Fe'] = []   #Exitacion force
    NEM_OUT['Ma'] = []   #Added Mass
    NEM_OUT['Bhyd'] = [] #Hydrodynamic Damping
    NEM_OUT['Rao'] = []  #Rao
    NEM_OUT['rad_vel'] = []
    NEM_OUT['inc_wave'] = []
    NEM_OUT['diff_wave'] = []
    NEM_OUT['rad_wave'] = []
    NEM_OUT['pert_wave'] = []
    NEM_OUT['tot_wave'] = []
    if NEM_ini['interpolation']:
        NEM_OUT['inc_waveN'] = []
        NEM_OUT['diff_waveN'] = []
        NEM_OUT['rad_waveN'] = []
        NEM_OUT['pert_waveN'] = []
        NEM_OUT['tot_waveN'] = []
        
    return(NEM_OUT)

def RAO(w,NEM_ini,NEM_BODY,Fe,Ma,Bhyd):
    #for w in w:#use in case one single NEMOH run is used for more than one frequency TODO not sure it will work
    Mass = np.diag(np.ones(NEM_BODY['nbody']))*NEM_BODY['Mass']
    Kh = np.diag(np.ones(NEM_BODY['nbody']))*NEM_BODY['Kh']
    Mextra = np.diag(np.ones(NEM_BODY['nbody']))*NEM_BODY['ptoProp'][0]
    Bpto = np.diag(np.ones(NEM_BODY['nbody']))*NEM_BODY['ptoProp'][1]
    Kextra = np.diag(np.ones(NEM_BODY['nbody']))*NEM_BODY['ptoProp'][2]
    #Den = w^2 * (Ma + Mass + Mextra) - i*w*(Bhyd+Bpto) + Kh + Kextra Denominator of the RAO equation
    #An alternative calculation would be RAO = np.dot(np.linalg.inv(Den),Fe) = Den^-1 * Fe
    Den = (-(w**2)*(Ma[0,0,:,:].squeeze() + Mass + Mextra) #Implict line continuation between brackets
    - 1j*w*(Bhyd[0,0,:,:].squeeze() + Bpto) + Kh + Kextra)#TODO first 0 is w iterator second 0 is for the dof iterator
    F0 = Fe[0,:,:].squeeze(axis = 0)
    Rao = np.linalg.solve(Den, F0) #Solving the RAO as a Matrix Equation System DEN * RAO = Fe
    #print(np.abs(Rao)) 
    return(Rao)
    
def rad_vel(w,Rao):
    #for w in w:#use in case one single NEMOH run is used for more than one frequency TODO not sure it will work
    rad_vel = -1j*w*Rao
    return(rad_vel)

#modify funcion for arbitrary number of arduments if nterp=true -> return also the the interpolated parameters 

def wave_field(NEM_ini,NEM_GRID,NEM_BODY,rad_vel,NEMOHDir,f,deg,depth):

    dirname = os.path.join(NEMOHDir,'Calculation','results')
    #for w in w:#use in case one single NEMOH run is used for more than one frequency TODO not sure it will work
    diff_wave = np.zeros((NEM_GRID['Ny'],NEM_GRID['Nx']),dtype=np.complex_)
    rad_wave = np.zeros((NEM_GRID['Ny'],NEM_GRID['Nx']),dtype=np.complex_)
    tot_wave = np.zeros((NEM_GRID['Ny'],NEM_GRID['Nx']),dtype=np.complex_)
    pert_wave = np.zeros((NEM_GRID['Ny'],NEM_GRID['Nx']),dtype=np.complex_)
    inc_wave = np.zeros((NEM_GRID['Ny'],NEM_GRID['Nx']),dtype=np.complex_)
    
    xFS = np.arange(-NEM_GRID['Lg']/2.0, NEM_GRID['Lg']/2.0+NEM_GRID['dx'],NEM_GRID['dx'])
    xFS = xFS[:NEM_GRID['Nx']]

    yFS = np.arange(-NEM_GRID['Wg']/2.0, NEM_GRID['Wg']/2.0+NEM_GRID['dy'],NEM_GRID['dy'])
    yFS = yFS[:NEM_GRID['Ny']]

    X,Y = np.meshgrid(xFS,yFS) #X = xvect and Y = yvect
    
    # Read free surface files
    modE = np.zeros((NEM_GRID['Ny'],NEM_GRID['Nx']))
    phaE = np.zeros((NEM_GRID['Ny'],NEM_GRID['Nx']))
    rad_wave_aux = np.zeros((NEM_BODY['nbody'],NEM_GRID['Ny'],NEM_GRID['Nx']),dtype=np.complex_)
    mR = np.zeros((NEM_BODY['nbody'],NEM_GRID['Ny'],NEM_GRID['Nx']))
    pR = np.zeros((NEM_BODY['nbody'],NEM_GRID['Ny'],NEM_GRID['Nx']))
    real_mw = np.zeros((NEM_BODY['nbody'],NEM_GRID['Ny'],NEM_GRID['Nx']),dtype=np.complex_)
    complex_mw = np.zeros((NEM_BODY['nbody'],NEM_GRID['Ny'],NEM_GRID['Nx']),dtype=np.complex_)
    #-----------------------------------------------------------------------------
    #                              4.1. Incident Wave
    #-----------------------------------------------------------------------------        
    Nx=int(NEM_GRID['Nx'])
    Ny=int(NEM_GRID['Ny'])
    xlength=NEM_GRID['Lg']
    ywidth=NEM_GRID['Wg']
    
    Lwav = wavlen(1.0,1/f,depth)
    theta = -deg *np.pi/180
            
    kr = 2.0*np.pi/Lwav
    for jj in range (0,Ny):
        for kk in range (0,Nx):
            xtheta = (-xlength/2 + kk*NEM_GRID['dx']) * np.cos(theta)
            ytheta = (+ywidth/2 - jj*NEM_GRID['dy']) * np.sin(theta) 
            inc_wave[jj,kk] = np.exp(1j*kr*(xtheta+ytheta)) 
    #-----------------------------------------------------------------------------
    #                              4.2. Diffracted Wave
    #-----------------------------------------------------------------------------
    with open(dirname+'/freesurface.    1.dat','r') as fid:
        allEta = fid.readlines()
    iL = 0
    for iX in range(int(NEM_GRID['Nx'])):
        for iY in range(int(NEM_GRID['Ny'])):
            line = [float(a) for a in allEta[2+iL].split()]
            modE[iY,iX] = line[2]   
            phaE[iY,iX] = line[3]
            iL += 2
                
    mD = modE
    pD = phaE
    real_d=mD[:,:]*np.cos(pD[:,:])#real value of the wave_field
    complex_d=mD[:,:]*np.sin(pD[:,:])#complex value of the wave_field
    diff_wave[:,:]=real_d+1j*complex_d  
    #-----------------------------------------------------------------------------
    #                              4.3 Radiated Wave
    #-----------------------------------------------------------------------------

    #for iw in range(Nw[0]):#use in case one single NEMOH run is used for more than one frequency
    for iB in range(NEM_BODY['nbody']):
        nrspaces = 5 - len(str(iB+2))#be careful this will only work with one frequency
        spaces = " "*nrspaces
        with open(dirname+'/freesurface.{0:s}{1:d}.dat'.format(spaces,iB+2),'r') as fid:#be careful this will only work with one frequency
            allEta = fid.readlines()
    
        iL = 0
        for iX in range(NEM_GRID['Nx']):
            for iY in range(NEM_GRID['Ny']):
                line = [float(a) for a in allEta[2+iL].split()]
                modE[iY,iX] = line[2]
                phaE[iY,iX] = line[3]
                iL += 2
        mR[iB,:,:] = modE
        pR[iB,:,:] = phaE 
        real_mw[iB,:,:] = mR[iB,:,:]*np.cos(pR[iB,:,:])#real value of the wave_field
        complex_mw[iB,:,:] = mR[iB,:,:]*np.sin(pR[iB,:,:])#complex value of the wave_field
        rad_wave_aux[iB,:,:] = (real_mw[iB,:,:] + 1j*complex_mw[iB,:,:])*rad_vel[iB]    
    rad_wave[:,:]=np.sum(rad_wave_aux,axis=0) 
    # Plotting #be careful with the round up of the grid size (I would impose dx and dy if known at the begining)

    #-----------------------------------------------------------------------------
    #                              4.4. Total Wave + Perturbed Wave
    #-----------------------------------------------------------------------------
    tot_wave[:,:]=rad_wave[:,:]+inc_wave[:,:]+diff_wave[:,:]
    pert_wave[:,:]=rad_wave[:,:]+diff_wave[:,:]  
    
    #-----------------------------------------------------------------------------
    #                              5.1. Interpolate NEMOH data
    #-----------------------------------------------------------------------------
    #TODO FIX INTERPOLATION SCRIPT  NF   = 1 because we are runnign for one freq at a time only
    if NEM_ini['interpolation']:
        (rad_waveN, diff_waveN, pert_waveN, inc_waveN, tot_waveN, xvectN, yvectN, NxN,
         NyN) = ne.nemoh_interp(NEM_ini,NEM_GRID,rad_wave,diff_wave,pert_wave,inc_wave,tot_wave,X,Y)
        #TODO check if when passing a rad_waveN x nfrequencies structure the interpolation is done for all the elements        
        #wave fields
        
        return(inc_waveN,diff_waveN,pert_waveN,rad_waveN,tot_waveN, xvectN, yvectN, NxN,NyN, inc_wave,diff_wave,pert_wave,rad_wave,tot_wave,X,Y)
        
    return(inc_wave,diff_wave,pert_wave,rad_wave,tot_wave,X,Y)
  

def wavlen(H,T,dpth):    
    g = 9.81
    dpi = 2*np.pi
    
    L0 = g*T**2/dpi
    L1 = g*T**2/dpi*np.tanh(dpi*dpth/L0)
    
    while (np.abs(L1-L0) > 0.001):
        L0 = L1
        L1 = g*T**2/dpi*np.tanh(dpi*dpth/L0)
    return L1