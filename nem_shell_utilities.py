# -*- coding: utf-8 -*-
"""
Created on Thu Oct 18 09:54:09 2018
@author: gveraofe
mod. PhilDog 18.02.19
rename nemoh_shell_utilities
this is used to run NEMOH from both the individual NEMOH thread asand from the MILDwave coupling thread which requyires NEMOH 
"""
import numpy as np
import os
import sys
import spectral_functions as sf
import math as math
#put the relative path to your NM_MW_shared folder
#path_shell = 'C:/Users/gveraofe/Documents/A_Ch_4_Simulations/Python'
#FUNCIONS FOR RUNNING NEMOH IMPORT FROM MW NEMOH SHARED
#sys.path.append(path_shell)#add the path where the files containing MW and NEM functions are located
sys.path.append(os.path.abspath(os.path.join(os.path.dirname("..//MW_NEM_shared"),".." )) )
from MW_NEM_shared import nem_post_utilities as pu
from MW_NEM_shared import nem_utilities as ne
import random

#Creates all the input variables needed for running NEMOH once the type of case and input wave conditions are defined
def create_input (NEM_ini):
    regCheck = NEM_ini['regular']
    irrCheck = NEM_ini['irregular']
    dirCheck = NEM_ini['directional']

    if regCheck:
        NEM_ini['f']=[]
        NEM_ini['w']=[]
        NEM_ini['amp']=[]
        for ii in range(len(NEM_ini['T'])):
            NEM_ini['f'].append(1/NEM_ini['T'][ii])
            print(1/NEM_ini['T'][ii])
        for ii in range(len(NEM_ini['f'])):
            NEM_ini['w'].append(2*np.pi*NEM_ini['f'][ii])
        for ii in range(len(NEM_ini['H'])):
            NEM_ini['amp'].append(NEM_ini['H'][ii]/2)

    elif irrCheck:
        NEM_ini['fp'] = [1/Tp for Tp in NEM_ini['Tp']]# Peak frequency in Hz
        NEM_ini['wp'] = [2*np.pi*fp for fp in NEM_ini['fp']] 
        NEM_ini['f']=[[] for jj in range(len(NEM_ini['Tp']))]
        NEM_ini['w']=[[] for jj in range(len(NEM_ini['Tp']))]
        NEM_ini['T']=[[] for jj in range(len(NEM_ini['Tp']))]
        for jj in range(len(NEM_ini['Tp'])):
            # GAEL are we sure we use this            
            for fp in NEM_ini['fp']:            
                NEM_ini['f'][jj]=(list(reversed(np.round(np.linspace(NEM_ini['fini'][jj],NEM_ini['fend'][jj],NEM_ini['Nf'],endpoint = True),3)))) #Linespacing the frequency maybe is better to do it with the period??                               #Limits should be adjusted depending on the peak period
                NEM_ini['T'][jj]=([round((1/f),3) for f in NEM_ini['f'][jj]  ])
                NEM_ini['w'][jj]=([2*np.pi*f for f in NEM_ini['f'][jj]  ])
            
        if NEM_ini['spectra'] == 'JS':
                # FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
            sf.JS(NEM_ini)
                # FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
        elif NEM_ini['spectra'] =='PM':
                # FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
            sf.PM(NEM_ini)
               # FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF            
       
        if dirCheck:
            #filename = 'dir_nemoh_input.txt'#Directiosn used in the MILDwave SC simulation
            #for ii in NEM_ini['dirfile']:
                #deg_MW = np.loadtxt(os.path.join(NEM_ini['dirfile'][ii],filename),skiprows=0)    
                #deg_MW = deg_MW[:,1]
                #NEM_ini['deg'] =list(reversed(deg_MW))   #TODO IMPLEMENT FOR DIFFERENT DIRECTIONAL WAVES CASES
           deg_aux = sf.DSM(NEM_ini)
        #CORRECTION OF THE ANGLES ACORDING TO MILDWAVE CORRECTION IN THE WAVE GENERATION LINE IN RELATION WITH GIRD SIZE AND WAVE GENERATION LINE SIZE  
           for ii in range(len(NEM_ini['deg_main'])):
               (L_T) = [wavlen(H,T,NEM_ini['depth'][ii]) for H,T in zip(NEM_ini['H'][ii],NEM_ini['T'][ii]) ]#L of each T
               k_T = [2*np.pi/L for L in L_T]#angular wave number
               #dir_name_Tp = 'T_'+'{:06.3f}'.format(NEM_ini['Tp'][ii])
               #tree = et.parse(os.path.join(MW_ini['dir_mw_irr'],dir_name_Tp,"MILDwave.xml"))
               jm1 = 777# int(tree.find(".//Ny").text)#TODO HOW TO RUN JUST SHORT CRESTED IN NEMOH
               #dx = 0.08#float(tree.find(".//dx").text )
               dy = 2.0#float(tree.find(".//dy").text)q
               sinRw = [math.sin(Theta) for Theta in deg_aux[ii]]
               rval = [np.round(k*sin*(jm1 - 1)*dy/(2*np.pi))*2*np.pi/((jm1-1)*dy*k) for sin,k in zip(sinRw,k_T)]               
               deg = np.zeros(NEM_ini['Nf'])
               #TODO INPUT ERROR MESSAGE IF THE PERIODIC LENGHT  IS NOT ENOUGH THE PROGRAM GIVES THE MAIN DIRECTION
               #     PERIODIC LENGHT HAS TO BE 3L AT LEAST
               for jj in range(len(rval)):
                  if (rval[jj] > 1.0):
                      deg[jj] = math.asin(rval[jj] - 2*np.pi/((jm1-1)*dy*k_T[jj]))*180/np.pi + NEM_ini['deg_main'][ii] #INPUT in MW AND NEMOH IS DEGREES and WE HAVE TO SHIFT FOR THE MAIN DIRECTION
                  elif (rval[jj] < -1.0):
                      deg[jj] = math.asin(rval[jj]+2*np.pi/((jm1-1)*dy*k_T[jj]))*180/np.pi + NEM_ini['deg_main'][ii]   #INPUT in MW AND NEMOH IS DEGREES and WE HAVE TO SHIFT FOR THE MAIN DIRECTION 
                  else:
                      deg[jj] = math.asin(rval[jj])*180/np.pi + NEM_ini['deg_main'][ii]     #INPUT in MW AND NEMOH IS DEGREES and WE HAVE TO SHIFT FOR THE MAIN DIRECTION
               NEM_ini['deg'].append(deg)
        else:
            NEM_ini['deg'] = [deg*NEM_ini['Nf'] for deg in NEM_ini['deg']]
        
    return(NEM_ini)
#GAEL was ist nrFreq? for possible running og NEMOH for manty frequencies suimultaneaously

def runNEM(nrFreq,NEM_ini,NEM_BODY,NEM_advOps,NEM_GRID,NEMOHDir):
    # FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
    NEM_OUT = pu.ini(NEM_ini) #TODO check in post_utliities
    # FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
    #!! zip canonly iterate over shortest length list #iterate over irr ~f~ separtely? TODO -> remove iteration over _f_ for now
    for ii in range(len(NEM_ini['case_dir'])):
    
        if NEM_ini['regular']:
            Fe_aux = []
            Ma_aux = []
            Bhyd_aux = []
            Rao_aux = []
            rad_vel_aux=[]
                
            if NEM_advOps['fsCheck']:
                #onlty creaete the emprty wave field cvectors if freee surface condition = true
                inc_wave_aux = [[] for jj in range(len(NEM_ini['T']))]
                diff_wave_aux = [[] for jj in range(len(NEM_ini['T']))]
                rad_wave_aux = [[] for jj in range(len(NEM_ini['T']))]
                pert_wave_aux = [[] for jj in range(len(NEM_ini['T']))]
                tot_wave_aux = [[] for jj in range(len(NEM_ini['T']))]
                
                if NEM_ini['interpolation']:
                    # if itnerpolation = true also create empty liosts for interp values
                    inc_wave_auxN = [[] for jj in range(len(NEM_ini['T']))]
                    diff_wave_auxN = [[] for jj in range(len(NEM_ini['T']))]
                    rad_wave_auxN = [[] for jj in range(len(NEM_ini['T']))]
                    pert_wave_auxN = [[] for jj in range(len(NEM_ini['T']))]
                    tot_wave_auxN = [[] for jj in range(len(NEM_ini['T']))]

            for jj in range(len(NEM_ini['T'])) : # iteration over T
                    
                Tjj = NEM_ini['T'][jj]
                deg, depth = NEM_ini['deg'][0],NEM_ini['depth'][0]#TODO various degrees
        #            for deg in NEM_ini['deg']:
                fmStart = 1./Tjj
                fmEnd = 1./Tjj
                NEM_advOps ['dirStart'] = deg
                NEM_advOps ['dirStop'] =  deg
                omega = [nrFreq,2.0*np.pi*fmStart,2.0*np.pi*fmEnd]#TODO RUN ALL FREQUENCIES IN A SINGLE NEMOH RUN
                
                ne.writeCalFile(NEM_advOps['rhoW'],depth,omega,NEM_BODY['cG'],NEM_BODY['dof'],aO=NEM_advOps,nbody=NEM_BODY['nbody'],
                                xG=NEM_BODY['xBody'],yG=NEM_BODY['yBody'])# Write nemoh.cal file
                    # FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
                ne.runNemoh(show_console = NEM_advOps['Show_Console'],nbody=NEM_BODY['nbody'])# RUN NEMOH FORTRAN
                    # FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
                #TODO IMPLEMET FOR Nw > 1 it should work as list comprehesion but problaby need to modify ne.readNemohResults
                #TODO check the outputs for the nemoh dir ii's
                (w,beta,Ma,Bhyd,Fe) = ne.readNemohResults(NEM_BODY['dof'],os.path.join(NEMOHDir,'Calculation','results'))

                (Rao) = pu.RAO(w,NEM_ini,NEM_BODY,Fe,Ma,Bhyd)
                (rad_vel) = pu.rad_vel(w,Rao)
                #append the non-freee surface paramters 
                Fe_aux.append(Fe)
                Ma_aux.append(Ma)
                Bhyd_aux.append(Bhyd)
                Rao_aux.append(Rao)
                rad_vel_aux.append(rad_vel)

                #TODO ideallywe need to put this in a prepro and call from main nem_shell.py -> too many subcalls to subcalls of funcitons
                     #onlty creaete the emprty wave field cvectors if freee surface condition = true
                if NEM_advOps['fsCheck']:
                    if NEM_ini['interpolation']:
                            # if itnerpolation = true also create empty liosts for interp values
                            #FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
                            (inc_waveN,diff_waveN,pert_waveN,rad_waveN,tot_waveN, xvectN, yvectN, NxN,NyN,
                             inc_wave,diff_wave,pert_wave,rad_wave,tot_wave,xvect,yvect) = pu.wave_field(NEM_ini,NEM_GRID,NEM_BODY,rad_vel,NEMOHDir,1/Tjj,deg,depth)
                            #FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
                            inc_wave_auxN[jj]=(inc_waveN)
                            diff_wave_auxN[jj]=(diff_waveN)
                            rad_wave_auxN[jj]=(rad_waveN)
                            pert_wave_auxN[jj]=(pert_waveN)
                            tot_wave_auxN[jj]=(tot_wave)
                    else :# no interpolation = defualt
                        #FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
                            (inc_wave,diff_wave,pert_wave,rad_wave,tot_wave,xvect,yvect
                             ) = pu.wave_field(NEM_ini,NEM_GRID,NEM_BODY,rad_vel,NEMOHDir,1/Tjj,deg,depth)
                        #FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
                    inc_wave_aux[jj]=(inc_wave)
                    diff_wave_aux[jj]=(diff_wave)
                    rad_wave_aux[jj]=(rad_wave)
                    pert_wave_aux[jj]=(pert_wave)
                    tot_wave_aux[jj]=(tot_wave)

                
                
        elif NEM_ini['irregular']: #irregular
            #!! zip canonly iterate over shortest length list #iterate over irr ~f~ separtely? TODO -> remove iteration over _f_ for now
            Fe_aux =  [[] for jj in range(len(NEM_ini['Tp']))]
            Ma_aux =  [[] for jj in range(len(NEM_ini['Tp']))]
            Bhyd_aux =  [[] for jj in range(len(NEM_ini['Tp']))]
            Rao_aux =  [[] for jj in range(len(NEM_ini['Tp']))]
            rad_vel_aux= [[] for jj in range(len(NEM_ini['Tp']))]
            if NEM_advOps['fsCheck']:
                    #onlty creaete the emprty wave field cvectors if freee surface condition = true
                    inc_wave_aux = [[] for jj in range(len(NEM_ini['Tp']))]
                    diff_wave_aux = [[] for jj in range(len(NEM_ini['Tp']))]
                    rad_wave_aux = [[] for jj in range(len(NEM_ini['Tp']))]
                    pert_wave_aux = [[] for jj in range(len(NEM_ini['Tp']))]
                    tot_wave_aux = [[] for jj in range(len(NEM_ini['Tp']))]
                    
                    if NEM_ini['interpolation']:
                        # if itnerpolation = true also create empty liosts for interp values
                        inc_wave_auxN = [[] for jj in range(len(NEM_ini['Tp']))]
                        diff_wave_auxN = [[] for jj in range(len(NEM_ini['Tp']))]
                        rad_wave_auxN = [[] for jj in range(len(NEM_ini['Tp']))]
                        pert_wave_auxN = [[] for jj in range(len(NEM_ini['Tp']))]
                        tot_wave_auxN = [[] for jj in range(len(NEM_ini['Tp']))]
    
            for jj in range(len(NEM_ini['Tp'])) : # iteration over Tp
                Tjj = NEM_ini['T'][jj]
                deg = NEM_ini['deg'][jj]
                depth = NEM_ini['depth'][0] #TODO various degrees
    
                for nn in range(len(Tjj)): # loop over subperiods
        #          for ii in range(NEM_ini['Nf']):
                    fmStart = 1./Tjj[nn]
                    fmEnd = 1./Tjj[nn]
                    NEM_advOps ['dirStart'] = deg[nn]
                    NEM_advOps ['dirStop'] =  deg[nn]
                    omega = [nrFreq,2.0*np.pi*fmStart,2.0*np.pi*fmEnd]#TODO RUN ALL FREQUENCIES IN A SINGLE NEMOH RUN
                    ne.writeCalFile(NEM_advOps['rhoW'],depth,omega,NEM_BODY['cG'],NEM_BODY['dof'],aO=NEM_advOps,nbody=NEM_BODY['nbody'],
                                    xG=NEM_BODY['xBody'],yG=NEM_BODY['yBody'])# Write nemoh.cal file
                    ne.runNemoh(show_console = NEM_advOps['Show_Console'],nbody=NEM_BODY['nbody'])# RUN NEMOH FORTRAN
                    #TODO IMPLEMET FOR Nw > 1 it should work as list comprehesion but problaby need to modify ne.readNemohResults
                    #TODO check the outputs for the nemoh dir ii's
                    (w,beta,Ma,Bhyd,Fe) = ne.readNemohResults(NEM_BODY['dof'],os.path.join(NEMOHDir,'Calculation','results'))
                    (Rao) = pu.RAO(w,NEM_ini,NEM_BODY,Fe,Ma,Bhyd)
                    (rad_vel) = pu.rad_vel(w,Rao)
                    #append the non-freee surface paramters 
                    Fe_aux[jj].append(Fe)
                    Ma_aux[jj].append(Ma)
                    Bhyd_aux[jj].append(Bhyd)
                    Rao_aux[jj].append(Rao)
                    rad_vel_aux[jj].append(rad_vel)
                        #TODO ideallywe need to put this in a prepro and call from main nem_shell.py -> too many subcalls to subcalls of funcitons
                    if NEM_advOps['fsCheck']:
                        if NEM_ini['interpolation']:
                            # if itnerpolation = true also create empty liosts for interp values
                            #FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
                            (inc_waveN,diff_waveN,pert_waveN,rad_waveN,tot_waveN, xvectN, yvectN, NxN,NyN,
                             inc_wave,diff_wave,pert_wave,rad_wave,tot_wave,xvect,yvect) = pu.wave_field(NEM_ini,NEM_GRID,NEM_BODY,rad_vel,NEMOHDir,1/Tjj[nn],deg[nn],depth)
                            #FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
                            inc_wave_auxN[jj].append(inc_waveN)
                            diff_wave_auxN[jj].append(diff_waveN)
                            rad_wave_auxN[jj].append(rad_waveN)
                            pert_wave_auxN[jj].append(pert_waveN)
                            tot_wave_auxN[jj].append(tot_wave)
                        else :# no interpolation = defualt
                        #FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
                            (inc_wave,diff_wave,pert_wave,rad_wave,tot_wave,xvect,yvect
                             ) = pu.wave_field(NEM_ini,NEM_GRID,NEM_BODY,rad_vel,NEMOHDir,1/Tjj[nn],deg[nn],depth)
                        #FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
                    inc_wave_aux[jj].append(inc_wave)
                    diff_wave_aux[jj].append(diff_wave)
                    rad_wave_aux[jj].append(rad_wave)
                    pert_wave_aux[jj].append(pert_wave)
                    tot_wave_aux[jj].append(tot_wave)

        # append the list for each peak period to the NEM_Out dict
    NEM_OUT['Fe']=(Fe_aux)
    NEM_OUT['Ma']=(Ma_aux)
    NEM_OUT['Bhyd']=(Bhyd_aux)
    NEM_OUT['Rao']=(Rao_aux)
    NEM_OUT['rad_vel']=(rad_vel_aux)
    
    if NEM_advOps['fsCheck']:
        NEM_OUT['inc_wave']=(inc_wave_aux)
        NEM_OUT['diff_wave']=(diff_wave_aux)
        NEM_OUT['rad_wave']=(rad_wave_aux)
        NEM_OUT['pert_wave']=(pert_wave_aux)
        NEM_OUT['tot_wave']=(tot_wave_aux)
        NEM_OUT['xvect']= xvect
        NEM_OUT['yvect']= yvect
        #pn.calcABKFe        
        if NEM_ini['interpolation']:   
            NEM_OUT['NxN']= NxN
            NEM_OUT['NyN']= NyN  
            NEM_OUT['xvectN']= xvectN
            NEM_OUT['yvectN']= yvectN  
            NEM_OUT['xvect']= xvect
            NEM_OUT['yvect']= yvect
            NEM_OUT['inc_waveN']=(inc_wave_auxN)
            NEM_OUT['diff_waveN']=(diff_wave_auxN)
            NEM_OUT['rad_waveN']=(rad_wave_auxN)
            NEM_OUT['pert_waveN']=(pert_wave_auxN)
            NEM_OUT['tot_waveN']=(tot_wave_auxN)    
    return(NEM_OUT)
        
#TODO add nm_dior and Tp ii iteration
def irr_wave(NEM_ini,NEM_GRID,NEM_OUT,time,dt):
    time = np.linspace(0,time,time/dt)
    kd_time = []
    #TODO KILL ZIP
    for T,a,w,Hs in zip(NEM_ini['T'],NEM_ini['amp'],NEM_ini['w'],NEM_ini['Hs_m01']):
        Nvar_aux = np.zeros((NEM_GRID['Ny'],NEM_GRID['Nx']))
        r_phase = np.zeros(len(T))
        for ii in range(len(T)):
            r_phase[ii] = random.uniform(0,2)*np.pi
        for tt in range(np.shape(time)[0]):
            nm_aux = np.zeros((NEM_GRID['Ny'],NEM_GRID['Nx']))
            for ii in range(len(T)):
                nm_tw_P = a[ii]*NEM_OUT['tot_wave'][ii] 
                nm_tw_P = np.real(np.abs(nm_tw_P) * np.exp(-1j*(w[ii]*time[tt]-r_phase[ii]+np.angle(nm_tw_P))))
                nm_aux += (nm_tw_P)
            Nvar_aux += ((nm_aux)**2)
    
        Nvar = Nvar_aux /time *dt
        kd = 4*np.sqrt(Nvar) / Hs
        kd[kd>=3]=np.nan
        kd_time.append(kd)
    return(kd_time)
                            
def wavlen(H,T,d):    
    g = 9.81
    dpi = 2*np.pi
    
    L0 = g*T**2/dpi
    L1 = g*T**2/dpi*np.tanh(dpi*d/L0)
    
    while (np.abs(L1-L0) > 0.001):
        L0 = L1
        L1 = g*T**2/dpi*np.tanh(dpi*d/L0)
        return L1