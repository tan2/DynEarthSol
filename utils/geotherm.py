#!/usr/bin/env python
import numpy as np
from math import sqrt
from scipy.special import erf
import sympy as sp

import matplotlib.pyplot as plt

def half_space_cooling_T(z, T0, Tm,  age_in_myrs, alpha):
    myrs2sec = 86400 * 365.2425e6

    T = T0 + (Tm - T0) * erf(z /
            sqrt(4 * alpha * age_in_myrs * myrs2sec) )
    return T

def continental_radiogenic_T(z,T0,hr,k,qm,rhoH0):
    T = T0 + qm/k*z + rhoH0*hr**2/k*(1-np.exp(-z/hr))
    return T

def continental_radiogenic_T2(z,T0,bdy_arr,hr,k_arr,qm,rhoH0_arr):
    '''
    No heat production below the moho
    '''
    dT_h_layer = dT_h(bdy_arr[:-1],bdy_arr[1:],k_arr,rhoH0_arr,hr)
    dtlayer = dT_h_layer + (bdy_arr[1:]-bdy_arr[:-1])*qm/k_arr
    
    T = np.zeros_like(z)
    for i, zt in enumerate(z):
        is_value = False
        for j in range(len(k_arr)):
            if zt <= bdy_arr[j+1]:
                dTh = dT_h(bdy_arr[j],zt,k_arr[j],rhoH0_arr[j],hr)
                T[i] = T0 + np.sum(dtlayer[:j]) + qm/k_arr[j]*(zt - bdy_arr[j]) + dTh
                is_value = True
                break
        if not is_value:
            raise ValueError(f'z={zt} is out of range')
    return T

def dT_h(zstart,zend,k,rhoH0,hr):
    diff = np.exp(-zstart/hr) - np.exp(-zend/hr)
    return diff*rhoH0*hr**2/k

def sum_dt_of_layers(bdy_arr,k_arr,rhoH0_arr,hr):
    dt = dT_h(bdy_arr[:-1],bdy_arr[1:],k_arr,rhoH0_arr,hr)
    return np.sum(dt)

def avg_k(k_arr,thick_arr):
    return np.sum(thick_arr)/np.sum(thick_arr/k_arr)

def get_hz_profile(z,hr,rhoH0_arr,bdy_arr):
    hz = np.zeros_like(z)
    for i, zt in enumerate(z):
        is_value = False
        for j in range(len(rhoH0_arr)):
            if zt <= bdy_arr[j+1]:
                hz[i] = rhoH0_arr[j]*np.exp(-zt/hr)
                is_value = True
                break
        if not is_value:
            raise ValueError(f'z={zt} is out of range')
    return hz

def main():
    # material properties
    mat_rho = [3300, 2800, 2800] # kg/m^3 density
    mat_k = [3.3, 2.5, 2.5] # W/mK thermal conductivity
    mat_cp = [1000, 1000, 1000] # J/kgK specific heat capacity
    mat_H0 = [0, 1e-10,6.0e-10] # W/kg heat production
    # 9.6e-10 W/kg for granite
    hr = 33e3 # km length_scale_for_the_decrease_of_heat_production

    mat_k = np.array(mat_k)
    mat_rho = np.array(mat_rho)
    mat_cp = np.array(mat_cp)
    mat_H0 = np.array(mat_H0)

    T0 = 273
    Tm = 1300 + 273

    Zulb = 22e3 # km
    Zmoho = 33e3 # km
    Zbot = 120e3 # km
    
    phase_arr = [2, 1, 0]
    bdy_arr = [0, Zulb, Zmoho, Zbot]
    bdy_arr = np.array(bdy_arr)
    thick_arr = [bdy_arr[i+1]-bdy_arr[i] for i in range(len(bdy_arr)-1)]
    thick_arr = np.array(thick_arr)
    
    k_arr = mat_k[phase_arr]
    rho_arr = mat_rho[phase_arr]
    H0_arr = mat_H0[phase_arr]
    nlayer = len(phase_arr)

    rhoH0_arr = rho_arr * H0_arr
    # Brune 2014
    rhoH0_arr[0] = 1.5e-6 # W/m^3
    rhoH0_arr[1] = 0.2e-6 # W/m^3
    
    z = np.linspace(0, Zbot, int(Zbot*0.1)+1)

    hz3 = get_hz_profile(z,hr,rhoH0_arr,bdy_arr)
    
    # calculate temperature increase by radiogenic heat production in the crust    
    dT = sum_dt_of_layers(bdy_arr,k_arr,rhoH0_arr,hr)
    qm = (Tm-T0-dT)/sum(thick_arr) * avg_k(k_arr,thick_arr)
    
    temp = continental_radiogenic_T2(z,T0,bdy_arr,hr,k_arr,qm,rhoH0_arr)

    print(f'temperature at surface: {temp[0]-273:.0f} degC')
    print(f'temperature at moho: {np.interp(Zmoho,z,temp)-273:.0f} degC')
    print(f'temperature at bottom: {temp[-1]-273:.0f} degC')

    # plot
    fig, ax = plt.subplots(figsize=(4,8))

    ax.vlines(Tm-273,0,Zbot/1e3,ls='--',color='grey',lw=1)
    ax.text(Tm-273+0.1,0.1,"Tm",fontsize=12,color='grey')
    
    for bdy in bdy_arr[1:-1]:
        ax.hlines(bdy*1e-3,0,1500,ls='--',color='grey',lw=1)

    color = 'k'
    ax.plot(temp-273, z/1e3,'-',color=color,label="Radiogenic")
    for bdy in bdy_arr[1:-1]:
        tmp = np.interp(bdy,z,temp)
        ax.plot([tmp-273],[bdy*1e-3],'o',color=color,ms=4,mfc='w')
        ax.text(tmp-273+10,bdy*1e-3-1,f"{tmp-273:.0f}"+'$^\circ{C}$',fontsize=12,ha='left',va='bottom',color=color)

    ax.plot([temp[0]-273],[0.],'o',color=color,ms=6,mfc='w')
    ax.text(100,0.5,f"{temp[0]-273:.0f}"+'$^\circ{C}$',fontsize=12,ha='left',va='top',color=color)
    ax.plot([temp[-1]-273],[Zbot*1e-3],'o',color=color,ms=4,mfc='w')
    ax.text(900,Zbot*1e-3-0.3,f"{temp[-1]-273:.0f}"+'$^\circ{C}$',fontsize=12,ha='left',va='bottom',color=color)

    ax.text(300,93,f'q={qm*1e3:.0f} mW/m$^2$',fontsize=12,ha='left',va='top',color=color)

    ax2 = ax.twiny()
    
    ax2.plot(hz3*1e6,z/1e3,'--',color='r',label="Heat production",lw=1)
    ax2.set_xlabel(r"Heat production ($\mu$W/m3)")
    ax2.legend(loc='upper right')
    ax2.set_xlim(0,4)

    ax.legend(loc='lower left')
    ax.set_xlabel(r"Temperature ($^\circ{C}$)")
    ax.set_ylabel("Depth (km)")
    ax.set_ylim(Zbot/1e3,0)
    ax.set_xlim(0,1500)
    ax.grid(ls='--')
    fig.tight_layout()
    filename = f'geo-{nlayer:d}-{rhoH0_arr[0]:.1e}-{hr/1e3:.0f}km.png'
    # filename = 'test.png'
    
    print(f'save figure to {filename}')
    fig.savefig(filename)
    
    return

if __name__ == "__main__":
    main()