#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 17 11:41:02 2022

@author: varo

V00 --> includes four photometric models
    lineal (vmag = H + alpha*beta)
    HG bowell (Bowell et al. 1989, in Asteroids II, pp. 524–556) 
        https://ui.adsabs.harvard.edu/abs/1989aste.conf..524B/abstract
    Shevshenko's (Shevshenko 1996, Lunar and Planetary Science 27, 1193)
        https://ui.adsabs.harvard.edu/abs/1996LPI....27.1193S/abstract
    HG1G2 & HG12* (Penttila et al. 2016, Planet. Space Sci. 123, 117)
        https://ui.adsabs.harvard.edu/abs/2010Icar..209..542M/abstract
        https://ui.adsabs.harvard.edu/abs/2016P%26SS..123..117P/abstract


"""
from sys import path
computer = 'mac'
#computer = 'bpz'

if computer == 'mac': 
    path.append('/Users/varo/Documents/Work/Codigos/Phase_curves/')
if computer == 'bpz': 
    path.append('/home/users/dss/alvaro/Codigos/Phase_curves/')

#import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy import optimize
from scipy import interpolate

if   computer == 'mac':
    aux = np.loadtxt('/Users/varo/Documents/Work/Phase_curves_catalogs/phi_func.dat')
elif computer == 'bpz':
    aux = np.loadtxt('/data/alvaro/Phases/phi_func.dat')

phi1 = interpolate.interp1d(aux[:,0], aux[:,1],kind='cubic')
phi2 = interpolate.interp1d(aux[:,0], aux[:,2],kind='cubic')
phi3 = interpolate.interp1d(aux[:,0], aux[:,3],kind='cubic')

def H_hg1g2(fpar,pha,vmag):

    h = fpar[0]
    g1 = fpar[1]
    g2 = fpar[2]
    
    model = h - 2.5*np.log10(g1*phi1(pha)+g2*phi2(pha)+(1.-g1-g2)*phi3(pha)) 
    
    return (model-vmag)**2.0

def H_hg12(fpar,pha,vmag):
    h  = fpar[0]
    g12 = fpar[1]
    
    model = h - 2.5*np.log10(0.84293649*g12*phi1(pha)+0.53513350*(1.-g12)*phi2(pha)+(1.0-0.84293649*g12-0.53513350*(1.-g12))*phi3(pha))

    return (model-vmag)**2.0

    #return abs(model-vmag)

def H_Shev(fpar,pha,vmag):
    a = fpar[0]
    b = fpar[1]
    h = fpar[2]
    
    model = h - a/(1.+pha)+b*pha
    
    return (model-vmag)**2.0

def HG_bowell(fpar,pha,vmag,flag):

# parameters following Bowell et al. 1989 @ asteroids II
    g=fpar[0]
    h=fpar[1]
    
    A1 = 3.332
    B1 = 0.631
    C1 = 0.986
    A2 = 1.862
    B2 = 1.218
    C2 = 0.238

    Phi1L = np.exp(-A1*(np.tan(np.pi*pha/360.0)**B1))
    Phi2L = np.exp(-A2*(np.tan(np.pi*pha/360.0)**B2))

    Phi1S = 1.-(C1*np.sin(np.pi*pha/180.)/(0.119+1.341*np.sin(np.pi*pha/180.)-0.754*np.sin(np.pi*pha/180.)**2.))
    Phi2S = 1.-(C2*np.sin(np.pi*pha/180.)/(0.119+1.341*np.sin(np.pi*pha/180.)-0.754*np.sin(np.pi*pha/180.)**2.))

    W = np.exp(-90.56*np.tan(np.pi*pha/360.0)**2.)

    Phi1 = W*Phi1S + (1.-W)*Phi1L
    Phi2 = W*Phi2S + (1.-W)*Phi2L

    mag = 10.**(-0.4*vmag)

    model = 10.**(-0.4*h) * ((1.-g)*Phi1 + g*Phi2)
  
    if flag == 1:
        return (model - mag)**2.0
    else:
        return np.real(-2.5*np.log10(model))


def main(x,y,model):
    # first guess
    aux = np.polyfit(x,y,1)
    bb = []

    if model == 'hg' :
        fpar = np.array([0.15,aux[1]])
        
        bb = scipy.optimize.leastsq(HG_bowell,fpar,args=(x,y,1))
        bb = bb[0]
        
        return bb[1],bb[0]
    
    if model == 'shev':
        fpar = np.array([0.5,aux[0],aux[1]])
        
        bb = scipy.optimize.leastsq(H_Shev,fpar,args=(x,y))
        bb = bb[0]
        
        if np.greater_equal(bb[0],0.0):
            pass
            
        else:
            bb=[-9.0,-9.0,-9.0]
            
        return bb[-1]-bb[0],bb[0],bb[1]
    
    if  model == 'lineal':
        bb = aux
        return bb[1],bb[0]
    
    if model == 'hg12':
        cte = (1.0-0.53513350)/(0.84293649-0.53513350)

        fpar = np.array([aux[1],0.6])
        bb = scipy.optimize.leastsq(H_hg12,fpar,args = (x,y))
        bb = bb[0]

        if np.greater_equal(bb[1],0.0)*np.less_equal(bb[1],cte):
            pass
        
        else:
            bb = [-9.,-9.]
            
        return bb
    
    if model == 'hg1g2':
        fpar = np.array([aux[1],0.4,0.4])
        bb = scipy.optimize.leastsq(H_hg1g2,fpar,args = (x,y))
        bb = bb[0]
        
        if np.greater_equal(bb[1],0.0)* np.greater_equal(bb[2],0.0)*np.greater_equal((1.0-bb[1]-bb[2]),0.0):
            pass
        
        else:
            bb = [-9.,-9.,-9.]

        return bb
