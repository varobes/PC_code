'''
@author: varo

version 00 --> creates the PDF of dmag using a first estimate of H from HG12 and the errors in the mags.
           --> this works, for now, with SLOAN filter set
version 01 --> takes into account (a random) rotational phase for the observation.
version 02 --> passing to python3
           --> giving option include effect of rotational phase angle cos(phi)
           --> using V mag or proxy
version 03 --> adapted to work with NIR data
           --> using solar color in the visible to estimate V mag from Y or J or H or K
'''

from sys import path
computer = 'mac'
#computer = 'bpz'

if computer == 'mac': 
    path.append('/Users/varo/Documents/Work/Codigos/Phase_curves/')
if computer == 'bpz': 
    path.append('/home/users/dss/alvaro/Codigos/Phase_curves/')

import scipy
import numpy as np
from numpy import random
import matplotlib.pyplot as plt
import importlib

import photmod00 as PM
importlib.reload(PM)

PMversion = '00'
PDMversion = '03'

#plt.rcParams.update({'axes.labelsize':25})
#plt.rcParams.update({'font.size':15})#18})
#plt.rcParams.update({'xtick.labelsize':22})
#plt.rcParams.update({'ytick.labelsize':22})

def main(id,xt,yt,data,np,method,model):
    # id     == object identification
    # xt     == phase angle / shape(xt) = np
    # yt     == (V mag & error or proxy) / shape 2xnp
    # data   == mags from list
    # np     == minimum number of points to consider in phase curves
    # method == use or not aleatory rotational phase
    #   0 = no (default 0)
    #   1 = yes
    
    vflag = 0
    flag  = 0
    iter = 10000
    check = np.arange(iter/5,iter+1,iter/5)
    
    if method == 0 :
        binmag = np.arange(0.0,3.01,0.02)
    elif method == 1 :
        binmag = np.arange(-1.5,1.51,0.01)
    
    print (' ')
    print ('   ...Estimating PDF delta_mag...')
    
    mmmm = []
    if method == 0 :
        print ('      No estimation of rotational phase')
        mmmm = 'without sin(phi)'
        
    elif method == 1 :
        print ('      Estimating rotational phase')
        mmmm = 'using sin(phi)'

    print ('   ', id)
    
    if len(np.shape(yt)) == 2:
        y = yt[:,0]
        s = yt[:,1]
    else:
        y = yt.copy()
        s = np.ones(len(y))*0.05
    
    ok = np.greater(y,-9.0)
    x = xt[ok]
    y = y[ok]
    s = s[ok]
    lent = len(x)
    vmed = np.median(y)
    vsig = max(s)
    out = np.c_[vmed,vsig]


    if lent > np-1: # how many points am I to consider to create the phase curves         
        #auxiliary arrays
        vacg=np.empty(iter)    

        ll = 0
        # in this step we'll just assume errors as gaussian
        for l in range(iter):
            yn = y + random.randn(lent)*s
            aux = PM.main(x,yn,model)
            #print(aux)
            vacg[ll] = aux[0]
            aux = vacg[ll]
            
            if l in check: print('      valid solutions =',ll, '   number of iterations = ',l)
            
            if np.isnan(aux) or aux < -4.:
                pass
            
            else:
                ll = ll + 1
        
        print('      valid solutions =',ll, '   after iterations = ',l+1)

        vacg = vacg[:ll]

        if ll == 0 : # no good result
            vflag = 1
      
    else: # number of points less than requested
        flag = 1
    
    
    if vflag ==1 or flag == 1 :
        pass
    
    else:
        vmag = vacg.copy()
        aux = np.histogram(vmag,bins=30)
      
        y = aux[0]*1.0
        y = y/y.sum()
        aux1 = aux[1]
    
        if method == 0:
            dm = np.zeros((150))
        
        elif method == 1:
            dm =np.zeros((300))
            
        #dm1 = np.zeros_like(dm)
        bin = binmag[1]-binmag[0]
        x = binmag[:-1]+bin/2.
    
        if method == 0:
            x = x/2.
            x = np.r_[-x[::-1],x]
            
            for jj in range(30):
                ok = np.greater_equal(data[:,0],aux1[jj])*np.less(data[:,0],aux1[jj+1])
                temp = data[ok,2]
                aux = np.histogram(temp,binmag)
                dm = aux[0]*y[jj] + dm
                
            dm = np.r_[dm[::-1],dm]

        elif method == 1:
            for jj in range(30):
                ok = np.greater_equal(data[:,0],aux1[jj])*np.less(data[:,0],aux1[jj+1])
                temp = data[ok,2]
                sin_theta=random.rand(len(temp))*2.0*np.pi
                sin_theta=np.sin(sin_theta)
                temp = temp*sin_theta
                aux = np.histogram(temp,binmag)
                dm = aux[0]*y[jj] + dm
                #dm1 = dm.copy()
    
        dm = dm / dm.sum()

        out = np.r_[out,scipy.c_[x,dm]]
        
    np.savetxt('dmags/'+id+'_dmag-'+PDMversion+'.dat',out,header=mmmm+' // PMversion='+PMversion+' // PMmethod = '+model)  
    #np.savetxt('test_pdm.dat',out,header=mmmm+' // PMversion='+PMversion+' // PMmethod = '+model)  

    return out