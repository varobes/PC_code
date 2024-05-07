"""
Created on Thu Mar 10 13:41:24 CET 2022

@author: varo
Conversion between photometric systems
V00 --> ugriz to VR
V01 --> NIR data to V, using solar colors
    #### Solar colors
    uses NIR solar colors from Casagrande et al. (2012)
        https://ui.adsabs.harvard.edu/abs/2012ApJ...761...16C/abstract
    uses transformation from 2MASS to VISTAS from Popescu et al. (2016) 
        https://ui.adsabs.harvard.edu/abs/2016A%26A...591A.115P/abstract

"""

from sys import path
computer = 'mac'
#computer = 'bpz'

if computer == 'mac': 
    path.append('/Users/varo/Documents/Work/Codigos/Phase_curves/')
if computer == 'bpz': 
    path.append('/home/users/dss/alvaro/Codigos/Phase_curves/')

import numpy as np

def main(inp,from_psys,to_psys):
    
    # y-j, y-k, j-h, j-k, h-k
    vista_c = [0.196,0.532,0.255,0.336,0.082] #vista colors
    vj = 0.220022 #Johnson's V - vista's J
    
    
    if from_psys == 'ugriz' and to_psys == 'VR':
        '''
        inp = shape Nx6
        '''
        
        mag1 = inp[:,0]-0.59*(inp[:,0]-inp[:,2])-0.01
        sig1 = (((1.-0.59)*inp[:,1])**2.0+(0.59*inp[:,3])**2.0)**0.5
    
        mag2 = mag1-(1.09*(inp[:,2]-inp[:,4]+0.22))
        sig2 = ( (1.09*inp[:,3])**2.0 + (1.09*inp[:,5])**2.0 + sig1**2.0 )**0.5
        
        return np.c_[mag1,sig1,mag2,sig2]
        
    #### All error budget in the transformations below comes from 2MASS. Popescu didn't compute error in their transformation.
    if from_psys =='yjkh' and to_psys =='V':
        
        j = 0 # first uses Y
        cte = vj - vista_c[0] # (v-j) - (y-j)
        mag1 = inp[:,j] + cte
        sig1 = (inp[:,j+1]**2.0+0.01**2)**0.5

        out = np.where(mag1<0)
        out = out[0]
        
        if len(out) > 0:
            j = 2 # Second uses J
            aux = inp[out,j]
            cte = vj  # (v-j)
            mag1[out] = inp[out,j] + cte
            sig1[out] = (inp[out,j+1]**2.0+0.01**2)**0.5

            out = np.where(mag1<0)
            out = out[0]

            if len(out) > 0:
                j = 4 # Third uses H
                aux = inp[out,j]
                cte = vj - vista_c[2] # (v-j) - (y-j)
                mag1[out] = inp[out,j] + cte
                sig1[out] = (inp[out,j+1]**2.0+0.01**2)**0.5

                out = np.where(mag1<0)
                out = out[0]

                if len(out) > 0:
                    j = 4 # Fourth uses H
                    aux = inp[out,j]
                    cte = vj - vista_c[3] # (v-j) - (y-j)
                    mag1[out] = inp[out,j] + cte
                    sig1[out] = (inp[out,j+1]**2.0+0.01**2)**0.5
    
                    out = np.where(mag1<0)
                    out = out[0]

                    if len(out) > 0:
                        mag1[out] = -10
                        sig1[out] = -10
        
        return np.c_[mag1,sig1]

        
    
