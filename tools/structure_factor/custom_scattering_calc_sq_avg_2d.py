#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  1 11:13:05 2021

@author: qinyu
adapted by pkawak for MCPC/OPWL config file (.xyz) usage
"""

import sys
import freud
import numpy as np
from numba import cuda, njit, float64
from timeit import default_timer as timer
from numba.cuda.random import create_xoroshiro128p_states, xoroshiro128p_uniform_float64
import matplotlib.pyplot as plt
import math

file_name = "config/config_00000.out"
if len(sys.argv) >= 2:
  file_name = sys.argv[1]

# load params
with open(file_name.split(':')[0]) as f:
  lines=f.readlines()[1].split(' ')
  Nc   = int(lines[1])
  Nb   = int(lines[2])
  Npo  = Nc*Nb
  lx   = float(lines[3])
  ly   = lx 
  lz   = lx 

q_range = 4.2
num_data = int(math.ceil(10*q_range))

min_q = 2/10 #2/lx
max_q = q_range*np.pi

max_nums = 200
# initiate the wave vector
sizex = min(max_nums, int(np.ceil(lx))*num_data)
sizey = min(max_nums, int(np.ceil(ly))*num_data)
sizez = min(max_nums, int(np.ceil(lz))*num_data)

if len(sys.argv) > 3:
  q_range = float(sys.argv[2])
  sizex = int(sys.argv[3])
  sizey = sizex
  sizez = sizex
      
qx = np.pi * np.linspace(-q_range, q_range, sizex)
qy = np.pi * np.linspace(-q_range, q_range, sizey)
qz = np.pi * np.linspace(-q_range, q_range, sizez)

Qx, Qy, Qz = np.meshgrid(qx, qy, qz, indexing = 'ij')
Qx = np.reshape(Qx, sizex*sizey*sizez)
Qy = np.reshape(Qy, sizex*sizey*sizez)
Qz = np.reshape(Qz, sizex*sizey*sizez)
Q  = np.vstack((Qx,Qy,Qz)).T
q = np.ascontiguousarray(Q) # for the purpose of GPU calculation
q_norm = np.sqrt(q[:,0]**2 + q[:,1]**2 + q[:,2]**2)
num_q  = np.shape(q)[0]

'''
Function to calculate Sq for single q
Get called by the calc_sq kernel
From: Li, D., & Zhang, K. (2021). Unifying the concepts of scattering and structure factor in ordered and disordered samples. Journal of Applied Crystallography, 54(2), 644–660. https://doi.org/10.1107/S1600576721001965
'''
@njit
def calc_single_Sq(q_xyz, R, Np_A):
  '''
  q_xyz: the wave vector
  R: particle coordinates
  Np_A: total number of particles
  '''
  sine   = 0
  cosine = 0
  for k in range(Np_A):
    x  = R[k,0]
    y  = R[k,1]
    z  = R[k,2]
    dot_product = q_xyz[0]*x + q_xyz[1]*y + q_xyz[2]*z
    sine   += math.sin(dot_product)
    cosine += math.cos(dot_product)
  sq = cosine**2 + sine**2 
  return(sq)

'''
Kernel for structure factor calculation
Each thread calculates the summation over one wave vector
'''
@cuda.jit
def calc_Sq(S_q, q, R, Np_A, num_q):
    '''
    Not the most efficient way to calculate Sq, but fast enough at the moment
    S_q   : the calculated structure factor
    q     : the wave vectors
    R     : particle coordinates
    Np_A  : total number of particles involved in the calculation
    num_q : number of wave vectors
    '''
    idx = cuda.threadIdx.x + cuda.blockIdx.x * cuda.blockDim.x
    if idx >= num_q:
        return
    q_xyz = q[idx]
    uq    = calc_single_Sq(q_xyz, R, Np_A)
    S_q[idx] = uq/Np_A


tpb = 512
bpg = int(np.ceil(num_q/tpb))

'''
Time averaged S(q)
'''
tpb_nc  = 512
bpg_nc  = int(np.ceil(num_q/tpb_nc))
Sq_sum = np.zeros(num_q) # this stores the summation of the calculated Sq at each frame
start = timer()
Np_A = Npo

for filename in file_name.split(':'):
  # load the coordinates
  R = np.loadtxt(filename, skiprows=2, usecols=(1,2,3))
  # Initiate Sq array and transfer the data to device memory
  S_q      = np.zeros(num_q)
  dS_q     = cuda.to_device(S_q)
  d_q      = cuda.to_device(q)
  d_R      = cuda.to_device(R)
  # Call the kernel for structure factor calculation
  calc_Sq[bpg,tpb](dS_q, d_q, d_R, Np_A, num_q)
  # Transfer the calculated Sq back to host memory
  S_q      = dS_q.copy_to_host()
  Sq_sum   = Sq_sum + S_q
    
end = (timer() - start)/60
print('Total time elapsed: ', end, ' minutes')
Sq_ave = Sq_sum/len(file_name.split(':')) # average the Sq over all the frames

'''
Print qx, qy, qz, Sq
'''
#np.savetxt(file_name.split(':')[0] + '.x' + str(dox) + '.y' + str(doy) + '.z' + str(doz) + '.Sq_3D.dat', np.c_[q, Sq_ave])

import pandas as pd
#create DataFrame
df = pd.DataFrame({'qx': q[:,0],
                   'qy': q[:,1],
                   'qz': q[:,2],
                   'Sq': Sq_ave
                   })
z_avg = df.groupby(['qx', 'qy']).mean().reset_index()
y_avg = df.groupby(['qx', 'qz']).mean().reset_index()
x_avg = df.groupby(['qy', 'qz']).mean().reset_index()
if len(file_name.split(':')) == 1:
  z_avg.to_csv(file_name + '.x1.y1.z.Sq_3D.dat', sep=' ', header=False, index=False)
  y_avg.to_csv(file_name + '.x1.y.z1.Sq_3D.dat', sep=' ', header=False, index=False)
  x_avg.to_csv(file_name + '.x.y1.z1.Sq_3D.dat', sep=' ', header=False, index=False)
else:
  z_avg.to_csv('avg.x1.y1.z.Sq_3D.dat', sep=' ', header=False, index=False)
  y_avg.to_csv('avg.x1.y.z1.Sq_3D.dat', sep=' ', header=False, index=False)
  x_avg.to_csv('avg.x.y1.z1.Sq_3D.dat', sep=' ', header=False, index=False)
