#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  1 11:13:05 2021

@author: qinyu
adapted by pkawak for MCPC/OPWL config file (.xyz) usage
"""

import sys
import os
import json
import numpy as np
from numba import cuda, njit, float64
from timeit import default_timer as timer
from numba.cuda.random import create_xoroshiro128p_states, xoroshiro128p_uniform_float64
import matplotlib.pyplot as plt
import math

file_name = "" #either name of file or number of files delimited by ':'
n_max     = 24 #determines the maximum q (wave vector) as n_max*pi/l
filt      = 0  #determines whether to filter by molecule id (zeroth column in xyz file)
method    = 1  #Uses either method for computing Sq (see functions below)
#        OR          all         #computes for all xyz files in directory with suffix defined in options file\n\
desired_ids = [""]
if len(sys.argv) >= 2:
  if sys.argv[1] == "options":
#read options file (options.calc_sq.json)
    with open("options.calc_sq.json", "r+") as file:
      data = json.load(file)
    file_name = data["file_name"]
    if file_name == "all":
      suffix = data["suffix"]
      files  = sorted(os.listdir())
      files  = [ x for x in files if x.endswith(suffix) ]
      file_name = ':'.join(files)
    n_max     = int(data["n_max"])
    filt      = int(data["filt"])
    method    = int(data["method"])
    if filt:
      desired_ids = data["ids"]
  else:
    file_name = sys.argv[1]
else:
  print("Usage: ./calc_sq.py options     #reads in options from options.calc_sq.json\n\
        OR          <file_name> #computes calc_sq.py on single file\
        ")
  sys.exit()

# load params
with open(file_name.split(':')[0]) as f:
  lines=f.readlines()[1].split(' ')
  Nc   = int(lines[1])
  Nb   = int(lines[2])
  Npo  = Nc*Nb
  lx   = float(lines[3])
  ly   = lx 
  lz   = lx 

# initiate the wave vector
sizex = 2*n_max+1#int(lx) 
sizey = sizex 
sizez = sizex

q_base = np.arange(-n_max, n_max+1)
qx = 2 * np.pi/lx * q_base
qy = 2 * np.pi/ly * q_base
qz = 2 * np.pi/lz * q_base

Qx, Qy, Qz = np.meshgrid(qx, qy, qz, indexing = 'ij')
Qx = np.reshape(Qx, sizex*sizey*sizez)
Qy = np.reshape(Qy, sizex*sizey*sizez)
Qz = np.reshape(Qz, sizex*sizey*sizez)
Q  = np.vstack((Qx,Qy,Qz)).T
q = np.ascontiguousarray(Q) # for the purpose of GPU calculation

#for binning later
q_norm = np.sqrt(q[:,0]**2 + q[:,1]**2 + q[:,2]**2)
del_q      = 2 * np.pi/lx
range_low  = del_q
range_high = del_q * n_max/2.
num_bins   = n_max/2.

#filter stuff that's not well sampled from now
range_min  = range_low
range_max  = range_high + del_q/2
to_remove = np.where((q_norm > range_max) | (q_norm < range_min))
q_norm = np.delete(q_norm, to_remove)
q      = np.delete(q     , to_remove, axis=0)
num_q  = np.shape(q)[0]

'''
method2
Function to calculate Sq for single q
Get called by the calc_sq kernel
From: Pan, J., Chen, C., Li, Y., Wang, L., Tan, L., Li, G., Tang, X., Xiao, L., Lu, J., & Zhuang, L. (2014). Constructing ionic highway in alkaline polymer electrolytes. Energy Environ. Sci., 7(1), 354–360. https://doi.org/10.1039/C3EE43275K
'''
@njit
def calc_single_Sq_Pan2014(q_xyz, R, Np_A, lx, ly, lz):
  '''
  q_xyz: the wave vector
  R: particle coordinates
  Np_A: total number of particles
  '''
  cosine = 0
  for k in range(Np_A):
    xk  = R[k,0]
    yk  = R[k,1]
    zk  = R[k,2]
    for j in range(k+1, Np_A):
      xj = R[j,0]
      yj = R[j,1]
      zj = R[j,2]
      dx = xk - xj
      dy = yk - yj
      dz = zk - zj
      if abs(dx) > lx/2.:
        if dx > lx/2.:
          dx -= lx
        else:
          dx += lx
      if abs(dy) > ly/2.:
        if dy > ly/2.:
          dy -= ly
        else:
          dy += ly
      if abs(dz) > lz/2.:
        if dz > lz/2.:
          dz -= lz
        else:
          dz += lz
      dot_product = q_xyz[0]*dx + q_xyz[1]*dy + q_xyz[2]*dz
      cosine += 2*math.cos(dot_product)
  sq = cosine 
  return(sq)

'''
method1
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
def calc_Sq(S_q, q, R, Np_A, num_q, lx, ly, lz):
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
  if method == 1:
    uq = calc_single_Sq(q_xyz, R, Np_A)
    S_q[idx] = uq/Np_A
  elif method == 2:
    uq = calc_single_Sq_Pan2014(q_xyz, R, Np_A, lx, ly, lz)
    S_q[idx] = uq/Np_A + 1


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
  if filt:
    ids = np.genfromtxt(filename, skip_header=2, usecols=(0,), dtype=str)
   
    ids_bool = np.array([], dtype=int)
    for desired_id in desired_ids:
      ids_bool = np.append(ids_bool, np.where(ids == desired_id)[0])
    R_old = R
    R = R[ids_bool]
    Np_A = len(R)
  # Initiate Sq array and transfer the data to device memory
  S_q      = np.zeros(num_q)
  dS_q     = cuda.to_device(S_q)
  d_q      = cuda.to_device(q)
  d_R      = cuda.to_device(R)
  # Call the kernel for structure factor calculation
  calc_Sq[bpg,tpb](dS_q, d_q, d_R, Np_A, num_q, lx, ly, lz)
  # Transfer the calculated Sq back to host memory
  S_q      = dS_q.copy_to_host()
  Sq_sum   = Sq_sum + S_q
    
end = (timer() - start)/60
print('Total time elapsed: ', end, ' minutes')
Sq_ave = Sq_sum/len(file_name.split(':')) # average the Sq over all the frames

'''
Bin and average over qnorm
'''
#bin stuff
q_means  = np.arange(1, n_max/2+1)*del_q
Sq_means = np.zeros(len(q_means))
count    = np.zeros(len(q_means), dtype=int)
for i in range(len(q_norm)):
  bin_ = int(round((q_norm[i]-2*np.pi/lx)/del_q))
  count[bin_] += 1
  Sq_means[bin_] += Sq_ave[i]
for i in range(len(Sq_means)):
  Sq_means[i] /= count[i]

# plot the curve
plt.figure(figsize=(6.5,4.88))
plt.rc('font', size = 12)
plt.scatter(q_means, Sq_means)
plt.yscale('log')
plt.xlabel(r'$|q|$')
plt.ylabel(r'$S(q)$')
plt.tight_layout()
plt.savefig('Structure_factor_vs_q.png',dpi=300)

#np.savetxt('q.dat', q, fmt='%10e')
#np.savetxt('Sq.dat', Sq_ave, fmt='%10e')
#np.savetxt('q_binned.dat', q_means, fmt='%10e')
#np.savetxt('Sq_binned.dat', Sq_means, fmt='%10e')
if len(file_name.split(':')) == 1:
  np.savetxt(file_name.split(':')[0] + '.Sq_q_binned.dat', np.c_[q_means, Sq_means])
else:
  np.savetxt('avg.Sq_q_binned.dat', np.c_[q_means, Sq_means])
