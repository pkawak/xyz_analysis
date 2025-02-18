#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
created on Tue June 16 08:41:00 2020

@author: pierrekawak
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import sys
#import matplotlib
#matplotlib.use('tkagg')

files = ["melt.x1.y1.z.3D.dat", "nematic.x1.y1.z.3D.dat", "cryst.x1.y1.z.3D.dat", \
         "melt.x1.y.z1.3D.dat", "nematic.x1.y.z1.3D.dat", "cryst.x1.y.z1.3D.dat", \
         "melt.x.y1.z1.3D.dat", "nematic.x.y1.z1.3D.dat", "cryst.x.y1.z1.3D.dat"]

num_files = len(files)
for num_file in range(num_files):
  file_name = files[num_file]
  fig, ax = plt.subplots(figsize=(1.15, 1.15))
  
  data = np.loadtxt(file_name)
  qx = data[:, 0] 
  qy = data[:, 1] 
  qz = data[:, 2] 
  
  strr = ""
  if '/' in file_name:
    strr = file_name.split('/')[1]
  else:
    strr = file_name
  strr = strr.split('.')
  for i in range(len(strr)):
    if strr[i] == 'x' or strr[i] == 'x1':
      str_loc = i
      break
  str_x = strr[str_loc]
  str_y = strr[str_loc+1]
  str_z = strr[str_loc+2]
  str_type = strr[str_loc-1]

  axis = ''
  if str_x == "x":
    axis = 'yz'
    ax.set_xticks([])
    ax.set_ylabel(r'$q_y$', fontsize=8, labelpad=1)
    ax.set_xlabel(r'$q_z$', fontsize=8, labelpad=2)
  if str_y == "y":
    axis = 'xz'
    ax.set_xticks([])
    #ax.set_ylabel(r'$q_x$ $(\sigma^{-1})$', fontsize=8, labelpad=2)
    ax.set_xlabel(r'$q_z$', fontsize=8, labelpad=2)
  if str_z == "z":
    axis = 'xy'
    if str_type == "melt":
      ax.set_ylabel(r'$q_x$', fontsize=8, labelpad=2)
    ax.set_xlabel(r'$q_y$', fontsize=8, labelpad=2)

  ax.set_xticks([])
  ax.set_yticks([])

  Sq = data[:, 3] 
  size1d = int(round(len(qx)**(1/2)))
  Sq = np.reshape(Sq, (-1, size1d))
  
  im = ax.imshow(Sq, interpolation='none', origin='lower', vmin=0., vmax=2.,
                 extent=[qy[0],qy[-1],qx[0],qx[-1]])#, norm=colors.LogNorm())
  plt.subplots_adjust(left=0.11, bottom=0.11, right=0.995, top=0.99, wspace=None, hspace=None)
  #plt.colorbar(im)
  fig.savefig("subfig-" + file_name + ".pdf", dpi = 900)
#plt.show()
