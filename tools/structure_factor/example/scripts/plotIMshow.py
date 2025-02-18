#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
created on Tue June 16 08:41:00 2020

@author: pierrekawak
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
#import matplotlib
#matplotlib.use('tkagg')

num_files = len(sys.argv)-1
for num_file in range(1,num_files+1):
  file_name = sys.argv[num_file]#.split(':')
  fig, ax = plt.subplots()
  
  data = np.loadtxt(file_name)
  qx = data[:, 0] 
  qy = data[:, 1] 
  qz = data[:, 2] 
  
  string = ""
  if '/' in file_name:
    string = file_name.split('/')[-1]
  else:
    string = file_name
  string = string.split('.')
  for i in range(len(string)):
    if string[i] == 'x' or string[i] == 'x1':
      string_loc = i
      break
  stringx = string[string_loc]
  stringy = string[string_loc+1]
  stringz = string[string_loc+2]

  axis = ''
  if stringx == "x":
    axis = 'yz'
    ax.set_xlabel('qz')
    ax.set_ylabel('qy')
    ax.tick_params(axis='both')
    ax.set_title("q = [0, qy, qz]")
  if stringy == "y":
    axis = 'xz'
    ax.set_xlabel('qz')
    ax.set_ylabel('qx')
    ax.tick_params(axis='both')
    ax.set_title("q = [qx, 0, qz]")
  if stringz == "z":
    axis = 'xy'
    ax.set_xlabel('qy')
    ax.set_ylabel('qx')
    ax.tick_params(axis='both')
    ax.set_title("q = [qx, qy, 0]")

  Sq = data[:, 3] 
  size1d = int(round(len(qx)**(1/2)))
  Sq = np.reshape(Sq, (-1, size1d))
  
  im = ax.imshow(Sq, interpolation='none', origin='lower', vmin=0., vmax=2.,
                 extent=[qy[0],qy[-1],qx[0],qx[-1]])#, norm=colors.LogNorm())
  plt.colorbar(im)
  fig.savefig(file_name + ".Sq_plot.png", dpi = 300)
#plt.show()
sys.exit()
