#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
created on Tue June 16 08:41:00 2020

@author: pierrekawak
"""

import sys
import numpy as np
import matplotlib.pyplot as plt

file_names = sys.argv[1].split(':')

fig = plt.figure(figsize=(3,2.8))
ax = fig.add_subplot(111)
ax.set_yscale('log')
for file_name in file_names:
  f = open(file_name, 'r')
  y = np.loadtxt(file_name, skiprows=0, usecols=(1,))
  x = np.loadtxt(file_name, skiprows=0, usecols=(0,))

  #ymean = np.mean(y)
  #ystd = np.std(y)
  #y = y/np.max(y)
  if len(file_names) == 1:
    ax.plot(x, y)
    ax.scatter(x, y)
  else:
    ax.plot(x, y, label=file_name.split('/')[0])
    ax.scatter(x, y, label=file_name.split('/')[0])
  #ax.axhline(ymean, x[0], x[-1], c='c')
  #ax.axhline(ymean+ystd*1.96/np.sqrt(len(y)), x[0], x[-1], c='r')
  #ax.axhline(ymean-ystd*1.96/np.sqrt(len(y)), x[0], x[-1], c='r')
  #ax.text(x[-1], ymean+ystd*1.96/np.sqrt(len(y)), "95% CI")
  #ax.text(0.7, 1, str(round(ymean, 3)) + " +- " + str(round(ystd*1.96/np.sqrt(len(y)), 3)), fontsize=12, transform=ax.transAxes)
if len(file_names)>1:
  plt.legend(fontsize=8)
ax.set_xlabel('q', fontsize=8)
ax.set_ylabel('S(q)', fontsize=8)
ax.tick_params(axis='both', labelsize=8)
plt.yticks(fontsize=8)
#ax.set_ylim(1e2,1e5)
plt.tight_layout()
fig.savefig("Sq_q.png", dpi = 300)
