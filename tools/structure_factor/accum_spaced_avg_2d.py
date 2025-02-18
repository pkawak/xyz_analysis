#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
created on Tue June 16 08:41:00 2020

@author: pierrekawak
"""

import sys
import numpy as np
import matplotlib.pyplot as plt

file_name = sys.argv[1]
outFile_name = sys.argv[2]

comment_lines = 0
col1 = 0
col2 = 1
col3 = 2
col4 = 3

i = 0
x = 0
y = 0
for filei in file_name.split(':'):
#  x = np.loadtxt(filei, skiprows=comment_lines, usecols=(col1,))
  if i == 0:
    x1 = np.loadtxt(filei, skiprows=comment_lines, usecols=(col1,))
    x2 = np.loadtxt(filei, skiprows=comment_lines, usecols=(col2,))
    x3 = np.loadtxt(filei, skiprows=comment_lines, usecols=(col3,))
    y = np.loadtxt(filei, skiprows=comment_lines, usecols=(col4,))
    i += 1
    continue
  y += np.loadtxt(filei, skiprows=comment_lines, usecols=(col4,))
  i += 1
y /= i

np.savetxt(outFile_name, np.c_[x1, x2, x3, y])
def line_prepender(filename, line):
  with open(filename, 'r+') as f:
    content = f.read()
    f.seek(0, 0)
    f.write(line.rstrip('\r\n') + '\n' + content)
#line_prepender(outFile_name, "s cosT")
