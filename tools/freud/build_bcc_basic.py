#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 11:50:17 2020

Uses freud (https://freud.readthedocs.io/en/stable/index.html) to make perfect body-centered cubic (BCC) configuration
@author: pkawak
"""

import sys
import freud

box, points = freud.data.UnitCell.bcc().generate_system(
    num_replicas=5, sigma_noise=0.0
)
N = len(points)
Lx = box.Lx
with open("bcc.freud.xyz", "w") as f:
  f.write(str(N) + "\n")
  f.write("// " + str(N) + " 1 " + str(Lx) + "\n")
  for i in range(N):
    f.write("CH3 %.12f %.12f %.12f 0\n" % (points[i][0], points[i][1], points[i][2]))
