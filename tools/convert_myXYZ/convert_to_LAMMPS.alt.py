#!/usr/bin/env python3

import numpy as np
import sys

is_bonds = 1
is_angles = 0

def line_prepender(filename, line):
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(line.rstrip('\r\n') + '\n' + content)

file_name = sys.argv[1]
data = np.loadtxt(file_name, skiprows=2, usecols=(1,2,3))
x = data[:,0]
y = data[:,1]
z = data[:,2]
with open(file_name) as f:
  comment_lines = f.readlines()[1]
  Nc = int(comment_lines.split()[1])
  Nb = int(comment_lines.split()[2])
  Lx = float(comment_lines.split()[3])

fout = open(file_name+".lammps", "w+")
fout.write("LAMMPS Description\n\n")
fout.write(str(Nc*Nb)+" atoms\n")
fout.write("8 atom types\n")
if is_bonds:
  fout.write(str(Nc*(Nb-1)) + " bonds\n")
  fout.write("4 bond types\n")
if is_angles:
  fout.write(str(Nc*(Nb-2)) + " angles\n")
  fout.write("1 angle types\n")
#fout.write(str(0) + " dihedrals\n")
#fout.write(str(0) + " impropers\n\n")
fout.write("\n")
fout.write("-" + str(0.) + " " + str(Lx) + " xlo xhi\n")
fout.write("-" + str(0.) + " " + str(Lx) + " ylo yhi\n")
fout.write("-" + str(0.) + " " + str(Lx) + " zlo zhi\n")
fout.write("\n")
fout.write("Masses\n")
fout.write("\n")
#fout.write("1 1.0 # BB\n")
fout.write("1 1.4764\n")
fout.write("2 1.4764\n")
fout.write("3 15.8333\n")
fout.write("4 1\n")
fout.write("5 1\n")
fout.write("6 1\n")
fout.write("7 1\n")
fout.write("8 1\n")
fout.write("\n")
fout.write("Atoms # bond\n")
fout.write("\n")
for i in range(Nc):
  fout.write(str(i*Nb+1) + " " + str(i+1) + " " + str(5) + " " + str(x[i*Nb]) + " " + str(y[i*Nb]) + " " + str(z[i*Nb]) + "\n")
  for j in range(1, Nb-1):
    fout.write(str(i*Nb+j+1) + " " + str(i+1) + " " + str(4) + " " + str(x[i*Nb+j]) + " " + str(y[i*Nb+j]) + " " + str(z[i*Nb+j]) + "\n")
  fout.write(str(i*Nb+Nb-1+1) + " " + str(i+1) + " " + str(5) + " " + str(x[i*Nb+Nb-1]) + " " + str(y[i*Nb+Nb-1]) + " " + str(z[i*Nb+Nb-1]) + "\n")
fout.write("\n")
if is_bonds:
  fout.write("Bonds\n")
  fout.write("\n")
  for i in range(Nc):
    for j in range(Nb-1): 
      fout.write(str(i*(Nb-1)+j+1) + " " + str(4) + " " + str(i*Nb+j+1) + " " + str(i*Nb+j+2) + "\n")
  fout.write("\n")
if is_angles:
  fout.write("Angles\n")
  fout.write("\n")
  for i in range(Nc):
    for j in range(Nb-2):
      fout.write(str(i*(Nb-2)+j+1) + " " + str(1) + " " + str(i*Nb+j+1) + " " + str(i*Nb+j+2) + " " + str(i*Nb+j+3) + "\n")
  fout.write("\n")
