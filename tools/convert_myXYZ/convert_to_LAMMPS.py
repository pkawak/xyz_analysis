#!/usr/bin/env python3

import numpy as np
import sys

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
fout.write(str(Nc*(Nb-1)) + " bonds\n")
fout.write(str(Nc*(Nb-2)) + " angles\n")
fout.write(str(0) + " dihedrals\n")
fout.write(str(0) + " impropers\n\n")
fout.write("1 atom types\n")
fout.write("1 bond types\n")
fout.write("1 angle types\n")
fout.write("\n")
fout.write("-" + str(Lx/2.) + " " + str(Lx/2.) + " xlo xhi\n")
fout.write("-" + str(Lx/2.) + " " + str(Lx/2.) + " ylo yhi\n")
fout.write("-" + str(Lx/2.) + " " + str(Lx/2.) + " zlo zhi\n")
fout.write("\n")
fout.write("Masses\n")
fout.write("\n")
fout.write("1 1.0 # BB\n")
fout.write("\n")
fout.write("Atoms # full\n")
fout.write("\n")
for i in range(Nc):
  for j in range(Nb):
    fout.write(str(i*Nb+j+1) + " " + str(i+1) + " " + str(1) + " " + str(0.0) + " " + str(x[i*Nb+j]) + " " + str(y[i*Nb+j]) + " " + str(z[i*Nb+j]) + "\n")
fout.write("\n")
fout.write("Bonds\n")
fout.write("\n")
for i in range(Nc):
  for j in range(Nb-1): 
    fout.write(str(i*(Nb-1)+j+1) + " " + str(1) + " " + str(i*Nb+j+1) + " " + str(i*Nb+j+2) + "\n")
fout.write("\n")
fout.write("Angles\n")
fout.write("\n")
for i in range(Nc):
  for j in range(Nb-2):
    fout.write(str(i*(Nb-2)+j+1) + " " + str(1) + " " + str(i*Nb+j+1) + " " + str(i*Nb+j+2) + " " + str(i*Nb+j+3) + "\n")
fout.write("\n")
