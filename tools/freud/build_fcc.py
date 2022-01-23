#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 11:50:17 2020

@author: pkawak
"""

import sys
from freud_functs import *

box, points = freud.data.UnitCell.fcc().generate_system(
    num_replicas=2, sigma_noise=0.0
)
box.Lx *= 2
box.Ly *= 2
box.Lz *= 2
points[:, 2] *= 2
points[:, 1] *= 2
points = points+0.1
points[:,0] -= 0.5 + 3.625
points[:,1] -= 0.5 + 2.625
points[:,2] -= 0.5 + 2.625
N = len(points)
Lx = box.Lx
with open("fcc.freud.xyz", "w") as f:
  f.write(str(N) + "\n")
  f.write("// " + str(N) + " 1 " + str(Lx) + "\n")
  for i in range(N):
    f.write("CH3 %.12f %.12f %.12f 0\n" % (points[i][0], points[i][1], points[i][2]))

#box_arr = [box]
#pos_arr = [points]
#ax = plot_rdf(box_arr, pos_arr, 'rdf', r_max=1.01, r_min=1.00001)
#print(freud.density.RDF(1000, 1.3))
#plt.show()


#r_max = 1.05
#ids = np.linspace(0, N-1, N, dtype=int)
#print(ids)
#nlist = system.query(points, {'r_max': r_max, 'num_neighbors': num_neighbors, 'exclude_ii': True}).toNeighborList()
#remove all entries with less than 6 neighs
#ToRemove = nlist.neighbor_counts < num_neighbors
#nlist.filter(ToRemove[nlist.query_point_indices]==False)

#print(ToRemove)
#print(nlist.neighbor_counts)
#print(np.count_nonzero( ToRemove ) )

#if np.count_nonzero( [not y for y in ToRemove] ):
#  L = 4
#  wl      = True
#  wl_norm = True
#  avg     = False
#  steinhardt = freud.order.Steinhardt(l=L, wl=wl, wl_normalize=wl_norm, average=avg)
#  steinhardt.compute(system, neighbors=nlist)
#  solid_q6s = np.nan_to_num(steinhardt.particle_order)
##  solid_q6s = steinhardt.particle_order[ [ not y for y in ToRemove ] ]
#  q4 = np.mean(solid_q6s)
#
#  L = 6
#  steinhardt = freud.order.Steinhardt(l=L, wl=wl, wl_normalize=wl_norm, average=avg)
#  steinhardt.compute(system, neighbors=nlist)
#  solid_q6s = np.nan_to_num(steinhardt.particle_order)
##  solid_q6s = steinhardt.particle_order[ [ not y for y in ToRemove ] ]
#  q6 = np.mean(solid_q6s)
#
#  L = 8
#  steinhardt = freud.order.Steinhardt(l=L, wl=wl, wl_normalize=wl_norm, average=avg)
#  steinhardt.compute(system, neighbors=nlist)
#  solid_q6s = np.nan_to_num(steinhardt.particle_order)
##  solid_q6s = steinhardt.particle_order[ [ not y for y in ToRemove ] ]
#  q8 = np.mean(solid_q6s)
#
#  L = 10
#  steinhardt = freud.order.Steinhardt(l=L, wl=wl, wl_normalize=wl_norm, average=avg)
#  steinhardt.compute(system, neighbors=nlist)
#  solid_q6s = np.nan_to_num(steinhardt.particle_order)
##  solid_q6s = steinhardt.particle_order[ [ not y for y in ToRemove ] ]
#  q10 = np.mean(solid_q6s)
#else:
#  q4 = 0.0
#  q6 = 0.0
#  q8 = 0.0
#  q10 = 0.0
#
#q4_str = "%.6f" % q4
#q6_str = "%.6f" % q6
#q8_str = "%.6f" % q8
#q10_str = "%.6f" % q10
#r_max_str = "%.4f" % r_max
#num_neighs_str = "%d" % num_neighbors
#num_str = "%.2f" % np.mean(nlist.neighbor_counts)
#
#print(file_name + " " + num_neighs_str + " " + r_max_str + " " + num_str + " " + q4_str + " " + q6_str + " " + q8_str + " " + q10_str)
#
#def get_ql_specific(num_particles, descriptors, nlist, l, weighted=False):
#    """Given a set of points and a LocalDescriptors object (and the
#    underlying NeighborList), compute the per-particle Steinhardt ql
#    order parameter for all :math:`l` values up to the maximum quantum
#    number used in the computation of the descriptors."""
##    print(descriptors.sph[0:int(descriptors.sph.shape[0]/num_particles), l**2:(l+1)**2])
#    qbar_lm = np.zeros((num_particles, descriptors.sph.shape[1]),
#                       dtype=np.complex128)
#    for i in range(num_particles):
#        indices = nlist.query_point_indices == i
#        Ylms = descriptors.sph[indices, :]
##        if weighted:
##            weights = nlist.weights[indices, np.newaxis]
##            weights /= np.sum(weights)
##            num_neighbors = 1
##        else:
#        weights = np.ones_like(Ylms)
#        num_neighbors = int(np.shape(Ylms)[0])
#        qbar_lm[i, :] = np.sum(Ylms * weights, axis=0)/num_neighbors
#
##    np.set_printoptions(precision=3)
#    np.set_printoptions(suppress=True)
##    print(qbar_lm[0, l**2:(l+1)**2])
#    ql = np.zeros(qbar_lm.shape[0])
#    for i in range(ql.shape[0]):
#      for k in range(l**2, (l+1)**2):
#          ql[i] += np.absolute(qbar_lm[i, k])**2
#     # print(ql[i])
#      ql[i] = np.sqrt(4*np.pi/(2*l + 1) * ql[i])
#
#    return ql

#L = 6
#l_max=6
#ld = freud.environment.LocalDescriptors(l_max, mode='global')
#ld.compute(system, neighbors=nlist);
#ld_ql_L = get_ql_specific(len(points), ld, nlist, L)
#print(np.mean(ld_ql_L))
# Get all vectors from central particles to their neighbors
#rijs = (points[nlist.point_indices] -
#       points[nlist.query_point_indices])
#rijs = box.wrap(rijs)
#for i, j in nlist[:]:
#  print(i, j, np.linalg.norm(points[i]-points[j]), np.linalg.norm(box.wrap(points[i]-points[j])))

#if np.allclose(steinhardt.ql, ld_ql[:, L]):
#    print("Our manual calculation matches the Steinhardt class!")
#with open(file_name + ".ste", "w") as f:
#  for qli in steinhardt.ql:
#    f.write("%s\n" % qli)


#plotSystem(system)

#nop = freud.order.Nematic()
#cl = freud.cluster.Cluster()
#cl.compute(system, neighbors={'r_max': r_max})

#analyzeClusterAndPlot(system, cl)

#clp = freud.cluster.ClusterProperties()
#clp.compute(system, cl.cluster_idx);

#COMandGyr(cl, system, clp)

#GyrationTensorPlot(cl, system, clp)

#print("Num of clusters: " + str(cl.num_clusters))
#print("Radii of gyration and sizes:")
#for cluster_id in range(cl.num_clusters):
#  print(str(clp.radii_of_gyration[cluster_id]) + " " + str(clp.sizes[cluster_id]))

#l_max = 12
#
#ld = freud.environment.LocalDescriptors(l_max, mode='global')
#ld.compute(system, neighbors=nlist);
def get_ql(num_particles, descriptors, nlist, weighted=False):
    """Given a set of points and a LocalDescriptors object (and the
    underlying NeighborList), compute the per-particle Steinhardt ql
    order parameter for all :math:`l` values up to the maximum quantum
    number used in the computation of the descriptors."""
    qbar_lm = np.zeros((num_particles, descriptors.sph.shape[1]),
                       dtype=np.complex128)
    for i in range(num_particles):
        indices = nlist.query_point_indices == i
        Ylms = descriptors.sph[indices, :]
        if weighted:
            weights = nlist.weights[indices, np.newaxis]
            weights /= np.sum(weights)
            num_neighbors = 1
        else:
            weights = np.ones_like(Ylms)
            num_neighbors = descriptors.sph.shape[0]/num_particles
        qbar_lm[i, :] = np.sum(Ylms * weights, axis=0)/num_neighbors

    ql = np.zeros((qbar_lm.shape[0], descriptors.l_max+1))
    for i in range(ql.shape[0]):
        for l in range(ql.shape[1]):
            for k in range(l**2, (l+1)**2):
                ql[i, l] += np.absolute(qbar_lm[i, k])**2
            ql[i, l] = np.sqrt(4*np.pi/(2*l + 1) * ql[i, l])

    return ql

#ld_ql = get_ql(len(points), ld, nlist)
