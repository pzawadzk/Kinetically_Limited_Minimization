import numpy as np
from scipy import optimize
import random
import re
from pylada.crystal import neighbors
from str_util import *
from time import time

def  calc_Dx(point, structure):
  '''Calculates the distnace to the nearest atom'''
  neighs = [n for n in neighbors(structure, 3, point)]
  return - neighs[0][2]

def  get_rnd_void_center(structure, seed):
    '''Finds a random void center'''
    seed = int(1000*time())*seed
    np.random.seed(seed)
    v = np.random.random(3)
    v0 =  np.dot(structure.cell, v)
    print v0
    result  = optimize.fmin(calc_Dx, v0, args = (structure,), maxiter=200, full_output=1, disp=0)
    #result  = optimize.fmin_powell(calc_Dx, v0, args = (), maxiter=30, full_output=1, disp=0)
    #result  = optimize.fmin_cg(calc_Dx, v0, args = (), maxiter=80, full_output=1, disp=1)
    #result  = optimize.fmin_bfgs(calc_Dx, v0, args = (), maxiter=80, full_output=1, disp=1)
    #result  = optimize.fmin(calc_Dx, v0, args = (structure), maxiter=10, full_output=1, disp=0)
    void_center  = result[0]
    void_value  = result[1]
    return void_center, void_value

def  get_void_centers(structure, mesh_width=2.0):
    '''Finds all void centers on a grid'''
    abc = [ v for v in structure.cell]
    ds  = [np.dot(v, v)**0.5 for v in abc]
    mesh_nbins  =  [int(d * structure.scale / mesh_width) for d in ds]
    #print mesh_nbins
    mesh_widths = [ds[i] * structure.scale / mesh_nbins[i] for i in range(3)]
    void_centers= []
    void_values= []
    bins = [np.linspace(0, 1, mesh_nbins[i]) for i in range(3)]
    for x in bins[0]:
      for y in bins[1]:
        for z in bins[2]:
          v0 =  np.dot(structure.cell, np.array([x,y,z]))
          result  = optimize.fmin(calc_Dx, v0, args = (structure,), maxiter=200, full_output=1, disp=0)
          #result  = optimize.fmin_powell(calc_Dx, v0, args = (), maxiter=30, full_output=1, disp=0)
          #result  = optimize.fmin_cg(calc_Dx, v0, args = (), maxiter=80, full_output=1, disp=1)
          #result  = optimize.fmin_bfgs(calc_Dx, v0, args = (), maxiter=80, full_output=1, disp=1)
          #result  = optimize.fmin(calc_Dx, v0, args = (structure), maxiter=10, full_output=1, disp=0)
          void_center_tmp = result[0]
          void_value_tmp = result[1]
          Niter = result[3]
          Niter = result[3]
          new = True
          for i, center in enumerate(void_centers):
            d = structure.scale * distance2_pos(structure.cell, center, void_center_tmp)**0.5
            if d <0.5:
              if void_values[i] > void_value_tmp:
                void_values[i]  = void_value_tmp
                void_centers[i] = void_center_tmp 
                new = False
                break
              else:
                new = False
                break
          if new:
            void_centers.append(void_center_tmp)
            void_values.append(void_value_tmp)
    zipped_voids = zip(void_values, void_centers)
    zipped_voids=sorted(zipped_voids, key=lambda void: void[0])
    return zipped_voids

def analyze_voids(zipped_voids, structure):
    ''' Analysis voids'''
    centers = []
    for value, point in zipped_voids:
      print '######', value*structure.scale
      centers.append(-value*structure.scale)
      neighs = [n for n in neighbors(structure, 7, point)]
      for n in neighs:
        print n[0].type, n[2]*structure.scale
    return centers
