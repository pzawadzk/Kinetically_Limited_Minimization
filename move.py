from pylada.crystal import neighbors
import pylada.periodic_table as pt
from scipy import optimize
import numpy as np
import random
from str_util import *
import void

def move_atom(structure, move_type = 'hop', seed=0, **kwargs):

    np.random.seed(seed)
    random.seed(seed)

    move_step = kwargs.get("move_step", 1.00)
    move_step /= float(structure.scale)

    move_width = kwargs.get("move_width", 2.00)
    move_width /= float(structure.scale)

    index = kwargs.get("index", None)

    min_dist_AA = kwargs.get("min_dist_AA", 1.0)
    min_dist_AB = kwargs.get("min_dist_AB", 1.0)

    relax = kwargs.get("relax", True)


    if move_type == 'hop':
      if index == None:
        index = random.randrange(0, len(structure))
      hop(structure, index, move_step, move_width)
    elif move_type == 'void_hop':
      void_hop(structure, seed)
    elif move_type == 'void_hop_type':
      void_hop_type(structure)
    else:
      pass

    if relax:
      relax_structure_v1(structure, min_dist_AA, min_dist_AB)

def hop(structure, index=0, move_step=0, move_width=0):

    dir = 2*np.random.random_sample((3))-1
    dir = dir/np.dot(dir,dir)**0.5 

    step = move_step + abs(np.random.normal(0, move_width))
    step *= dir
    structure[index].pos += step

def void_hop(structure,seed):
    point, value  = void.get_rnd_void_center(structure, seed=seed)
    neighs = [n for n in neighbors(structure, 1, point)]
    n[0].pos = point

def void_hop_type(structure,seed):
    point, value  = void.get_rnd_void_center(structure, seed=seed)
    neighs = [n for n in neighbors(structure, 12, point)]
    n0 = neighs[0]
    type0 = n0[0].type
    n_in = 0
    n_o = 0 
    for n in neighs[1:]:
      type_tmp =  n[0].type
      if (n[2] - n0[2])*structure.scale <0.4:
        if type_tmp != 'O':
          n_in +=1 
        if type_tmp == 'O':
          n_o +=1

    for n in neighs:
      type = n[0].type
      if n_in < n_o:
        if type != 'O':
            dp = point - n[0].pos
            ddp = structure.scale*np.dot(dp,dp)**0.5
            dp = impose_pbc(structure.cell, dp)
            ddp1 = structure.scale*np.dot(dp,dp)**0.5
            print type, ddp, ddp1, structure.scale*n[2]
            n[0].pos = point
            break
      else:
          if type == 'O':
            dp = point - n[0].pos
            ddp = structure.scale*np.dot(dp,dp)**0.5
            dp = impose_pbc(structure.cell, dp)
            ddp1 = structure.scale*np.dot(dp,dp)**0.5
            print type, ddp, ddp1, structure.scale*n[2]
            n[0].pos = point
            break

    #print '######moved', n_in, n_o, type, index, n[2]*structure.scale
    return structure



def print_to_xyz(structure, tag):
    str = '%d\n  %s\n' % (len(structure), tag)
    for atom in structure:
      str += '%s \t %.4f \t %.4f \t %.4f\n ' % ( atom.type, atom.pos[0]*structure.scale,  atom.pos[1]*structure.scale, atom.pos[2]*structure.scale)
    return str 

def print_to_file(structure, tag):
    f = open('data.xyz', 'a')
    print >>f, len(structure)
    print >>f, 'itner', tag
    for atom in structure:
      print >>f, ' ', atom.type, atom.pos[0] * structure.scale,  atom.pos[1] * structure.scale, atom.pos[2] * structure.scale
    f.close()

