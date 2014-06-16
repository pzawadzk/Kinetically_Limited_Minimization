from pylada.crystal import neighbors
import pylada.periodic_table as pt
import numpy as np
import random
import quantities as pq

def get_index(structure, atom):
  N = len(structure)
  for i in range(N):
    if atom == structure[i]:
      return i
  return -1

def distance2_pos(cell, posA, posB):
  dr = posB - posA
  dr  = impose_pbc(cell, dr)
  return np.dot(dr,dr)

def get_disp(structure, structure_disp):
    r2 = 0
    for index in range(len(structure)):
      dp = structure[index].pos - structure_disp[index].pos 
      dp = impose_pbc(structure.cell, dp)
      r2 += np.dot(dp,dp)**0.5
    return r2 * structure.scale

def impose_pbc(cell, v):
    dp = 1e6
    for i in [-1, 0, 1]:
      for j in [-1, 0, 1]:
        for k in [-1, 0, 1]:
          p_temp = v + np.dot(cell, np.array([i,j,k]))
          dp_temp = np.dot(p_temp, p_temp)
          if dp_temp < dp:
            p = p_temp
            dp = dp_temp
    return p

def same_kind(structure,index_A, index_B):
    if structure[index_A].type == structure[index_B].type:
      return True
    else:
      return False  

def relax_pair_cm(structure, index_A, index_B, move_by, min_dist):
    dir = get_vector(structure, index_A, index_B)
    m_A = getattr(pt, structure[index_A].type).atomic_weight
    m_B = getattr(pt, structure[index_A].type).atomic_weight

    ddir = np.dot(dir,dir)**0.5

    if ddir > min_dist:
      dir = -dir/ddir
    else:
      dir = dir/ddir

    x = m_A/(m_A + m_B)
    structure[index_A].pos -= dir * move_by *  (1-x)
    structure[index_B].pos += dir * move_by * x

def get_vector(structure, index_A, index_B):
    d = structure[index_B].pos - structure[index_A].pos  
    return  d

def move_to_cell(structure, index):
    pos = structure[index].pos
    pos = impose_pbc(structure.cell, pos)
    structure[index].pos = pos

def relax_pair(structure,  index_A, index_B, move_by, min_dist):
    dir = get_vector(structure, index_A, index_B)
    ddir = np.dot(dir,dir)**0.5

    if ddir > min_dist:
      dir = -dir/ddir
    else:
      dir = dir/ddir
    structure[index_B].pos += dir * move_by 

def check_dist(structure, min_dist_AA=1, min_dist_AB=1):
    for index_A in range(len(structure)):
      for index_B in range(index_A+1, len(structure)):
        dr = structure[index_B].pos -  structure[index_A].pos
        dr  = impose_pbc(structure.cell, dr)
        dist = float(structure.scale) * np.dot(dr,dr)**0.5
        if not same_kind(structure, index_A, index_B):
          min_dist =  min_dist_AB
          type = 'AB'
        else:
          min_dist =  min_dist_AA
          type = 'AA'
        if min_dist - dist > 0.02:
          result = 'FAIL'
          print index_A, index_B, 'dist=%.2f'%dist, 'min_dist=%.2f'%min_dist, result
        else:
          result = 'PASS'
    return 1 

def relax_structure_v1(structure,  min_dist_AA=1., min_dist_AB=1., Niter = 500):
    from itertools import  combinations
    from pylada.crystal import neighbors
    indexes = range(len(structure))
    for iter in range(Niter):
      succes = 1
      np.random.shuffle(indexes)
      for index_A in indexes:
            neighs = [n for n in neighbors(structure, 1, structure[index_A].pos)]
            index_B = get_index(structure, neighs[0][0])
            
            if not same_kind(structure, index_A, index_B):
              min_dist =  min_dist_AB /float(structure.scale)
            else:
              min_dist =  min_dist_AA /float(structure.scale)
                
            if neighs[0][-1] < min_dist:
            
                move_by = min_dist - neighs[0][-1]
                if move_by*structure.scale > 0.02:
                  succes = 0
                  relax_pair_cm(structure, index_A, index_B, 1.2*move_by, min_dist) 
                  move_to_cell(structure, index_A)
                  move_to_cell(structure, index_B)
      if succes:
        return 1
    return 0 

'''
def relax_structure_v1(structure,  min_dist_AA=1., min_dist_AB=1., Niter = 500):
    for iter in range(Niter):
      succes = 1
      for index_A in range(len(structure)):
        for index_B in range(index_A+1, len(structure)):
            dr = structure[index_B].pos -  structure[index_A].pos
            dr  = impose_pbc(structure.cell, dr)
            dist = np.dot(dr,dr)**0.5
            if not same_kind(structure, index_A, index_B):
              min_dist =  min_dist_AB /float(structure.scale)
            else:
              min_dist =  min_dist_AA /float(structure.scale)
            move_by = min_dist - dist
            if move_by*structure.scale > 0.02:
              succes = 0
              relax_pair_cm(structure, index_A, index_B, move_by, min_dist) 
              move_to_cell(structure, index_A)
              move_to_cell(structure, index_B)
      if succes:
        return 1
    return 0 


def relax_structure_v2(structure,  min_dist_AA=1., min_dist_AB=1., Niter = 500):
    succes_AA = 1
    succes_AB = 1
    for iter in range(Niter):
      if succes_AA and succes_AB:
       iter2 = 0
      else:
       iter2 = 1
      succes_AA = 1
      succes_AB = 1
      if iter % 2: 
        min_dist =  min_dist_AB /float(structure.scale)
        for index_A in range(len(structure)):
          B_range = range(len(structure))
          B_range.pop(index_A)
          for index_B in B_range:
            if not same_kind(structure,index_A, index_B):
              dr = structure[index_B].pos -  structure[index_A].pos
              dr  = impose_pbc(structure.cell, dr)
              dist = np.dot(dr,dr)**0.5
              move_by = (min_dist - dist)
              if move_by*structure.scale > 0.02:
                succes_AB = 0
                relax_pair_cm(structure, index_A, index_B, move_by, min_dist) 
                move_to_cell(structure, index_A)
                move_to_cell(structure, index_B)
      else:
        min_dist =  min_dist_AA /float(structure.scale)
        for index_A in range(len(structure)):
          B_range = range(len(structure))
          B_range.pop(index_A)
          for index_B in B_range:
            if same_kind(structure, index_A, index_B):
              dr = structure[index_B].pos -  structure[index_A].pos
              dr  = impose_pbc(structure.cell, dr)
              dist = np.dot(dr,dr)**0.5
              move_by = (min_dist - dist)
              if move_by*structure.scale > 0.02:
                succes_AA = 0
                relax_pair_cm(structure, index_A, index_B, move_by, min_dist) 
                move_to_cell(structure, index_A)
                move_to_cell(structure, index_B)
      if succes_AA and succes_AB and iter2:
        return 1
    return 0 

'''
