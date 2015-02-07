import numpy as np
import re
import sys
import pickle

def impose_pbc(v, cell):
  '''Imposes periodic boundary condtions on a vector'''
  dp = np.dot(v,v)
  p = v.copy()
  for i in [-1, 0, 1]:
      for j in [-1, 0, 1]:
        for k in [-1, 0, 1]:
          dv = np.dot(np.array([i,j,k]), cell)
          p_temp = v + dv 
          dp_temp = np.dot(p_temp, p_temp)
          if dp_temp < dp:
            p = p_temp
            dp = dp_temp
  return p

def get_velocities(path=None, time_step=1):
  '''Reads velocities from XDATCAR''' 
  if path == None:
    path = 'XDATCAR'

  f  = open(path, 'r')

  pattern = re.compile("\s*([-+]?\d{0,2}\.\d{0,8})\s+([-+]?\d{0,2}\.\d{0,8})\s+([-+]?\d{0,2}\.\d{0,8})\s*")
  lines = f.readlines()
  name = lines[0][:-1]
  scale = float(lines[1][:-1])
  cell = np.zeros([3,3])
  for i in range(3):
    m = pattern.search(lines[2+i])
    cell[i]= np.array([ float(m.group(1)), float(m.group(2)), float(m.group(3))])

  cell *= scale
  atoms = lines[5].split()
  Natoms = []
  for atom in lines[6].split():
    Natoms.append(int(atom))

  Nat = sum(Natoms)
  print Nat
  N = (len(lines) - 7)/(Nat + 1)
  count = 7

  p0 = np.zeros([Nat, 3], dtype=np.dtype('d'))
  p = np.zeros([Nat, 3], dtype=np.dtype('d'))
  pos = np.zeros([N,Nat, 3], dtype=np.dtype('d'))

  for n in range(N):
    count +=1
    for i in range(Nat):
        m = pattern.search(lines[count])
        p[i] = np.dot(np.array([float(m.group(1)), float(m.group(2)), float(m.group(3))]), cell)
        count +=1
    dp = p-p0
    if n == 0:
      pos[n] = p0
    else:
      pos[n] = dp
      for j in range(Nat):
        dist = np.dot(pos[n][j], pos[n][j])
        if dist>10.:
          pos[n][j] = impose_pbc(pos[n][j], cell)
        dist = np.dot(pos[n][j], pos[n][j])
    p0 = p.copy()
  print time_step
  V = pos/time_step
  return V


def str_template(name, scale, cell):
  from pylada.crystal import Structure
  structure = Structure()
  structure.name   = name
  structure.scale  = scale
  structure.cell = cell
  return structure


def get_structures(path=None):
  '''Reads structures from XDATCAR'''
  if path == None:
    path = 'XDATCAR'

  f  = open(path, 'r')

  pattern = re.compile("\s*([-+]?\d{0,2}\.\d{0,8})\s+([-+]?\d{0,2}\.\d{0,8})\s+([-+]?\d{0,2}\.\d{0,8})\s*")
  lines = f.readlines()
  name = lines[0][:-1]
  scale = 1#float(lines[1][:-1])
  cell = np.zeros([3,3])
  for i in range(3):
    m = pattern.search(lines[2+i])
    cell[i]= np.array([ float(m.group(1)), float(m.group(2)), float(m.group(3))])


  atoms = lines[5].split()
  Natoms = []
  for atom in lines[6].split():
    Natoms.append(int(atom))

  structure = str_template(name, scale, cell)

  Structures = []
  N = (len(lines) - 7)/(sum(Natoms)+1)
  count = 7
  for n in range(N):
    count +=1
    structure = str_template(name, scale, cell)
    for i, atom in enumerate(atoms):
      for j in range(Natoms[i]):
        m = pattern.search(lines[count])
        pos = np.array([float(m.group(1)), float(m.group(2)), float(m.group(3))])
        pos = np.dot(pos, cell)
        structure.add_atom(pos[0], pos[1], pos[2], atom)
        count +=1
    #print structure
    Structures.append(structure)

  return Structures
