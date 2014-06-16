import numpy as np
def get_vacuum(structure):
  z = structure.cell[2][2]
  zpos = get_zpos(structure)
  vacuum = float((z - (zpos.max() - zpos.min())) * structure.scale)
  return vacuum

def strip_layer(structure, threshold = 0.8):
  dist_list, index_sorted = get_zdist(structure)

  dist_tot = 0
  for i in range(len(dist_list)):
    dist_tot += dist_list[i]
    if dist_tot > 4.0:
      break
  max_dl = max(dist_list[:i])

  indexes = []
  for i, dist in enumerate(dist_list):
    indexes.append(index_sorted[i])
    if dist > max_dl - threshold:
      break
  layer = []
  zmax = structure[indexes[-1]].pos[2]

  indexes.sort(reverse=True)
  for index in indexes:
    layer.append(structure.pop(index))
  center_structure(structure)

  return layer, dist, zmax

def collect_layers(structure, layers):
  for layer, dist, zmax in layers[::-1]:
    zpos = get_zpos(structure) 
    dz = zpos.min() - zmax - dist
    for atom in layer:
      atom.pos[2] += dz
      structure.append(atom)
    center_structure(structure)

def center_structure(structure):
  z = structure.cell[2][2]
  zpos = get_zpos(structure)
  dz = z/2. - zpos.mean()
  for atom in structure:
    atom.pos[2] += dz

def get_mforce(structure):
  mf = 0
  for atom in structure: 
    if getattr(atom, 'freeze', 'xyz') != '':
      f = abs(atom.force).max()
      if mf < f: 
        mf = f
  return mf

def freeze_atoms(structure, indexes):
  from pylada.crystal import neighbors
  for atom in structure:
    atom.freeze = 'xyz'
  for index in indexes:
    structure[index].freeze = ''

def get_zpos(structure):
  zpos = []
  for atom in structure:
    zpos.append(atom.pos[2])
  return np.array(zpos)

def get_zdist(structure):
  Nat = len(structure)
  zpos = get_zpos(structure)

  index_sorted = sorted(range(Nat), key=lambda k: zpos[k])
  dist_list = []

  for i in range(Nat-1):
    dist_list.append(float((zpos[index_sorted[i+1]] - zpos[index_sorted[i]]) * structure.scale)) 
  return dist_list, index_sorted

def get_freeze_indexes(structure, d_frz): 
  dist_list, index_sorted = get_zdist(structure)
  dist_tot = 0 
  freeze_indexes = []
  for i, dist in enumerate(dist_list):
    if dist_tot < d_frz: 
      freeze_indexes.append(index_sorted[i])
      dist_tot += dist
    else:
      break
  return freeze_indexes


def deposit_atoms(structure, atoms, d_frz, Nclimb = 0):
  from time import time
  seed = int(1000*time())
  np.random.seed(seed)
  id = len(structure) - 1
  indexes = []
  mzs = []
  for atom in structure:
      mzs.append(atom.pos[2])
  mzs.sort()
  mz = np.sum(mzs[-4:])/4.

  for atom in atoms:
    pos = np.dot(structure.cell, np.random.random(3))
    pos[2] = mz - 0.0/float(structure.scale)
    pos = min_Dz(structure, pos, d0=2.5)
    #for nclimb in range(Nclimb):
    #  pos[2] -= 2.0/float(structure.scale)
    #  pos[2] = min_Dz(structure, pos[2], pos[0], pos[1], d0=2.0)
    structure.add_atom(pos[0], pos[1], pos[2], atom)

    id +=1
    indexes.append(id)
  indexes_report = order_atoms(structure, indexes)
  freeze_indexes = get_freeze_indexes(structure, d_frz)
  freeze_atoms(structure, freeze_indexes)

  center_structure(structure)

  return indexes_report

def move_atoms(structure, indexes, d_frz):
  from time import time
  seed = int(1000*time())
  np.random.seed(seed)
  indexes.sort(reverse=True)
  atoms_to_move = []
  for id in indexes:
    at = structure.pop(id)
    atoms_to_move.append(at.type)
  deposit_atoms(structure, atoms_to_move, d_frz)

#def  Dz(z, x, y, structure, k=1, d0=1.5):
def  Dz(point, structure, k=1, d0=1.5):
  from pylada.crystal import neighbors
  neighs = [n for n in neighbors(structure, k, point)]
  d0 /= float(structure.scale)
  value = 0
  for i in range(k):
    value += (abs(neighs[k-1][-1]-d0))**2
  return   value

def  min_Dz(structure, point, d0):
  print 'starting', point
  from scipy import optimize
  result  = optimize.fmin_powell(Dz, point, args = (structure,1,d0), maxiter=200, full_output=1, disp=1)
  center = result[0]
  distance = result[1]
  print 'finish', center
  return center

def order_atoms(structure, indexes=[]):
  from pylada.crystal import  specieset
#  indexes = [i[0] for i in sorted(enumerate(structure), key=lambda a: a[1].type)]
#  structure[:] = sorted(structure, key=lambda a: a.type)
#  return indexes
  species = specieset(structure)
  result = []
  indexes_return = []
  count = 0
  for s in species: 
    for index, atom in enumerate(structure):
      if atom.type != s: continue
      result.append(atom)
      if index in indexes:
        indexes_return.append(count)
      count+=1 
  structure[:] = result
  return indexes_return

'''
  for s in species:
    string += "{0} ".format(len([0 for atom in structure if atom.type == s]))
  inv_cell = matrix(structure.cell).I
  selective_dynamics =\
      any([len(getattr(atom, 'freeze', '')) != 0 for atom in structure])
  if selective_dynamics: string += "\nselective dynamics\ndirect\n"
  else: string += '\ndirect\n'
  for s in species: 
    for atom in structure:
      if atom.type != s: continue
      string += "  {0[0]} {0[1]} {0[2]}"\
                .format(dot(inv_cell, atom.pos).tolist()[0])
      freeze = getattr(atom, 'freeze', '')
      if selective_dynamics:
        string += "  {1} {2} {3}\n"\
                    .format( 'T' if 'x' in freeze  != 0 else 'F',
                             'T' if 'y' in freeze  != 0 else 'F',
                             'T' if 'z' in freeze  != 0 else 'F' )
      else: string += '\n'
'''
