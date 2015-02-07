import numpy as np
from scipy import optimize
import random
import re
from pylada.crystal import neighbors
from str_util import *

def get_types(structure):
    Natoms = len(structure)
    N = {}
    atom_types = [] 
    for i in range(Natoms):
      at_type = structure[i].type
      if at_type not in atom_types:
        atom_types.append(at_type)
      if at_type not in N.keys():
        N[at_type] = 1
      else:
        N[at_type] += 1

    return atom_types, N

def get_pdf_types(atom_types):
    types_pdf = []
    for i, atom_a in enumerate(atom_types):
      for j, atom_b in enumerate(atom_types):
        types_pdf.append(atom_a  + atom_b)
    return types_pdf

def get_adf_types(atom_types):
    types_ang = []
    for i, atom_a in enumerate(atom_types):
      for j, atom_b in enumerate(atom_types):
        for k, atom_c in enumerate(atom_types):
          types_ang.append(atom_a  + atom_b + atom_c)
    return types_ang

##################################################################################
#PDF
def calc_correlation(structure,  rad, **kwargs):

    scale = float(structure.scale)
    Natoms = len(structure)

    cell = structure.cell.T
    abc = [ v for v in structure.cell]
    ds  = [np.dot(v, v)**0.5 for v in abc]
    Vol =  np.linalg.det(cell)

    atom_types, N = get_types(structure)
    print atom_types, N 
    types_pdf = kwargs.get("types_pdf", get_pdf_types(atom_types))
    Npdfs = len(types_pdf)


    pdf_cut = kwargs.get("pdf_cut", Vol**(1/3.) * scale)/scale
    pdf_nbins  =  len(rad)

    cut2 = 6.25 / (scale**2)
    N_in_rcut = kwargs.get("N_in_rcut", int(1.2 *Natoms*4*np.pi*pdf_cut**3/(3*Vol)))

    pdfs = {}
    count = {}
    for type in types_pdf:
      [A, B] = re.findall('[A-Z][^A-Z]*', type)
      pdfs[A+B]  = np.zeros(pdf_nbins)
      count[A+B] = 0
    if 1:
      A_range = range(Natoms)
      for index_A in A_range:
        A = structure[index_A].type  
        pos_A = structure[index_A].pos
        neighs = [n for n in neighbors(structure, N_in_rcut, pos_A)]
        for i, nB in enumerate(neighs):
          index_B = get_index(structure, nB[0])
          if index_B > index_A:
            B = structure[index_B].type  
            dist = nB[2]*scale
            for Nbin, d in enumerate(rad[1:]):
              if dist < d:
                pdf_width = rad[Nbin-1] - d
                break
            if Nbin < pdf_nbins:
              pdf_width
              weigth0 = 4 * np.pi * pdf_width**3 / Vol
              val = np.zeros(pdf_nbins)
              val[Nbin] = weigth0
              pdfs[A+B] += val
              pdfs[B+A] += val
    for type in types_pdf:
      [A, B] = re.findall('[A-Z][^A-Z]*', type)
      pdfs[A+B] /=  N[B] * N[A]
      for i in range(pdf_nbins):
        pdfs[A+B][i] /=  (i + 1e-5)**2
    return pdfs

def calc_pdfs(structure,  template = False, index = -1, **kwargs):

    scale = float(structure.scale)
    Natoms = len(structure)

    cell = structure.cell.T
    abc = [ v for v in structure.cell]
    ds  = [np.dot(v, v)**0.5 for v in abc]
    Vol =  np.linalg.det(cell)

    atom_types, N = get_types(structure)
    print atom_types, N 
    types_pdf = kwargs.get("types_pdf", get_pdf_types(atom_types))
    Npdfs = len(types_pdf)

    pdf_width = kwargs.get("pdf_width", 0.05)
    smear_flag = kwargs.get("smear_flag", True)
    pdf_sigma = kwargs.get("pdf_sigma", pdf_width)

    pdf_nsigma  = pdf_sigma/pdf_width
    pdf_nsmear = int(3*pdf_nsigma+0.5)

    pdf_cut = kwargs.get("pdf_cut", Vol**(1/3.) * scale)/scale
    pdf_nbins  =  2 *  int(pdf_cut * scale / (2* pdf_width) ) + 1
    pdf_width = pdf_cut / pdf_nbins
    cut2 = 6.25 / (scale**2)
    N_in_rcut = kwargs.get("N_in_rcut", int(1.2 *Natoms*4*np.pi*pdf_cut**3/(3*Vol)))
    weight0 = 4 * np.pi * pdf_width**3 / Vol


    pdfs = {}
    count = {}
    for type in types_pdf:
      [A, B] = re.findall('[A-Z][^A-Z]*', type)
      pdfs[A+B]  = np.zeros(pdf_nbins)
      count[A+B] = 0
    if 1:
      if index >=0:
        A_range = [index]
        weight0 = Natoms * weight0 
      else:
        A_range = range(Natoms)
        weight0 = weight0 
      for index_A in A_range:
        A = structure[index_A].type  
        pos_A = structure[index_A].pos
        neighs = [n for n in neighbors(structure, N_in_rcut, pos_A)]
        for i, nB in enumerate(neighs):
          index_B = get_index(structure, nB[0])
          if index_B > index_A:
            B = structure[index_B].type  
            Nbin =  np.around(nB[2]/pdf_width)
            if Nbin < pdf_nbins:
              #val = gauss1D( 1., Nbin, pdf_nbins, 0, pdf_nsigma)
              if smear_flag:
                val = gauss1D( 1., Nbin, pdf_nbins, pdf_nsmear, pdf_nsigma)
              else:
                val = np.zeros(pdf_nbins)
                val[Nbin] = 1
              pdfs[A+B] += val
              pdfs[B+A] += val
    for type in types_pdf:
      [A, B] = re.findall('[A-Z][^A-Z]*', type)
      pdfs[A+B] /= weight0 * N[B] * N[A]
      for i in range(pdf_nbins):
        pdfs[A+B][i] /=  (i + 1e-5)**2
    r =  pdf_width * np.arange(pdf_nbins) * scale
    return pdfs, r 

def plot_pdfs(pdfs = None, r=None, count=0):
    for type  in pdfs:
        name = type 
        plot_pdf(pdfs[type], r, count, name)

  
def plot_pdf(pdf, r, count, name):
    from matplotlib import rc
    import matplotlib.pyplot as plt
    fig=plt.figure(num=None, figsize=(8,8), dpi=80,  edgecolor='k')

    ax = fig.add_subplot(111)
    ax.set_ylabel('g(r)', fontsize=22)
    ax.set_xlabel('r', fontsize=22)
    ax.set_title('iteration = %d'%count, fontsize=22)

    yy=ax.yaxis.get_ticklabels()
    for i in yy:
      i.set_fontsize(19)

    xx=ax.xaxis.get_ticklabels()
    for i in xx:
      i.set_fontsize(19)
    plt.plot(r, pdf)
#    plt.show()
#    plt.savefig('pdf_%05d_%s.eps'%(count, name))

##################################################################################
def gauss1D( weight, coord, Nbins, Nsmear=1, Nsigma=1):
    G  = np.zeros([Nbins])
    Norm = 1. 
    g = []
    for dx in range(-Nsmear, Nsmear+1):
      vx = coord + dx 
      dist = dx**2
      if dist <= Nsmear**2 and vx>=0 and vx<Nbins:
        g.append(Norm * np.exp(-dist/(2.* Nsigma**2)))

    for dx in range(-Nsmear, Nsmear+1):
      vx = coord + dx 
      dist = dx**2
      if dist <= Nsmear**2 and vx>=0 and vx<Nbins:
        G[vx] = g[Nsmear+dx] * weight/np.sum(g)
    return G

##################################################################################
#CN
def int_cn_from_pdfs(structure, pdfs, **kwargs):
    pdf_width = kwargs.get("pdf_width", 0.05)
    Vol =  np.linalg.det(structure.cell*structure.scale)
    weight0 = 4 * np.pi * pdf_width**3 / Vol
    int_pdf = {}

    atom_types, N = get_types(structure)
    print N
    for key in pdfs.keys():
      [A, B] = re.findall('[A-Z][^A-Z]*', key)
      pdf = pdfs[key]
      int_pdf[key] = np.zeros(len(pdf))
      v = 0
      for i, g in enumerate(pdf):
        v += pdf[i] * weight0 * i**2 * N[B]
        int_pdf[key][i] =  v
    return int_pdf

def get_cn(structure,  cut_off = 2.8, index = -1, **kwargs):

    atom_types, N = get_types(structure)
    types_pdf = kwargs.get("types_pdf", get_pdf_types(atom_types))
    scale = float(structure.scale)

    cn = {}
    max_coord = 15
    for type in types_pdf:
      cn[type] = np.zeros(max_coord)
    for index_A in range(len(structure)):
      A = structure[index_A].type  
      pos_A = structure[index_A].pos
      neighs = [n for n in neighbors(structure, max_coord-1, pos_A)]
      d2_AB_prev = 10.
      count = {}
      for i, nB in enumerate(neighs):
        B = nB[0].type  
        d2_AB = (scale *nB[2])
        delta_dist = d2_AB - d2_AB_prev
        d2_AB_prev = d2_AB
        if d2_AB > cut_off:
          break
        #if delta_dist>0.5:
        #  break
        else:
          if A+B not in count.keys():
            count[A+B] = 0
          count[A+B] += 1
      for key in count.keys():  
        cn[key][count[key]] += 1./N[A]
    return cn

##################################################################################
#ADFS
def calc_adfs( structure, template = False, index = -1, **kwargs):

    atom_types, N = get_types(structure)
    scale = float(structure.scale)
    cell = structure.cell

    adf_width = kwargs.get("adf_width", 1.0)
    adf_sigma = kwargs.get("adf_sigma", 5*adf_width)

    adf_nsigma  = adf_sigma/adf_width #does not have to be integer
    adf_nsmear = int(3*adf_nsigma+0.5)
    adf_nbins  =  int(180/adf_width)

    adf_cut = kwargs.get("adf_cut", 3.0)
    adf_cut2 = (adf_cut/scale)**2
    types_adf = kwargs.get("types_adf", get_adf_types(atom_types))

    Natoms = len(structure)

    adfs = {}
    for type in types_adf:
      [A, B, C] = re.findall('[A-Z][^A-Z]*', type)
      adfs[type] = np.zeros(adf_nbins)
    N = len(structure)
    if template == False:
      A_range = range(N)
      for index_A in A_range:
        A = structure[index_A].type
        B_range=range(index_A+1,N)
        for index_B in B_range:
            B = structure[index_B].type
            C_range=range(index_B+1,N)
            d2_AB = distance2_pos(cell, structure[index_A].pos, structure[index_B].pos)
            if d2_AB <  4*adf_cut2:
              for index_C in C_range:
                C = structure[index_C].type
                d2_AC = distance2_pos(cell, structure[index_A].pos, structure[index_C].pos)
                if d2_AC <  4*adf_cut2:
                  d2_BC = distance2_pos(cell, structure[index_B].pos, structure[index_C].pos)
                  if d2_BC <  4*adf_cut2:
                    cos_A = 0.9999 * (-d2_BC + d2_AB + d2_AC)/(2*(d2_AB*d2_AC)**0.5)
                    cos_B = 0.9999 * (-d2_AC + d2_AB + d2_BC)/(2*(d2_AB*d2_BC)**0.5) 
                    alpha= adf_nbins * np.arccos(cos_A)/np.pi
                    beta= adf_nbins * np.arccos(cos_B)/np.pi
                    gamma = adf_nbins - alpha  - beta
                    if d2_AC < adf_cut2 and d2_AB < adf_cut2:
                      adfs[B+A+C] +=  gauss1D( 1, int(alpha), adf_nbins, adf_nsmear, adf_nsigma)
                    if d2_BC < adf_cut2 and d2_AB < adf_cut2:
                      adfs[A+B+C] += gauss1D( 1, int(beta), adf_nbins, adf_nsmear, adf_nsigma)
                    if d2_AC < adf_cut2 and d2_BC < adf_cut2:
                      adfs[A+C+B] +=  gauss1D( 1, int(gamma), adf_nbins, adf_nsmear, adf_nsigma)
        return adfs
      else:
        return adfs 



'''
  class  Analyzer(object):
  """ 
  paramters:


  """

  def __init__(self, structure, **kwargs):
    self.structure = structure
    self.scale = float(structure.scale)
    self.Natoms = len(self.structure)

    self.cell = structure.cell.T
    abc = [ v for v in structure.cell]
    self.ds  = [np.dot(v, v)**0.5 for v in abc]
    self.Vol =  np.linalg.det(self.cell)

##################################################################################
    self.pdf_mode =  kwargs.get("pdf_mode", False)
    self.adf_mode =  kwargs.get("adf_mode", False)
    Nadfs = len(types_ang)
##################################################################################
#Types

##################################################################################
#ADF

    if self.adf_mode:
      self.adfs =  self.get_adfs(structure, template = False)
      if self.start_structure:
        self.adfs_template =  self.get_adfs(self.start_structure)
      else:
        self.adfs_template =  self.get_adfs(structure, template = True)
##################################################################################
#Mesh
    mesh_width = kwargs.get("mesh_width", 2.0)
    self.mesh_nbins  =  [int(d * self.scale / mesh_width) for d in self.ds]
    print self.mesh_nbins
    print self.scale, self.cell
    self.mesh_widths = [self.ds[i] * self.scale / self.mesh_nbins[i] for i in range(3)]
   # self.zipped_voids = self.get_void_centers(self.structure)
   # self.void_dict = self.analyze_voids(self.zipped_voids)
##################################################################################

#METHODS
#WRITE METHOD THAT WRITES DOWN ALL THE PARAMTERS


  def get_avg_adf(self):
    Nadf = {}
    adf = {}
    for type in self.types_ang:
      adf[type] = []
      Nadf[type] = 0
    for index_A in range(self.Natoms):
      A = self.structure[index_A].type  
      pos_A = self.structure[index_A].pos
      neighs = [n for n in neighbors(self.structure, 8, pos_A)]
      for i, nB in enumerate(neighs):
        B = nB[0].type  
        for j, nC in enumerate(neighs[i+1:]):
          C = nC[0].type  
          d2_AB = nB[2]**2
          d2_AC = nC[2]**2
          d2_BC = self.distance2(nB[0].pos, nC[0].pos)
          cos_A = 0.9999 * (-d2_BC + d2_AB + d2_AC)/(2*(d2_AB*d2_AC)**0.5)
          cos_B = 0.9999 * (-d2_AC + d2_AB + d2_BC)/(2*(d2_AB*d2_BC)**0.5) 
          alpha = 180*np.arccos(cos_A)/np.pi
          beta= 180*np.arccos(cos_B)/np.pi
          gamma = 180  - alpha  - beta
          if d2_AC < self.cut2 and d2_AB < self.cut2:
            adf[A+B+C].append(alpha)
            adf[A+C+B].append(alpha)
            Nadf[A+B+C] +=1
            Nadf[A+C+B] +=1
          if d2_BC < self.cut2 and d2_AB < self.cut2:
            adf[B+A+C].append(beta)
            adf[B+C+A].append(beta)
            Nadf[B+A+C] +=1
            Nadf[B+C+A] +=1
          if d2_AC < self.cut2 and d2_BC < self.cut2:
            adf[C+A+B].append(gamma)
            adf[C+B+A].append(gamma)
            Nadf[C+A+B] +=1
            Nadf[C+B+A] +=1
    for type in self.types_ang:
        adf[type] = np.array(adf[type])
    return adf 


##################################################################################
  #def get_adfs(self, structure, template = False, index = -1):
  #  adfs = {}
  #  for type in self.types_ang:
  #    [A, B, C] = re.findall('[A-Z][^A-Z]*', type)
  #    adfs[type] = self.get_adf(structure, A, B, C, template, index = index)
  #  adfs = self.get_avg_adf()
  #  return adfs

  #def get_adf(self, structure, type1='A', type2 = 'B', type3 = 'C', template = False, index = -1):

##################################################################################
#CN
  def calc_pdf_extrema(self, pdfs, delta=1.0):
    minima = {}
    for key in pdfs.keys():
      pdf = pdfs[key]
      N = len(pdf)
      minima[key] = [0]
      for i in range(N-2):
        if pdf[i] > pdf[i+1] and pdf[i+2] > pdf[i+1]:
          minima[key].append( (i+1) *  self.scale * self.pdf_width )
    return minima

  def int_cn_form_pdfs(self, pdfs, delta=1.0):
    int_pdf = {}
    for key in pdfs.keys():
      [A, B] = re.findall('[A-Z][^A-Z]*', key)
      pdf = pdfs[key]
      N = len(pdf)
      int_pdf[key] = np.zeros(N)
      v = 0 
      for i, g in enumerate(pdf):
        dr =   self.pdf_width * self.scale
        v += pdf[i] * self.weight0 * i**2 *  self.N[B]
        int_pdf[key][i] =  v
    return int_pdf

  def get_cn2(self):
    cn = {}
    for type in self.types_pdf:
      cn[type] = 0
    for index_A in range(self.Natoms):
      A = self.structure[index_A].type  
      for index_B in range(index_A+1, self.Natoms):
          d2_AB = self.scale**2 * self[2]2(self.structure, index_A, index_B)
          if d2_AB<2.7**2:
            B = self.structure[index_B].type  
            cn[A+B] += 1./(self.N[A])
            cn[B+A] += 1./(self.N[B])
    return cn

  def get_cn(self):
    cn = {}
    for type in self.types_pdf:
      cn[type] = 0
    for index_A in range(self.Natoms):
      A = self.structure[index_A].type  
      pos_A = self.structure[index_A].pos
      neighs = [n for n in neighbors(self.structure, 30, pos_A)]
      for i, nB in enumerate(neighs):
        B = nB[0].type  
        d2_AB = (self.scale *nB[2])
        if d2_AB < 2.7:
          cn[A+B] += 1./self.N[A]
    return cn

  def get_avg_cn(self):
    cn = {}
    for type in self.types_pdf:
      cn[type] = np.zeros(self.Natoms)
    for index_A in range(self.Natoms):
      A = self.structure[index_A].type  
      pos_A = self.structure[index_A].pos
      neighs = [n for n in neighbors(self.structure, 30, pos_A)]
      for i, nB in enumerate(neighs):
        B = nB[0].type  
        d2_AB = nB[2]**2
        if d2_AB < (5*self.scale)**2:
          cn[A+B][index_A] += .5
          cn[B+A][index_B] += .5
    return cn



##################################################################################
#PLOT
  def plot_adfs(self, count=0, template = False):
    if template == False:
      for type  in self.adfs:
        name = type 
        self.plot_adf(self.adfs[type], name, count)
    else:
      for type  in self.adfs_template:
        name = type + '_temp'
        self.plot_adf(self.adfs_template[type], name, count)

  


  def plot_adf(self, adf, name, count):
    from matplotlib import rc
    import matplotlib.pyplot as plt
    fig=plt.figure(num=None, figsize=(8,8), dpi=80,  edgecolor='k')
    rc('text', usetex=True)
    rc('font', serif='Times')

    ax = fig.add_subplot(111)
    ax.set_ylabel('g(theta)', fontsize=22)
    ax.set_xlabel('theta', fontsize=22)
    ax.set_title('iteration = %d'%count, fontsize=22)
    ax.axis(xmin = 0, xmax = 180 , ymin = None, ymax=None)

    yy=ax.yaxis.get_ticklabels()
    for i in yy:
      i.set_fontsize(19)

    xx=ax.xaxis.get_ticklabels()
    for i in xx:
      i.set_fontsize(19)
    x =  180*np.arange(len(adf)) / len(adf)
    plt.plot(x, adf)
    plt.show()
    plt.savefig('adf_%05d_%s.png'%(count, name))


##################################################################################
'''
