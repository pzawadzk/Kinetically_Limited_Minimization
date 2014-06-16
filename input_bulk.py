import basin_hoping
import rnd_structure
#from In2O3 import In2O3_160

#in2o3_160 = In2O3_160()

mode = 'rls'
move_type = 'void_hop'
move_n = 1
T = 1

vasp = Vasp()
import os
vasp.program = os.path.expanduser('/home/pzawadzk/pyladaExec.sh')

""" VASP functional """
vasp.prec = 'Normal'
vasp.lplane     = True
vasp.addgrid    = True
vasp.ismear = 0
vasp.ediff = 1e-6
vasp.ediffg = -0.05
#vasp.encut = 220 #change
vasp.nelm = 30
vasp.sigma = 0.2 
vasp.lreal = "A"
vasp.nelmin  = 4
vasp.algo = "fast"
vasp.nbands = 920
vasp.kpoints    = "Kpoints\n0\nM\n2 2 2\n0 0 0"

vasp.add_param = "ncore", 16

pseudoDir = "/nopt/nrel/ecom/cid/vasp.pseudopotentials.a/pseudos"

vasp.add_specie = "In", pseudoDir + "/In"
vasp.add_specie = "Sn", pseudoDir + "/Sn"
vasp.add_specie = "O", pseudoDir + "/O"

programs = [
{"add_param": {"potim": 0.20, "ncore": 16}, "vasp_param": {"ediff": 1e-5, "nsw": 30,  "isif": 1, "prec": "low","ibrion": 2,  "kpoints": "Kpoints\n1\nreciprocal\n0.25 0.25 0.25 1\n"}, "prog_param": {"move_atoms": True,"Nstep": 150,} },
{"add_param": {"ncore": 16},"vasp_param": {"ediff": 1e-6,"ibrion": 2, "nsw": 50,  "isif": 1, "ediffg": -0.05, "kpoints": "Kpoints\n0\nM\n1 1 1\n0. 0. 0."}, "prog_param":{"Nstep": 1, "move_atoms": False}},
{"add_param": {"ncore": 16},"vasp_param": {"ediff": 1e-6,"ibrion": 1, "nsw": 50,  "isif": 7, "ediffg": -0.05, "kpoints": "Kpoints\n0\nM\n2 2 2\n0. 0. 0."}, "prog_param":{"Nstep": 1, "move_atoms": False}},
{"add_param": {"ncore": 16},"vasp_param": {"ediff": 1e-6,"ibrion": 1, "nsw": 100,  "isif": 1, "ediffg": -0.05, "kpoints": "Kpoints\n0\nM\n2 2 2\n0. 0. 0."}, "prog_param":{"Nstep": 1, "move_atoms": False}},
]

calc = basin_hoping.functional

func_param = {'mode': mode,  'restart': True,
'min_dist_AA': 2.7, 'min_dist_AB': 1.7,
'move_type':  move_type,
'T': T,
}

Nat = 150
matLatPairs =[]


#n_A = [15,]
n_A = [0,3,5,7,10]
#n_A = [0,4,8,12,16]

key1 = 'In'
key2 = 'Sn'
A = [{'In':2, 'O':3}, 7.18]
B = [{'Sn':1, 'O':2}, 6.95]
import  numpy as np
NA = np.sum([A[0][key] for key in A[0].keys()])
NB = np.sum([B[0][key] for key in B[0].keys()])


sc_vol = 1/(1.02)**3

for i_A in n_A:
  ni_A = i_A*NB
  ni_B = (Nat - ni_A * NA)/NB

  atoms = []
  for key in A[0].keys(): 
    atoms += ni_A * A[0][key] * [ key, ]

  for key in B[0].keys(): 
    atoms += ni_B * B[0][key] * [ key, ]

  rho = (NA * ni_A * A[1] +  NB * ni_B * B[1])/(NA * ni_A + NB * ni_B)
  ratio = 1.0 * A[0][key1]*ni_A / (A[0][key1]*ni_A +  B[0][key2]*ni_B)

  print  i_A,  ratio, rho

  for seed  in [1,2,3,4]:#range(1,11):
    str = rnd_structure.rand(1, atoms, rho=rho*sc_vol, seed=seed, name = 'A_%d_B_%d_T0'%(ni_A, ni_B))
    print len(str)
    matLatPairs.append( ("bulk", str),)


