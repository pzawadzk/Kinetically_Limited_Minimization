from md import SimulatedAnnealing
import rnd_structure
from In2O3 import In2O3_160

in2o3_160 = In2O3_160()
cell = in2o3_160.cell

vasp = Vasp()
""" VASP functional """

vasp.ediff = 1e-6
vasp.lplane     = True
vasp.encut = 320
vasp.nelm = 30
vasp.nelmin  = 3
vasp.addgrid    = True
vasp.lreal = "A"
vasp.algo = "veryfast"
vasp.kpoints    = "Kpoints\n0\nM\n2 2 2\n0 0 0"

vasp.add_param = 'npar', 8
vasp.add_specie = "In", "../pseudos/In"
vasp.add_specie = "O", "../pseudos/O"
vasp.add_specie = "Zn", "../pseudos/Zn"

#melt time [ps]
M = 15
#step time [fs]
ts = 4.0
#quench time [ps]
Q1 = 40.
programs =[ 
{"ibrion": 0, "nsw": int(1000 * M / ts),  "tebeg": 3000, "teend": 3000, "nblock": 50, "smass": 0, "potim": ts, "kpoints": "Kpoints\n1\nreciprocal\n0.25 0.25 0.25 1\n", "encut": 250},
{"ibrion": 0, "nsw": int(1000 * Q1 / ts),  "tebeg": 3000, "teend": 500, "nblock": 50, "smass": 0, "potim": ts, "kpoints": "Kpoints\n1\nreciprocal\n0.25 0.25 0.25 1\n", "encut": 250},
{"ibrion": 1, "nsw": 60,  "isif": 7,"encut": 320, "ediffg": -0.05, "kpoints": "Kpoints\n1\nreciprocal\n0.25 0.25 0.25 1\n"},
{"ibrion": 1, "nsw": 30,  "isif": 1,"encut": 320, "ediffg": -0.05, "kpoints": "Kpoints\n0\nM\n2 2 2\n0.01 0.01 0.01", },
{"ibrion": 1, "nsw": 30,  "isif": 7,"encut": 320, "ediffg": -0.05, "kpoints": "Kpoints\n0\nM\n2 2 2\n0.01 0.01 0.01", },
{"potim": 0.3, "ibrion": 2, "nsw": 100,  "isif": 1,"encut": 320, "ediffg": -0.05, "kpoints": "Kpoints\n0\nM\n2 2 2\n0.01 0.01 0.01", },
]

calc = SimulatedAnnealing(vasp, programs, restart = True)

materials = ["In2O3_ZnO"]

lattices =[]
N_units = 1
x_in2o3 = [0, 8, 16, 24, 32]
#densities
rho_zno = 5.61
rho_in2o3 = 7.18

sc_vol = 1/(1.02)**3
for i_in2o3 in x_in2o3: 
  ni_in2o3 = i_in2o3/N_units
  ni_zno = (160 - i_in2o3*5)/(2*N_units)
  atoms=['In', ] * 2 *ni_in2o3 + ['Zn',]*ni_zno + ['O',] * (3*ni_in2o3 + ni_zno)
  rho = (5 * ni_in2o3 * rho_in2o3 + 2 * ni_zno * rho_zno)/(5 * ni_in2o3 +  2 * ni_zno)
  print ni_in2o3, ni_zno, len(atoms), rho
  for seed  in [1,]:
    lattices.append(rnd_structure.rand(N_units, atoms, rho=rho*sc_vol, seed=seed, name = '%d_%d_' %(M, Q1), cell = cell))

""" Lattice for which to create structures. """
