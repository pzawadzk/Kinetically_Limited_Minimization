
__docformat__ = "restructuredtext en"
__all__ = ['functional', 'simulated_annealing']
def Extract(outdir=None):
  """ An extraction function for a siumlated annealing functional """
  from os.path import exists
  from os import getcwd
  from collections import namedtuple
  from pickle import load
  from pylada.misc import Changedir

  if outdir == None: outdir = getcwd()
  Extract = namedtuple('Extract', ['success', 'directory', 'energy', 'structure', 'value', 'functional'])
  if not exists(outdir): return Extract(False, outdir, None, None, None, None)
  with Changedir(outdir) as pwd:
    if not exists('OUTCAR'): return Extract(False, outdir, None, None, None, None)
    with open('OUTCAR', 'r') as file:
      structure, energy, value, functional = load(file)
      return Extract(True, outdir, energy, structure, value, functional)
  
def functional(vasp, structure, programs, restart = False, outdir=None, value=False, **kwargs):
  """ A simulated annealing functional """
  from copy import deepcopy
  from pickle import dump
  from random import random
  from pylada.misc import Changedir
  from os.path import join
  from pylada.crystal import vasp_ordered
  from pylada.misc import  RelativePath

  structure = deepcopy(structure)
  vasp = deepcopy(vasp)
  outdir = getcwd() if outdir is None else RelativePath(outdir).path
  #structure = vasp_ordered(structure)
  if programs  is None: programs = [{}]

  # TODO: restart
  step_todo = 0
  output = None

   
  for step, program in enumerate(programs[step_todo:]): 
      params = kwargs.copy()
      params.update(program)
      if 'add_param' in params.keys():
        for key in params['add_param'].keys():
          if key not in vasp.__dict__['_input'].keys():
            print key
            vasp.add_keyword(key, params['add_param'][key])

      output = vasp\
               (\
                structure,
                 outdir = join(outdir, str(step+step_todo)),
                 restart = output,
                 **params
               )
      output = vasp.Extract(join(outdir, str(step+step_todo)))
      structure = output.structure
      assert output.success, RuntimeError("VASP run did not succeed.")

    # performs final calculation
  print programs[0].keys() 
  output = vasp\
             (\
               structure, \
               outdir = outdir,\
               relaxation = "static",\
               restart = output, \
               **kwargs\
             )

  return Extract(outdir)

functional.Extract = Extract

