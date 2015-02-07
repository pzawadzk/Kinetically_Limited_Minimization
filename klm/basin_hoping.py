
__docformat__ = "restructuredtext en"
__all__ = ['functional', 'basin_hoping']
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
      #structure, energy, value, functional = load(file)
      return Extract(True, outdir, None, None, None, None)
      #return Extract(True, outdir, energy, structure, value, functional)
  
def functional(vasp, structure, programs, func_param, outdir=None, value=False, **kwargs):
  """ A simulated annealing functional """
  from copy import deepcopy
  from pickle import dump
  from random import random
  from pylada.misc import Changedir
  from os.path import join
  from pylada.crystal import vasp_ordered
  from pylada.misc import  RelativePath
  import numpy as np
  import move
  import str_util
  import deposit
  import os
  import shutil

  structure_tmp = deepcopy(structure)
  vasp = deepcopy(vasp)
  outdir = getcwd() if outdir is None else RelativePath(outdir).path
  #structure = vasp_ordered(structure)
  if programs  is None: programs = [{}]

  # TODO: restart feature does not work in pylada
  step_todo = 0
  output = None

  mode  = func_param.get("mode", "rls")
  min_dist_AA = func_param.get("min_dist_AA", 1.0)
  min_dist_AB = func_param.get("min_dist_AB", 1.0)
  T = func_param.get("T", 1.)
  move_type =   func_param.get("move_type", "gauss")
  print 'MOVE type', move_type
  restart = func_param.get("restart", False)

  beta = 1./(8.617e-5 * T)


  #write_parameters(outdir, **kwargs)
  Nat = len(structure)
  if restart:
      for step_todo, program in enumerate(programs): 
        restartdir = join(outdir, str(step_todo))
        with Changedir(restartdir) as pwd:
          last_iter = 0
          if os.path.exists('energy.dat'):
            with open("energy.dat", "r") as file:
              lines = file.readlines()
              last_iter = int(lines[-1][:10])

          iter_todo =  programs[step_todo]["prog_param"]["Nstep"] - last_iter  
          print 'restart info: step_todo, iter_todo', step_todo, iter_todo

          if iter_todo > 0:
            if os.path.exists('accept.dat'):
              with open("accept.dat", "r") as file:
                lines = file.readlines()
                last_accept_iter =  int(lines[-1][:10])
                last_accept_iter_dir = join(outdir, join(join(str(step_todo), "data"), str(last_accept_iter) +'_0' ) )
                with Changedir(last_accept_iter_dir) as pwd:
                  structure_tmp = vasp.Extract().structure
                  E = float(vasp.Extract().energy_sigma0)
                N_min_e = last_accept_iter 
            else:
              N_min_e = 0
              E = 0.0

            programs[step_todo]["prog_param"]["Nstep"] = iter_todo
            with open("tot_iter", "w") as file:
              file.write( "last_iter %d\n" %  last_iter)
              file.write( "iter_todo %d\n" %  iter_todo)
              file.write( "step %d\n" %  step_todo)
              file.write( "Energy %.3f\n" %  E)
            break
          else:
            if os.path.exists('accept.dat'):
              with open("accept.dat", "r") as file:
                lines = file.readlines()
                last_accept_iter =  int(lines[-1][:10])
                last_accept_iter_dir = join(outdir, join(join(str(step_todo), "data"), str(last_accept_iter) +'_0' ) )
                with Changedir(last_accept_iter_dir) as pwd:
                  structure_tmp = vasp.Extract().structure
            with open("tot_iter", "w") as file:
              file.write( "success")
              file.write( "last_iter %d\n" %  last_iter)
              file.write( "iter_todo %d\n" %  iter_todo)
              file.write( "step %d\n" %  step_todo)
  else:
    step_todo = 0
    last_iter = 0
    N_min_e = 0
    E = 0.0
    str_util.relax_structure_v1(structure_tmp, min_dist_AA, min_dist_AB)

  '''
  params = kwargs.copy()
  vasp_tmp = deepcopy(vasp)
  output = vasp_tmp\
               (\
                 structure_tmp,
                 outdir = join(outdir, join(str(step_todo), join("init", str(last_iter) ) )),\
                 relaxation = "static",\
                 **params
               )
  output = vasp_tmp.Extract(join(outdir, join(str(step_todo), join("init", str(last_iter) ) )))
  structure = output.structure
  assert output.success, RuntimeError("VASP run did not succeed.")

  structure_tmp = structure.copy()
  E = float(output.energy_sigma0)
  '''
  disp = 0

  for step, program in enumerate(programs[step_todo:]): 
      print 'entering step', step
      params = kwargs.copy()
      params.update(program["vasp_param"])
      vasp_tmp = deepcopy(vasp)
      if 'add_param' in program.keys():
        for key in program['add_param'].keys():
          if key not in vasp_tmp.__dict__['_input'].keys():
            print key
            vasp_tmp.add_keyword(key, program['add_param'][key])
      if 'vasp_param_final' in program.keys():
        vasp_final = deepcopy(vasp)
        params_final = kwargs.copy()
        params_final.update(program["vasp_param_final"])
        if 'add_param_final' in program.keys():
          for key in program['add_param_final'].keys():
            if key not in vasp_final.__dict__['_input'].keys():
              print key
              vasp_final.add_keyword(key, program['add_param_final'][key])


      for N in range(1, program["prog_param"]["Nstep"]+1):
        dir_out = join(outdir, join(str(step+step_todo), join("data", str(N+last_iter)+'_0' ) ))
        if os.path.isdir(dir_out): shutil.rmtree(dir_out, ignore_errors=True)
        print 'entering iter', N
        output = vasp_tmp\
               (\
                 structure_tmp,\
                 outdir = dir_out,\
                 #relaxation = 'ionic',\
                 #relaxation = 'static',\
                 **params\
               )
        output = vasp_tmp.Extract(dir_out)
        assert output.success, RuntimeError("VASP run did not succeed.")
        if 'vasp_param_final' in program.keys():
          output = vasp_final\
               (\
                 output.structure,\
                 outdir = join(outdir, join(str(step+step_todo), join("data", str(N+last_iter) ) )),\
                 #relaxation = 'ionic',\
                 #relaxation = 'static',\
                 **params_final\
               )
          output = vasp_tmp.Extract(join(outdir, join(str(step+step_todo), join("data", str(N+last_iter)) )) )
          assert output.success, RuntimeError("VASP run did not succeed.")

        fmax = deposit.get_mforce(output.structure)
        E_new = float(output.energy_sigma0)
        dE = E_new - E
        disp_final = str_util.get_disp(structure_tmp, output.structure)
        type_run = 'i' 
        if 1:
          dE = E_new - E
          if dE <0:
            E = float(E_new)
            accept = 1
            structure = output.structure.copy()
            structure_tmp = output.structure.copy()
            structure_min_e = output.structure.copy()
            N_min_e = N+last_iter
          else:
            A_move = np.exp( -beta * dE)
            if A_move > np.random.random():
              E = float(E_new)
              accept = A_move
              structure = output.structure.copy()
              structure_tmp = output.structure.copy()
            else:
              structure_tmp = structure.copy()
              accept = 0
          with Changedir(join(outdir, str(step+step_todo)) ) as pwd:
            #output.structure.write_poscar('poscar_%d'%N )
            with open("energy.dat", "a") as file:
              tag = '%10d %10.5f %10.5f %10.3f %10.3f %10.3f %.2f %s %d\n' %(N + last_iter, E_new/Nat, dE/Nat, disp, disp_final, fmax, accept, type_run, N_min_e)
              file.write(tag) # need to print structure new
            if accept:
              with open("accept.dat", "a") as file:
                tag = '%10d %10.5f %10.5f %10.3f\n' %(N+last_iter, E_new/Nat, dE/Nat, disp_final)
                file.write(tag) # need to print structure new
              with open("traj.xyz", "a") as file:
                tag = 'disp  %.3f' %(disp)
                file.write( move.print_to_xyz(structure, tag)) # need to print structure new
            else:
              #rmtree(join(outdir, join("data", str(N))))
              pass
            with open("accept_all.dat", "a") as file:
              tag = '%10d %10.5f\n' %(N+last_iter, E/Nat)
              file.write(tag) # need to print structure new
          if program["prog_param"]["move_atoms"]:
            if mode == 'rls':
                move.move_atom(structure_tmp, move_type=move_type, seed=N+last_iter, min_dist_AA=min_dist_AA, min_dist_AB=min_dist_AB, relax = True)
                disp = str_util.get_disp(structure_tmp, structure)
            else:
              pass
          else:
            structure_tmp = output.structure.copy()

        if N == program["prog_param"]["Nstep"] and program["prog_param"]["move_atoms"]:
          structure = structure_min_e.copy()
          with Changedir(join(outdir, str(step+step_todo)) ) as pwd:
            with open("N_min_e.dat", "w") as file:
              tag = '%d' %(N_min_e)
              file.write( tag )
        assert output.success, RuntimeError("VASP run did not succeed.")


    # performs final calculation
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

'''
def write_parameters(outdir, **kwargs):
    from lada.opt.changedir import Changedir
    from os.path import join
    import time

    localtime = time.asctime( time.localtime(time.time()) )
    tag = '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n'
    tag += 'Local current time :' + localtime + '\n'
    tag += '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n'
    with Changedir(outdir) as pwd:
      with open("Parameters.dat", "a") as file:
        tag += 'restart =\t'+ str(restart) + '\n'
        tag += 'mode  =\t'+ .mode + '\n'
        tag += 'beta =\t'+ str(beta) + '\n'
        tag += 'T (=1/(beta *k_B) =\t'+ str(T) + 'K\n'
        for key in rls_params:
          tag += key + '=\t'+ str(rls_params[key]) + '\n'
        for i, program in enumerate(programs):
          tag += '############# Program %d #############\n' % i
          for key in program:
            tag += key + '=\t'+ str(program[key]) + '\n'
        file.write(tag) 
'''
