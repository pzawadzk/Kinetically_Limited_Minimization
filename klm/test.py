""" High-Thoughput of A2BO4 structures. """
__docformat__ = "restructuredtext en"
__all__ = ['simulated_annealing']

def save_job(path, inputpath="input.py", **kwargs):
  """ Jobs to explore possible ground-states. 
  
      :Parameters:
        path 
          Path where the job-folder will be saved. Calculations will be
          performed in the parent directory of this file. Calculations will be
          performed in same directory as this file.
        inputpath
          Path to an input file. Defaults to input.py. 
        kwargs
          Any keyword/value pair to take precedence over anything in the input
          file.

      Creates a high-throughput job-folder to compute the non-magnetic
      ground-state of a host-material.  The new job-folder is loaded into
      memory automatically. No need to call explore. It is save to the path
      provided on input.
  """
  import re
  from IPython.core.interactiveshell import InteractiveShell
  from copy import deepcopy
  from pylada.vasp import read_input
  from pylada.jobfolder import JobFolder
  from pylada import interactive

  from pylada.misc import setBugLev
  setBugLev(0)   # set global debug level
  from pylada.misc import bugLev

  # reads input.
  input = read_input(inputpath)
  input.update(kwargs)

  # Job dictionary.
  jobfolder = JobFolder()

  # loop over material-lattice pairs.
  for (material,lattice) in input.matLatPairs:

    structure = deepcopy(lattice)

    # job folder for this lattice.
    lat_jobfolder = jobfolder / material 

    job = lat_jobfolder / lattice.name / "basin_hoping"
    job.functional = input.calc
    job.params["structure"] = structure
    job.params["ispin"] = 1
    job.params["vasp"] = input.vasp
    job.params["programs"] = input.programs
    job.params["func_param"] = input.func_param
    # saves some stuff for future reference.
    job.material = material
    job.lattice  = lattice
    #print '    test/hi/test: job: ', job
    #print '    test/hi/test: === job.functional ===\n%s\n=== end functional === ' % (job.functional,)


  interactive.jobfolder = jobfolder
  InteractiveShell.instance().magic("savefolders " + path)

