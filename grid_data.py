"""
Reads a SpecFile and creates a data cube
"""

# Python Imports	
import numpy as np		
import matplotlib.pyplot as pl
import subprocess		 
import netCDF4			 

# Line Data Reduction Imports
from lmtslr.spec.spec import *
from lmtslr.ifproc.ifproc import *

from lmtslr.reduction.line_reduction import *

from lmtslr.utils.parser import HandleGridOptions


def main(argv):

    Opts = HandleGridOptions()
    Opts.parse_options(argv,'grid_data',1,True)
    
    # check to see whether output file exists and remove it if it does
    if os.path.isfile(Opts.output_file_name) == True:
        os.remove(Opts.output_file_name) 

    with open('out.txt','w+') as outputfile:
        with open('err.txt','w+') as errorfile:
            exit_code=subprocess.call([Opts.program_path,
                                       '-i',Opts.input_file_name,
                                       '-o',Opts.output_file_name,
                                       '-f',str(Opts.otf_select),
                                       '-l',str(Opts.resolution),
                                       '-c',str(Opts.cell),
                                       '-z',str(Opts.rms_cut),
                                       '-n',str(Opts.n_samples),
                                       '-r',str(Opts.rmax),
                                       '-0',str(Opts.otf_a),
                                       '-1',str(Opts.otf_b),
                                       '-2',str(Opts.otf_c),
                                       '-x',str(Opts.x_extent),
                                       '-y',str(Opts.y_extent)],
                                      stdout=outputfile,
                                      stderr=errorfile)
            # reset stdout file to read from it
            outputfile.seek(0)
            # save output (if any) in variable
            standard_output=outputfile.read()
            print('STDOUT ****************************')
            print(standard_output)

            # reset stderr file to read from it
            errorfile.seek(0) 
            # save errors (if any) in variable
            standard_error = errorfile.read()
            print('STDERR ****************************')
            print(standard_error)

    print('Exit Code: %d'%(exit_code))



main(sys.argv[1:])

