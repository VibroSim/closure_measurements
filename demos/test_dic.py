import sys
import os

import multiprocessing
import numpy as np
#from matplotlib import pyplot as pl 
from function_as_script import scriptify
from closure_measurements.perform_dic import execute_dic as execute_dic_function
execute_dic=scriptify(execute_dic_function)


# Probably want to run view_dic_input on the same data file
# prior to running this to set TipCoords1 and 2 and XRange.
#
# using view_dic_input you can click on points in the plots
# and it will print out the coordinates in meters,
# suitable for use in

processpool = multiprocessing.Pool(multiprocessing.cpu_count()/2+1)
#processpool=None

if __name__=="__main__":

    dgdfilename = "/tmp/C14-UTCB-004F_tortuosity_2014-09-17_collect_optical_data-0027.dgd"
    
    dic_scalefactor=5
    dic_radius=20 # measured (I think) in the unscaled original pixels
    
    
    #dic_span=20 # formerly step... this is measured in the scaled piexels
    #dic_smoothing_window=3  # formerly window... This is measured in the scaled pixels
    
    TipCoords1=(0.000339087,0.00317911) # should have smaller value of y
    TipCoords2=(0.000375043,0.00690454) # Should have larger value of y
    XRange=(.15e-3,.8e-3)
    
    tmpdir='/tmp'

    
    
    dgdbasename=os.path.split(dgdfilename)[1]
    dgs_outfilepart=os.path.splitext(dgdbasename)[0]+"_dic.dgs"
    dgs_outfilename=os.path.join(tmpdir,dgs_outfilepart)
    
    execute_dic(dgdfilename,dgs_outfilename,dic_scalefactor,dic_radius,TipCoords1,TipCoords2,XRange,n_threads=4,processpool=processpool,debug=True)

    pass
