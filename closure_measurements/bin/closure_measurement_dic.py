import sys
import os
import ast

import multiprocessing
import numpy as np

import numpy as np
from closure_measurements.perform_dic import execute_dic 

# Probably want to run closure_measurement_coords on the same data file
# prior to running this to set TipCoords1 and 2 and YRange.
#
# using closure_measurement_coords you can click on points in the plots
# and it will print out the coordinates in meters,
# suitable for use in these parameters

processpool = multiprocessing.Pool(multiprocessing.cpu_count()/2+1)
#processpool=None

def main(args=None):
    if args is None:
        args=sys.argv
        pass

    dic_scalefactor=5
    dic_radius=20 # measured (I think) in the unscaled original pixels
    debug = False

    if len(args) < 5:
        print("Usage: closure_measurement_dic <dgd_file> <tipcoords1> <tipcoords2> <yrange> [ debug ] [ dic_scalefactor ] [ dic_radius ]")
        print("Perform closure measurement on given .dgd optical microscopy file")
        print("Must specify crack tip locations and yrange as tuples")
        print("debug, dic_scalefactor and dic_radius are optional parameters")
        print("defaulting to %s, %d and %d respectively" % (str(debug),dic_scalefactor, dic_radius))
        print(" ")
        print("See closure_measurement_coords script for finding crack tip locations")
        sys.exit(0)
        pass

    
    dgdfilename = args[1]
    
    
    #dic_span=20 # formerly step... this is measured in the scaled piexels
    #dic_smoothing_window=3  # formerly window... This is measured in the scaled pixels
    
    #TipCoords1=(0.000339087,0.00317911) # should have smaller value of x
    TipCoords1 = ast.literal_eval(args[2])
    #TipCoords2=(0.000375043,0.00690454) # Should have larger value of x
    TipCoords2 = ast.literal_eval(args[3])

    if TipCoords1[0] > TipCoords2[0]:
        raise ValueError("Second tip coordinate should have larger value of x")
    
    #YRange=(.15e-3,.8e-3)
    YRange = ast.literal_eval(args[4])

    if len(args) > 5:
        debug = bool(ast.literal_eval(args[5]))
        pass

    if len(args) > 6:
        dic_scalefactor = int(args[6])
        pass
    
    if len(args) > 7:
        dic_radius = int(args[7])
        pass
    
    #tmpdir='/tmp'

    
    
    #dgdbasename=os.path.split(dgdfilename)[1]
    dgs_outfilename=os.path.splitext(dgdfilename)[0]+"_dic.dgs"

    if os.path.exists(dgs_outfilename):
        raise ValueError("Output file \"%s\" exists; will not overwrite" % (dgs_outfilename))
    
    execute_dic(dgdfilename,dgs_outfilename,dic_scalefactor,dic_radius,TipCoords1,TipCoords2,YRange,n_threads=4,processpool=processpool,debug=debug)

    pass
