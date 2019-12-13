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
        print("Usage: closure_measurement_dic <dgd_file> <tipcoords1> <centercoords> <tipcoords2> <yrange> [ debug ] [ dic_scalefactor ] [ dic_radius ]")
        print("Perform closure measurement on given .dgd optical microscopy file")
        print("Must specify crack tip locations, center, and yrange as tuples. ")
        print("A non-present tip can be specified as \"None\". ")
        print("If centercoords is specified as \"None\" it will be estimated from")
        print("the average of the two tip positions. ")
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
    CrackCenterCoords = ast.literal_eval(args[3])
    #TipCoords2=(0.000375043,0.00690454) # Should have larger value of x
    TipCoords2 = ast.literal_eval(args[4])

    if CrackCenterCoords is None:
        CrackCenterCoords=((TipCoords1[0]+TipCoords2[0])/2.0,(TipCoords1[1]+TipCoords2[1])/2.0)
        pass
    if TipCoords1 is not None:

        if TipCoords1[0] > CrackCenterCoords[0]:
            raise ValueError("First tip coordinate should have lower value of x than crack center")
        pass


    if TipCoords2 is not None:

        if TipCoords2[0] < CrackCenterCoords[0]:
            raise ValueError("Second tip coordinate should have larger value of x than crack center")
        pass
        
    #YRange=(.15e-3,.8e-3)
    YRange = ast.literal_eval(args[5])

    if len(args) > 6:
        debug = bool(ast.literal_eval(args[6]))
        pass

    if len(args) > 7:
        dic_scalefactor = int(args[7])
        pass
    
    if len(args) > 8:
        dic_radius = int(args[8])
        pass
    
    #tmpdir='/tmp'

    
    
    #dgdbasename=os.path.split(dgdfilename)[1]
    dgs_outfilename=os.path.splitext(dgdfilename)[0]+"_dic.dgs"

    if os.path.exists(dgs_outfilename):
        raise ValueError("Output file \"%s\" exists; will not overwrite" % (dgs_outfilename))
    
    execute_dic(dgdfilename,dgs_outfilename,dic_scalefactor,dic_radius,TipCoords1,CrackCenterCoords,TipCoords2,YRange,n_threads=4,processpool=processpool,debug=debug)

    pass
