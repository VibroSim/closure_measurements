import sys
import os

import numpy as np

from perform_dic import dic_raw_plots as dic_raw_plots


#processpool=None

def main(args=None):
    if args is None:
        args=sys.argv
        pass

    dic_scalefactor=5
    dic_radius=20 # measured (I think) in the unscaled original pixels
    

    if len(args) < 3:
        print("Usage: closure_measurement_coords <dgd_file>")
        print("Show plots of given dgd file for identifying tip coordinates and XRange")
        sys.exit(0)
        pass

    
    dgdfilename = args[1]
    
    dic_raw_plots(dgdfilename)
    pl.show()
    
    pass
