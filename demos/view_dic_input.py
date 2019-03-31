import sys

import numpy as np
from matplotlib import pyplot as pl 
from function_as_script import scriptify

from perform_dic import dic_raw_plots as dic_raw_plots_function
dic_raw_plots=scriptify(dic_raw_plots_function)

# NOTE: You can click on points in the plots and it will print out
# the coordinates in meters

dgdfilename = "/tmp/C14-UTCB-004F_tortuosity_2014-09-17_collect_optical_data-0027.dgd"
dic_raw_plots(dgdfilename)
pl.show()
