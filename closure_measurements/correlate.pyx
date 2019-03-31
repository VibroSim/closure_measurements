import multiprocessing
import sys 
import glob #Error mesesage if not included is undeclared name not builtin

from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool
##from libcpp.ptrdiff_t 

# Calling import_array() is required prior to
# using functions in ncorr.pxd
cimport numpy as np
np.import_array()


cimport ncorr    


def correlate(input1,input2,ROI,scalefactor,r,n_threads=multiprocessing.cpu_count(),debug=False):
    cdef ncorr.DIC_analysis_input DIC_input
    cdef ncorr.DIC_analysis_output DIC_output
    #    cdef ncorr.strain_analysis_input strain_input
    #cdef strain_analysis_output strain_output
    cdef vector[ncorr.Image2D] imagevec 
    cdef ncorr.Disp2D disp
    cdef ncorr.difference_type cutoff
    cdef ncorr.Array2D[double] v_c_array
    cdef ncorr.Array2D[double] u_c_array
    cdef ncorr.ROI2D ROI_out
    cdef ncorr.Array2D[bool] ROI_out_c_array
    
    
    imagevec.push_back(ncorr.Image2D(ncorr.numpy2doublearray(input1)))
    imagevec.push_back(ncorr.Image2D(ncorr.numpy2doublearray(input2))) 
     
    
    cutoff = 0
    
    DIC_input = ncorr.DIC_analysis_input(imagevec,
                                         ncorr.ROI2D(ncorr.numpy2boolarray(ROI),cutoff),
                                         int(scalefactor),
                                         ncorr.QUINTIC_BSPLINE_PRECOMPUTE,
                                         ncorr.CIRCLE,
                                         int(r),
                                         n_threads,
                                         ncorr.NO_UPDATE,
                                         debug)
     
    DIC_output = ncorr.DIC_analysis(DIC_input)
    
    # Get the displacement field 
    disp = DIC_output.disps.back()

    # Get the corresponding v and u
    v_c_array = disp.get_v().get_array()
    u_c_array = disp.get_u().get_array()
    ROI_out = disp.get_u().get_roi()
    v_array = ncorr.doublearray2numpy(v_c_array)
     
    u_array = ncorr.doublearray2numpy(u_c_array)

    ROI_out_c_array = ROI_out.get_mask()
    ROI_out_array = ncorr.boolarray2numpy(ROI_out_c_array)
    
    
    return (v_array,u_array,ROI_out_array)
 
