import sys
import os

import multiprocessing
import numpy as np

from matplotlib import pyplot as pl 
from function_as_script import scriptify

from closure_measurements.process_dic import load_dgs
from closure_measurements.process_dic import Calc_CTODs as calc_CTODs_function
from closure_measurements.process_dic import CalcInitialModel as CalcInitialModel_function
from closure_measurements.process_dic import InitializeFullModel as InitializeFullModel_function
from closure_measurements.process_dic import CalcFullModel as CalcFullModel_function

import pyopencl as cl



#Calc_CTODs=scriptify(calc_CTODs_function)
#CalcInitialModel=scriptify(CalcInitialModel_function)
#CalcFullModel=scriptify(CalcFullModel_function)
Calc_CTODs=calc_CTODs_function
CalcInitialModel=CalcInitialModel_function
InitializeFullModel=InitializeFullModel_function
CalcFullModel=CalcFullModel_function


# Probably want to run view_dic_input on the same data file
# prior to running this to set TipCoords1 and 2 and XRange.
#
# using view_dic_input you can click on points in the plots
# and it will print out the coordinates in meters,
# suitable for use in

if __name__=="__main__":

    #dgsfilename = "/tmp/C18-AFVT-018J_optical_collect_optical_data_dic.dgs"
    dgsfilename = "/tmp/C18-AFVT-011X_optical_collect_optical_data_dic.dgs"
    
    dic_span=20 # formerly step... this is measured in the scaled piexels
    dic_smoothing_window=3  # formerly window... This is measured in the scaled pixels


    nominal_length=2e-3 # nominal crack length, for nondimensional normalization
    nominal_modulus=100.0e9 # nominal modulus
    nominal_stress=50e6 # nominal stress

    tip_tolerance = 100e-6 # 100 microns

    Symmetric_COD=True # assume a symmetric form for the COD -- appropriate when the data is from surface cracks of length 2a where the center is (roughly) a symmetry point
    
    
    ctx = cl.create_some_context()  # set ctx and dev equal to None in order to disable OpenCL acceleration
    dev = ctx.devices[0]

    (dic_dx,dic_dy,dic_nx,dic_ny,XRangeSize,nloads,Xinivec,Xposvecs,load1,load2,u_disps,v_disps,ROI_out_arrays,CrackCenterX,TipCoords1,TipCoords2,ROI_dic_yminidx,ROI_dic_ymaxidx) = load_dgs(dgsfilename)


    #print(TipCoords1)
    #print(TipCoords1[1])
    #print(TipCoords2)
    #print(TipCoords2[1])




    CTODs = Calc_CTODs(dic_nx,nloads,XRangeSize,Xposvecs,v_disps,ROI_out_arrays,ROI_dic_yminidx,ROI_dic_ymaxidx,dic_span,dic_smoothing_window)

    (InitialModels_side1,
     InitialCoeffs_side1,
     Error_side1,
     npoints_side1,
     XPositions_side1,
     CTODValues_side1) = CalcInitialModel(nloads,CTODs,load1,load2,Xposvecs,CrackCenterX,dic_dy,dic_span,Symmetric_COD,side=1,nominal_length=nominal_length,nominal_modulus=nominal_modulus,nominal_stress=nominal_stress,doplots=True)

    

    (InitialModels_side2,
     InitialCoeffs_side2,
     Error_side2,
     npoints_side2,
     XPositions_side2,
     CTODValues_side2) = CalcInitialModel(nloads,CTODs,load1,load2,Xposvecs,CrackCenterX,dic_dy,dic_span,Symmetric_COD,side=2,nominal_length=nominal_length,nominal_modulus=nominal_modulus,nominal_stress=nominal_stress,doplots=True)

    
    (minload_side1,maxload_side1,seed_param_side1) = InitializeFullModel(load1,load2,TipCoords1,TipCoords2,InitialCoeffs_side1,Error_side1,npoints_side1,XPositions_side1,CTODValues_side1,InitialModels_side1,CrackCenterX,tip_tolerance,Symmetric_COD,side=1,doplots=True)

    (minload_side2,maxload_side2,seed_param_side2) = InitializeFullModel(load1,load2,TipCoords1,TipCoords2,InitialCoeffs_side2,Error_side2,npoints_side2,XPositions_side2,CTODValues_side2,InitialModels_side2,CrackCenterX,tip_tolerance,Symmetric_COD,side=2,doplots=True)

    
    #(full_model_params_side1,full_model_result_side1) = CalcFullModel(load1,load2,InitialCoeffs_side1,Error_side1,npoints_side1,YPositions_side1,CTODValues_side1,InitialModels_side1,CrackCenterY,Symmetric_COD,side=1,minload=minload_side1,maxload=maxload_side1,seed_param=seed_param_side1,nominal_length=nominal_length,nominal_modulus=nominal_modulus,nominal_stress=nominal_stress,doplots=True,opencl_ctx=ctx,opencl_dev=dev)


    #(full_model_params_side2,full_model_result_side2) = CalcFullModel(load1,load2,InitialCoeffs_side2,Error_side2,npoints_side2,YPositions_side2,CTODValues_side2,InitialModels_side2,CrackCenterY,Symmetric_COD,side=2,minload=minload_side2,maxload=maxload_side2,seed_param=seed_param_side2,nominal_length=nominal_length,nominal_modulus=nominal_modulus,nominal_stress=nominal_stress,doplots=True,opencl_ctx=ctx,opencl_dev=dev)


    
    
    ## Plot diagnostics...
    # Should have at least one plot that evaluates
    # overall performance in fitting the entire data set.
    #
    # ... how to collapse it down to 2D?
    # (Could add additional lines to the InitialModel plots...)
    #pl.figure()
    
    pl.show()
    pass
