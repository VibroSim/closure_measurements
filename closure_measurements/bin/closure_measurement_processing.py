import sys
import os
import ast

import multiprocessing
import numpy as np

from matplotlib import pyplot as pl

from closure_measurements.process_dic import load_dgs
from closure_measurements.process_dic import Calc_CTODs
from closure_measurements.process_dic import CalcInitialModel
from closure_measurements.process_dic import InitializeFullModel

# Probably want to run closure_measurement_coords on the same data file
# prior to running this to set TipCoords1 and 2 and XRange.
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
    
    Symmetric_COD=True # assume a symmetric form for the COD -- appropriate when the data is from surface cracks of length 2a where the center is (roughly) a symmetry point

    dic_span=20 # formerly step... this is measured in the scaled pixels
    dic_smoothing_window=3  # formerly window... This is measured in the scaled pixels
    tip_tolerance = 100e-6 #Tolerance for the visual measurement of the crack tip location.  This parameter controls how far beyond the manually selected tip location the code is willing to identfy a closure point
    
    # Non-adjustable parameters
    nominal_length=2e-3 # nominal crack length, for nondimensional normalization
    #nominal_modulus=100.0e9 # nominal modulus
    nominal_stress=50e6 # nominal stress

    if len(args) < 3:
        print("Usage: closure_measurement_processing <dgs_file> <Symmetric_COD> <YoungsModulus> [DIC_span] [DIC_smoothing_window] [Tip tolerance]")
        print(" ")
        print("Process DIC closure measurement on given .dgs DIC output file")
        print(" ")
        print("TIP: To use the clickable plots, run through ipython, e.g.:")
        print("      ipython qtconsole &")
        print(" (then within the new qtconsole:)")
        print("      %%run %s <dgs_file> <Symmetric_COD> <YoungsModulus> [DIC_span] [DIC_smoothing_window] [Tip tolerance]" % (args[0]))
        print(" ")
        print("Parameters:")
        print("  Symmetric_COD:    True or False indicating whether a symmetric form")
        print("                    should be used for the crack opening displacement.")
        print("                    Symmetric form is appropriate for surface cracks")
        print("                    that are actually of length 2a, but where we are")
        print("                    analyzing half of the crack at a time")
        print("  YoungsModulus:    Value for Young's modulus, in Pascals. Typically")
        print("                    either 113.8e9 for Ti-6-4 or 200e9 for In718")
        print("  DIC_span:         (optional, default %d)" % (dic_span))
        print("                    Span of DIC perpendicular-to-crack step analysis.")
        print("                    This is measured in the scaled-up pixels used in the")
        print("                    .dgs file and must be large enough combined with the")
        print("                    smoothing window to span both roughness extremes")
        print("                    across the crack horizontal line")
        print("  DIC_smoothing_window: (optional, default %d)" % (dic_smoothing_window))
        print("                    Size of the smoothing window for analysis. This is")
        print("                    measured in the scaled-up pixels used in the .dgs")
        print("                    file.")
        print("  Tip_Tolerance: (optional, unit meters default %g)" %(tip_tolerance))
        print("                    Tolerance for the visual measurement of the crack tip")
        print("                    location.  This parameter controls how far beyond the")
        print("                    manually selected tip location the code is willing to")
        print("                    identfy a closure point")
        sys.exit(0)
        pass

    
    dgsfilename = args[1]
    Symmetric_COD=bool(ast.literal_eval(args[2]))
    YoungsModulus = float(args[3])
    
    if len(args) > 4:
        dic_span=int(args[4])
        pass
    if len(args) > 5:
        dic_smoothing_window=int(args[5])
        pass
    if len(args) > 6:
        tip_tolerance=float(args[6])
        pass
    

    (dic_dx,dic_dy,
     dic_nx,dic_ny,
     XRangeSize,
     nloads,
     Xinivec,Xposvecs,
     load1,load2,u_disps,v_disps,
     ROI_out_arrays,
     CrackCenterX,TipCoords1,TipCoords2,
     ROI_dic_yminidx,ROI_dic_ymaxidx,
     relshift_firstimg_lowerleft_corner_x_ref,
     relshift_firstimg_lowerleft_corner_x_diff,
     relshift_firstimg_lowerleft_corner_y_ref,
     relshift_firstimg_lowerleft_corner_y_diff) = load_dgs(dgsfilename)


    CTODs = Calc_CTODs(dic_nx,nloads,XRangeSize,Xposvecs,v_disps,ROI_out_arrays,ROI_dic_yminidx,ROI_dic_ymaxidx,dic_span,dic_smoothing_window)


    
    (InitialModels_side1,
     InitialCoeffs_side1,
     Error_side1,
     npoints_side1,
     XPositions_side1,
     CTODValues_side1) = CalcInitialModel(nloads,CTODs,
                                          load1,load2,
                                          Xposvecs,CrackCenterX,
                                          dic_dy,dic_span,
                                          Symmetric_COD,1,YoungsModulus,
                                          relshift_firstimg_lowerleft_corner_x_ref=relshift_firstimg_lowerleft_corner_x_ref,
                                          nominal_length=nominal_length,nominal_stress=nominal_stress,
                                          doplots=True)


    (InitialModels_side2,
     InitialCoeffs_side2,
     Error_side2,
     npoints_side2,
     XPositions_side2,
     CTODValues_side2) = CalcInitialModel(nloads,CTODs,
                                          load1,load2,
                                          Xposvecs,CrackCenterX,
                                          dic_dy,dic_span,
                                          Symmetric_COD,2,YoungsModulus,
                                          relshift_firstimg_lowerleft_corner_x_ref=relshift_firstimg_lowerleft_corner_x_ref,
                                          nominal_length=nominal_length,nominal_stress=nominal_stress,
                                          doplots=True)

    
    (minload_side1,maxload_side1,seed_param_side1,fm_plots_side1) = InitializeFullModel(load1,load2,TipCoords1,TipCoords2,InitialCoeffs_side1,Error_side1,npoints_side1,XPositions_side1,CTODValues_side1,InitialModels_side1,CrackCenterX,tip_tolerance,Symmetric_COD,side=1,doplots=True)
    (fitplot_side1,pickableplot_side1,c5plot_side1)=fm_plots_side1

    (minload_side2,maxload_side2,seed_param_side2,fm_plots_side2) = InitializeFullModel(load1,load2,TipCoords1,TipCoords2,InitialCoeffs_side2,Error_side2,npoints_side2,XPositions_side2,CTODValues_side2,InitialModels_side2,CrackCenterX,tip_tolerance,Symmetric_COD,side=2,doplots=True)

    
    print("seed_param_side1=%s" % (str(seed_param_side1)))
    print(" ")
    print("seed_param_side2=%s" % (str(seed_param_side2)))

    pl.show()
    pass
