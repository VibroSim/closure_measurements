import sys
import os

import multiprocessing
import numpy as np

from closure_measurements.process_dic import load_dgs
from closure_measurements.process_dic import Calc_CTODs
from closure_measurements.process_dic import CalcInitialModel
from closure_measurements.process_dic import CalcFullModel

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

    
    # Non-adjustable parameters
    nominal_length=2e-3 # nominal crack length, for nondimensional normalization
    nominal_modulus=100.0e9 # nominal modulus
    nominal_stress=50e6 # nominal stress

    if len(args) < 5:
        print("Usage: closure_measurement_processing <dgs_file> <Symmetric_COD> [DIC_span] [DIC_smoothing_window]")
        print(" ")
        print("Process DIC closure measurement on given .dgs DIC output file")
        print(" ")
        print("Parameters:")
        print("  Symmetric_COD:    True or False indicating whether a symmetric form")
        print("                    should be used for the crack opening displacement.")
        print("                    Symmetric form is appropriate for surface cracks")
        print("                    that are actually of length 2a, but where we are")
        print("                    analyzing half of the crack at a time")
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
        sys.exit(0)
        pass

    
    dgsfilename = args[1]
    Symmetric_COD=bool(args[2])

    if len(args) > 3:
        dic_span=int(args[3])
        pass
    if len(args) > 4:
        dic_smoothing_window=int(args(4))
        pass
    
    
    (dic_dy,dic_dx,dic_ny,dic_nx,YRangeSize,nloads,Yinivec,Yposvecs,load1,load2,u_disps,v_disps,ROI_out_arrays,CrackCenterY,TipCoords1,TipCoords2,ROI_dic_xminidx,ROI_dic_xmaxidx) = load_dgs(dgsfilename)

    CTODs = Calc_CTODs(dic_ny,nloads,YRangeSize,Yposvecs,u_disps,ROI_out_arrays,ROI_dic_xminidx,ROI_dic_xmaxidx,dic_span,dic_smoothing_window)


    
    (InitialModels_side1,
     InitialCoeffs_side1,
     Error_side1,
     npoints_side1,
     YPositions_side1,
     CTODValues_side1) = CalcInitialModel(nloads,CTODs,load1,load2,Yposvecs,CrackCenterY,Symmetric_COD,side=1,nominal_length=nominal_length,nominal_modulus=nominal_modulus,nominal_stress=nominal_stress,doplots=True)


    (InitialModels_side2,
     InitialCoeffs_side2,
     Error_side2,
     npoints_side2,
     YPositions_side2,
     CTODValues_side2) = CalcInitialModel(nloads,CTODs,load1,load2,Yposvecs,CrackCenterY,Symmetric_COD,side=2,nominal_length=nominal_length,nominal_modulus=nominal_modulus,nominal_stress=nominal_stress,doplots=True)

    pass
