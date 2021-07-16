import sys
import os
import inspect
import json

import multiprocessing
import numpy as np

import scipy.optimize
import scipy as sp
import scipy.interpolate

import pandas as pd

from matplotlib import pyplot as pl 
from function_as_script import scriptify

from closure_measurements.process_dic import load_dgs
from closure_measurements.process_dic import Calc_CTODs as calc_CTODs_function
from closure_measurements.process_dic import CalcInitialModel as CalcInitialModel_function
from closure_measurements.process_dic import InitializeFullModel as InitializeFullModel_function
from closure_measurements.process_dic import CalcFullModel as CalcFullModel_function
from closure_measurements.process_dic import calculate_closureprofile as calculate_closureprofile_function
from closure_measurements.process_dic import TestRegistration

#from crackclosuresim2 import inverse_closure,solve_normalstress
#from crackclosuresim2 import Tada_ModeI_CircularCrack_along_midline

from crackclosuresim2.crackclosure import inverse_closure2
from crackclosuresim2.crackclosure import solve_normalstress
from crackclosuresim2.crackclosure import Tada_ModeI_CircularCrack_along_midline
from crackclosuresim2.crackclosure import ModeI_throughcrack_CODformula
from crackclosuresim2.crackclosure import save_closurestress

import pyopencl as cl
#Sets parameters for plotting 
#pl.rcParams['text.usetex']=True
pl.rcParams.update({'font.size' : 20})
pl.rcParams.update({"text.usetex": True,"font.family": "sans-serif","font.sans-serif": ["Helvetica"]})
pl.rcParams['lines.linewidth'] = 3.0
#Calc_CTODs=scriptify(calc_CTODs_function)
#CalcInitialModel=scriptify(CalcInitialModel_function)
#CalcFullModel=scriptify(CalcFullModel_function)
Calc_CTODs=calc_CTODs_function
CalcInitialModel=CalcInitialModel_function
InitializeFullModel=InitializeFullModel_function
CalcFullModel=CalcFullModel_function
calculate_closureprofile=calculate_closureprofile_function


# Probably want to run view_dic_input on the same data file
# prior to running this to set TipCoords1 and 2 and XRange.
#
# using view_dic_input you can click on points in the plots
# and it will print out the coordinates in meters,
# suitable for use in

if __name__=="__main__":

        #new DIC data for C18-AFVT-013L which is the sample used to generate most of the plots used in the crack closure paper
	#dgsfilename = "/databrowse/AFRLvibro2016/fatigue/Thermal_Optical_data/C18-AFVT-013L_optical_files/0001-C18-AFVT-013L_optical_collect_optical_data_dic.dgs.bz2"
        
        #Location for synthetic DIC data, there are three different files in two directories, in modeling_3d the files are tension for the surface crack and ideal for the circular penny-shaped crack.  In modeling the file partial provides the best data 
	dgsfilename = "/home/cgiuffre/crackclosure_abaqus_modeling_3d/closure_measurement/synthetic_FEA_data_tension.dgs"
	#dgsfilename = "/home/cgiuffre/crackclosure_abaqus_modeling/closure_measurement/synthetic_FEA_data_partial.dgs"
        
	#Cannot currently remember why these two files are in a seperate group
	#dgsfilename =  "/databrowse/AFRLvibro2016/fatigue/Optical_DIC_data/C18-AFVT-007R_optical_files/0001-C18-AFVT-007R_optical_collect_optical_data_dic.dgs"
	#dgsfilename = "/databrowse/AFRLvibro2016/fatigue/Optical_DIC_data/C18-AFVT-008M_optical_files/0000-C18-AFVT-008M_optical_collect_optical_data-1_dic.dgs.bz2"	
        #Short Crack Titanium Samples grown from FIB notches
	#dgsfilename = "/databrowse/AFRLvibro2016/fatigue/Optical_DIC_data/C18-AFVT-013L_optical_files/0001-C18-AFVT-013L_optical_collect_optical_data_dic.dgs.bz2"
	#dgsfilename = "/databrowse/AFRLvibro2016/fatigue/Optical_DIC_data/C18-AFVT-010W_optical_files/0001-C18-AFVT-010W_optical_collect_optical_data_dic.dgs.bz2"
	#dgsfilename = "/databrowse/AFRLvibro2016/fatigue/Optical_DIC_data/C18-AFVT-011X_optical_files/0001-C18-AFVT-011X_optical_collect_optical_data_dic.dgs"
	#dgsfilename = "/databrowse/AFRLvibro2016/fatigue/Optical_DIC_data/C18-AFVT-014C_optical_files/0000-C18-AFVT-014C_optical_collect_optical_data_dic.dgs.bz2"
	#dgsfilename = "/databrowse/AFRLvibro2016/fatigue/Optical_DIC_data/C18-AFVT-016K_optical_files/0000-C18-AFVT-016K_optical_collect_optical_data_dic.dgs.bz2"
	#dgsfilename = "/databrowse/AFRLvibro2016/fatigue/Optical_DIC_data/C18-AFVT-017G_optical_files/0000-C18-AFVT-017G_optical_collect_optical_data_dic.dgs.bz2"
	#dgsfilename = "/databrowse/AFRLvibro2016/fatigue/Optical_DIC_data/C18-AFVT-018J_optical_files/0000-C18-AFVT-018J_optical_collect_optical_data_dic.dgs"
	
        #Short Crack Inconel Samples grown from FIB notches
	#dgsfilename = "/databrowse/AFRLvibro2016/fatigue/Optical_DIC_data/C18-AFVN-001K_optical_files/0000-C18-AFVN-001K_optical_collect_optical_data_dic.dgs.bz2"
	#dgsfilename = "/databrowse/AFRLvibro2016/fatigue/Optical_DIC_data/C18-AFVN-006X_optical_files/0000-C18-AFVN-006X_optical_collect_optical_data_dic.dgs.bz2"
	#dgsfilename = "/databrowse/AFRLvibro2016/fatigue/Optical_DIC_data/C18-AFVN-008N_optical_files/0000-C18-AFVN-008N_optical_collect_optical_data_dic.dgs.bz2"
	#dgsfilename = "/databrowse/AFRLvibro2016/fatigue/Optical_DIC_data/C18-AFVN-009T_optical_files/0000-C18-AFVN-009T_optical_collect_optical_data_dic.dgs.bz2"
	#dgsfilename = "/databrowse/AFRLvibro2016/fatigue/Optical_DIC_data/C18-AFVN-013X_optical_files/0001-C18-AFVN-013X_optical_collect_optical_data-7_dic.dgs.bz2"
	#dgsfilename = "/databrowse/AFRLvibro2016/fatigue/Optical_DIC_data/C18-AFVN-015G_optical_files/0000-C18-AFVN-015G_optical_collect_optical_data_dic.dgs.bz2"
	#dgsfilename = "/databrowse/AFRLvibro2016/fatigue/Optical_DIC_data/C18-AFVN-016C_optical_files/0001-C18-AFVN-016C_optical_collect_optical_data_dic.dgs.bz2"

    
    if relshift_middleimg_lowerleft_corner_x_ref is not None:
        # Only applies to DIC dgs files generated through dc_process that have additional registration info added!
        TestRegistration(nloads,Xposvecs,u_disps,v_disps,
                         ROI_out_arrays,
                         relshift_middleimg_lowerleft_corner_x_ref=relshift_middleimg_lowerleft_corner_x_ref,
                         relshift_middleimg_lowerleft_corner_x_diff=relshift_middleimg_lowerleft_corner_x_diff,
                         relshift_middleimg_lowerleft_corner_y_ref=relshift_middleimg_lowerleft_corner_y_ref,
                         relshift_middleimg_lowerleft_corner_y_diff=relshift_middleimg_lowerleft_corner_y_diff)
        pass
    
    
    
    ## Plot diagnostics...
    # Should have at least one plot that evaluates
    # overall performance in fitting the entire data set.
    #
    # ... how to collapse it down to 2D?
    # (Could add additional lines to the InitialModel plots...)
    #pl.figure()
    
    pl.show()
    pass
