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

        #3 Seperate stages of crack growth in a titanium sample
	#dgsfilename = "/databrowse/AFRLvibro2016/fatigue/crack_growth_thermography/Crack_Growth_data/C18-AFVT-015E_560micron_optical_files/0000-C18-AFVT-015E_560micron_optical_collect_optical_data_dic.dgs.bz2"
	#dgsfilename = "/databrowse/AFRLvibro2016/fatigue/crack_growth_thermography/Crack_Growth_data/C18-AFVT-015E_750micron_optical_files/0000-C18-AFVT-015E_750micron_optical_collect_optical_data_dic.dgs.bz2"
	#dgsfilename = "/databrowse/AFRLvibro2016/fatigue/crack_growth_thermography/Crack_Growth_data/C18-AFVT-015E_1000micron_optical_files/0000-C18-AFVT-015E_1000micron_optical_collect_optical_data_dic.dgs.bz2"

        #3 Seperate stages of crack growth in an inconel sample
	#dgsfilename = "/databrowse/AFRLvibro2016/fatigue/crack_growth_thermography/Crack_Growth_data/C18-AFVN-003C_560micron_optical_files/0001-C18-AFVN-003C_560micron_optical_collect_optical_data-1_dic.dgs.bz2"
	#dgsfilename = "/databrowse/AFRLvibro2016/fatigue/crack_growth_thermography/Crack_Growth_data/C18-AFVN-003C_750micron_optical_files/0001-C18-AFVN-003C_750micron_optical_collect_optical_data_dic.dgs.bz2"
	#dgsfilename = "/databrowse/AFRLvibro2016/fatigue/crack_growth_thermography/Crack_Growth_data/C18-AFVN-003C_1000micron_optical_files/0000-C18-AFVN-003C_1000micron_optical_collect_optical_data_dic.dgs.bz2"
	
	#Long Cracks grown from a FIB notch in Inconel
	#dgsfilename = "/databrowse/AFRLvibro2016/fatigue/Thermal_Optical_data/C18-AFVN-018A_optical_files/0000-C18-AFVN-018A_optical_collect_optical_data_dic.dgs.bz2"
	#dgsfilename ="/databrowse/AFRLvibro2016/fatigue/Thermal_Optical_data/C18-AFVN-019F_optical_files/0000-C18-AFVN-019F_optical_collect_optical_data_dic.dgs.bz2"
	
        #File location for long cracks grown during the previous effort
	#dgsfilename ="/databrowse/AFRLvibro2016/UTC_specimens/Optical_DIC_data/C14-UTCA-008F_optical_files/0000-C14-UTCA-008F_optical_collect_optical_data-1_dic.dgs.bz2"

	#dgsfilename = "/databrowse/AFRLvibro2016/UTC_specimens/Optical_DIC_data/C14-UTCA-011G_optical_files/0000-C14-UTCA-011G_optical_collect_optical_data_dic.dgs.bz2"
	#Names for the output files used to save the closure stress distribution
	profile_savefile_side1 = 'closure_profiles/UTCA_008F_side1.csv'
	profile_savefile_side2 = 'closure_profiles/UTCA_008F_side2.csv'

	dic_fullmodel_optimization=True #Run the functions to calculate the full integral model.  For an unknown reason this model calls an optimziation that requires the kernel to be reset after every run
	YoungsModulus=117.66e9  # 113.8 GPa for Ti-6-4
	PoissonsRatio=0.32 #Assumned to be the same in both materials
	#YoungsModulus=200.0e9 # 200 GPa for In718 
	
	SampleHalfWidth = 15.16e-3 #12.88e-3 Long #12.88e-3 AFVN #12.98e-3 AFVT
	SampleThickness = 12.8e-3 #12.8e-3 Long #12.8e-3  AFNV #12.9e-3 AFVT

	#sys.modules["__main__"].__dict__.update(globals())
	#sys.modules["__main__"].__dict__.update(locals())

	dic_span=5 # formerly step... this is measured in the scaled piexels 20 for DIC, 5 for FEA
	dic_smoothing_window=3  # formerly window... This is measured in the scaled pixels 5 for DIC, 3 for FEA
	print('Debugging')

	min_dic_points_per_meter=1 #For FEA
	#min_dic_points_per_meter = 40000 #For 

	nominal_length=10e-3 # nominal crack length, for nondimensional normalization
	if dic_fullmodel_optimization:
		nominal_modulus=100.0e9 # nominal modulus
		pass

	nominal_stress=50e6 # nominal stress

	tip_tolerance = 100e-6 # 100 microns Sets the limits of the spline fit of the closure point analysis

	Symmetric_COD=True # assume a symmetric form for the COD -- appropriate when the data is from surface cracks of length 2a where the center is (roughly) a symmetry point.  For any circular crack this is the only apprioprate analysis to run
	
	"""
	#Crack Closure Parameters
	sigma_yield = 880e6
	a=5.0e-3  # half-crack length (m)
	xmax = 6.5e-3 # as far out in x as we are calculating (m)
	xsteps = 200
	crack_middle = 7.5e-03 #coordinates of the crack cetner

	# x_bnd represents x coordinates of the boundaries of
	# each mesh element 
	x_bnd_1=np.linspace(crack_middle,crack_middle-xmax,xsteps,dtype='d')
	x_bnd_2=np.linspace(0,xmax,xsteps,dtype='d')

	
	dx=x_bnd_2[1]-x_bnd_2[0]
	x2 = (x_bnd_2[1:]+x_bnd_2[:-1])/2.0  # x represents x coordinates of the centers of each mesh element
	x1 = (x_bnd_1[1:]+x_bnd_1[:-1])/2.0  # x represents x coordinates of the centers of each mesh element
	crack_model = Tada_ModeI_CircularCrack_along_midline(YoungsModulus,PoissonsRatio)
	"""
	
	if dic_fullmodel_optimization:
		ctx = cl.create_some_context()  # set ctx and dev equal to None in order to disable OpenCL acceleration
		dev = ctx.devices[0]
		print("Using accelerator \"%s\" for fullmodel optimization" % (dev.name))
		pass
	else:
		
		ctx = None
		dev = None
		pass


	(dic_dx,dic_dy,
	 dic_nx,dic_ny,
	 XRangeSize,
	 nloads,
	 Xinivec,Xposvecs,
	 load1,load2,u_disps,v_disps,
	 ROI_out_arrays,
	 CrackCenterCoords,TipCoords1,TipCoords2,
	 ROI_dic_yminidx,ROI_dic_ymaxidx,
	 relshift_middleimg_lowerleft_corner_x_ref,
	 relshift_middleimg_lowerleft_corner_x_diff,
	 relshift_middleimg_lowerleft_corner_y_ref,
	 relshift_middleimg_lowerleft_corner_y_diff) = load_dgs(dgsfilename)


#Crack Closure Parameters
	sigma_yield = 880e6 #Yield Stress of the Material
	a=(TipCoords2[0]-TipCoords1[0])/2  # half-crack length (m)
	xmax = a + 300e-6 # as far out in x as we are calculating (m)
	xsteps = 200
	crack_middle = CrackCenterCoords[0] #coordinates of the crack cetner

	# x_bnd represents x coordinates of the boundaries of
	# each mesh element 
	#x_bnd_1=np.linspace(crack_middle,crack_middle-xmax,xsteps,dtype='d')
	x_bnd_1=np.linspace(0,xmax,xsteps,dtype='d')
	x_bnd_2=np.linspace(0,xmax,xsteps,dtype='d')

	
	dx=x_bnd_2[1]-x_bnd_2[0]
	x2 = (x_bnd_2[1:]+x_bnd_2[:-1])/2.0  # x represents x coordinates of the centers of each mesh element
	x1 = (x_bnd_1[1:]+x_bnd_1[:-1])/2.0  # x represents x coordinates of the centers of each mesh element
	crack_model = Tada_ModeI_CircularCrack_along_midline(YoungsModulus,PoissonsRatio)
	#crack_model = ModeI_throughcrack_CODformula(YoungsModulus,True)
	#print(TipCoords1)
	#print(TipCoords1[1])
	#print(TipCoords2)
	#print(TipCoords2[1])



        #Script to Function Call for the calculation of the crack opeing displacement using virtual displacement gages
	CTODs = Calc_CTODs(dic_nx,nloads,XRangeSize,Xposvecs,v_disps,ROI_out_arrays,ROI_dic_yminidx,ROI_dic_ymaxidx,dic_span,dic_smoothing_window)

	#Calculates the initial/differential model for the left side of the crack as the data is split around the location of the crack center
	(InitialModels_side1,
	 InitialCoeffs_side1,
	 Error_side1,
	 npoints_side1,
	 XPositions_side1,
	 CTODValues_side1) = CalcInitialModel(nloads,CTODs,load1,load2,Xposvecs,CrackCenterCoords,dic_dy,dic_span,Symmetric_COD,1,YoungsModulus,relshift_middleimg_lowerleft_corner_x_ref=relshift_middleimg_lowerleft_corner_x_ref,nominal_length=nominal_length,nominal_stress=nominal_stress,doplots=False)
	

        #Calculates the initial/differential model for the right side of the cracks as the data is split around the locatin of the crack center
	(InitialModels_side2,
	 InitialCoeffs_side2,
	 Error_side2,
	 npoints_side2,
	 XPositions_side2,
	 CTODValues_side2) = CalcInitialModel(nloads,CTODs,load1,load2,Xposvecs,CrackCenterCoords,dic_dy,dic_span,Symmetric_COD,2,YoungsModulus,relshift_middleimg_lowerleft_corner_x_ref=relshift_middleimg_lowerleft_corner_x_ref,nominal_length=nominal_length,nominal_stress=nominal_stress,doplots=False)

	#Generates the spline coefficients to describre the closure point as a funcitonof the average applied load.  Follows the same form as the previous steps using seperate calls for the left and right crack tips.  In this function, the crack depth opening is also determined and a spline fit is attempted.  This has only successfully been done on the FEA data as the experimetnal data lacks the fedility to provide consistent data
	(minload_side1,maxload_side1,seed_param_side1_avg,seed_param_side1_depth,lowest_avg_load_used_side1,fm_plots,fm_plotdata_side1_avg,fm_plotdata_side1_max) = InitializeFullModel(load1,load2,TipCoords1,TipCoords2,InitialCoeffs_side1,YoungsModulus,Error_side1,npoints_side1,XPositions_side1,CTODValues_side1,InitialModels_side1,CrackCenterCoords,tip_tolerance,SampleHalfWidth,SampleThickness,min_dic_points_per_meter,Symmetric_COD,side=1,doplots=False)
	
	(minload_side2,maxload_side2,seed_param_side2_avg,seed_param_side2_depth,lowest_avg_load_used_side2,fm_plots,fm_plotdata_side2_avg,fm_plotdata_side2_max) = InitializeFullModel(load1,load2,TipCoords1,TipCoords2,InitialCoeffs_side2,YoungsModulus,Error_side2,npoints_side2,XPositions_side2,CTODValues_side2,InitialModels_side2,CrackCenterCoords,tip_tolerance,min_dic_points_per_meter,Symmetric_COD,side=2,doplots=True)

        #Refines the spline fit to the closure point measurements
	#fm_plotdata_side2 contains the point data for the initial fit of the data
	#[0] and [1] are all of the data and [2] [3] are the valid data points
	#(output_loads1,tippos_side1_avg,output_loads2,tippos_side2_avg) = calculate_closureprofile(minload_side1,maxload_side1,minload_side2,maxload_side2,nloads,seed_param_side1_avg,seed_param_side2_avg,TipCoords1[0],TipCoords2[0])
	(output_loads1,tippos_side1_avg,output_loads2,tippos_side2_avg) = calculate_closureprofile(minload_side1,maxload_side1,minload_side2,maxload_side2,50,seed_param_side1_avg,seed_param_side2_avg,TipCoords1[0],TipCoords2[0])
	
	(output_loads1_depth,tippos_side1_depth,output_loads2_depth,tippos_side2_depth) = calculate_closureprofile(minload_side1,maxload_side1,minload_side2,maxload_side2,50,seed_param_side1_depth,seed_param_side2_depth,2*TipCoords1[0],2*TipCoords2[0])	
	#(output_loads1,tippos_side1_max,output_loads2,tippos_side2_max) = calculate_closureprofile(minload_side1,maxload_side1,minload_side2,maxload_side2,nloads,seed_param_side1_max,seed_param_side2_max,TipCoords1[0],TipCoords2[0])	
	
        #Plot of the right crack tip closure point behavior as a function of applied load
	pl.figure(2001)
	#pl.plot(fm_plotdata_side2_avg[0]*1.0e6-crack_middle*1e6,fm_plotdata_side2_avg[1]/1.0e6,'x',label="Differential Model Closure Points")
	pl.plot(fm_plotdata_side1_avg[2]*1.0e6-crack_middle*1e6,fm_plotdata_side1_avg[3]/1.0e6,'o')	
	pl.plot(fm_plotdata_side2_avg[2]*1.0e6-crack_middle*1e6,fm_plotdata_side2_avg[3]/1.0e6,'o')	
	#pl.plot(np.abs(fm_plotdata_side1_avg[2]*1.0e6-crack_middle*1e6),fm_plotdata_side1_avg[3]/1.0e6,'rx',label="YZ Plane Closure Points")
	#pl.plot(tippos_side2_avg*1.0e6-crack_middle*1e6,output_loads2/1.0e6,'r-',label='Differential Model Spline')


	#pl.plot(tippos_side2_max*1.0e6-crack_middle*1e6,output_loads2/1.0e6,'-',label="Initial Spline Fit")	

	pl.grid()
        #Statements needed to run the full/integral model.  As with the initial/differential model the crack tips are studied independtly.
	if dic_fullmodel_optimization:
		(full_model_params_side1,full_model_result_side1,full_model_optim_plots_side1) = CalcFullModel(load1,load2,InitialCoeffs_side1,Error_side1,npoints_side1,XPositions_side1,CTODValues_side1,InitialModels_side1,CrackCenterCoords,Symmetric_COD,side=1,minload=minload_side1,maxload=maxload_side1,seed_param=seed_param_side1_avg,nominal_length=nominal_length,nominal_modulus=nominal_modulus,nominal_stress=nominal_stress,doplots=False,fm_plotdata=fm_plotdata_side1_avg,opencl_ctx=ctx,opencl_dev=dev)

		
		(full_model_params_side2,full_model_result_side2,full_model_optim_plots_side2) = CalcFullModel(load1,load2,InitialCoeffs_side2,Error_side2,npoints_side2,XPositions_side2,CTODValues_side2,InitialModels_side2,CrackCenterCoords,Symmetric_COD,side=2,minload=minload_side2,maxload=maxload_side2,seed_param=seed_param_side2_avg,nominal_length=nominal_length,nominal_modulus=nominal_modulus,nominal_stress=nominal_stress,doplots=False,fm_plotdata=fm_plotdata_side2_avg,opencl_ctx=ctx,opencl_dev=dev)
		#Refines the spline fit to follow the monotoncially increaseing rule.
		(output_loads_full_1,tippos_side1_full,output_loads_full_2,tippos_side2_full) = calculate_closureprofile(0,maxload_side1,0,maxload_side2,50,full_model_params_side1,full_model_params_side2,TipCoords1[0],TipCoords2[0])

		#while tippos_side1_full[0]-crack_middle > 0 or  tippos_side2_full[0]-crack_middle < 0:
                        
                        #(output_loads_full_1,tippos_side1_full,output_loads_full_2,tippos_side2_full) = calculate_closureprofile(minload_side1,maxload_side1,minload_side2,maxload_side2,50,full_model_params_side1,full_model_params_side2,TipCoords1[0],TipCoords2[0])
			#minload_side1 = minload_side1 + 1e3
                        #minload_side2 = minload_side2 + 1e3
                        #pass
		
		#pl.plot(tippos_side2_full*1.0e6-crack_middle*1e6,output_loads_full_2/1.0e6,'k-', label="Integral Model Spline")
		
		
		pass
	#pl.axvline(x=TipCoords2[0]*1e6-crack_middle*1e6,label="Crack Tip Position")
	pl.xlabel(r'$X (\mu m)$')
	pl.ylabel(r'$Applied \ Load \ (MPa)$')
	#pl.xlim([0, (TipCoords2[0]-crack_middle+tip_tolerance)*1e6])
	#pl.xlim([-500, (3e-3)*1e6])
	#pl.legend(loc="upper left", prop={'size':16})
	pl.grid()
	pl.savefig('Sample_integral.png',dpi=600)

	

	
	#observed_reff = np.array([ 0.0,  1e-3, 1.5e-3, 2e-3  ],dtype='d')
	#observed_seff = np.array([ 10e6, 15e6, 30e6, 150e6  ],dtype='d')


        #Calculation of the closure stress distirubtions.  This is done using the refined spline fits of the closure points.  The left (side1) and right (side2) are examined seperatly
	sigma_closure_side1_avg = inverse_closure2(-(tippos_side1_avg[np.argwhere(tippos_side1_avg<crack_middle)]-crack_middle),output_loads1[np.argwhere(tippos_side1_avg<crack_middle)],x1,x_bnd_1,dx,a,sigma_yield,crack_model)
	
	sigma_closure_side2_avg = inverse_closure2(tippos_side2_avg[np.argwhere(tippos_side2_avg>crack_middle)]-crack_middle,output_loads2[np.argwhere(tippos_side2_avg>crack_middle)],x2,x_bnd_2,dx,a,sigma_yield,crack_model)
	#sigma_closure_side2_max = inverse_closure(tippos_side2_max-crack_middle,output_loads2,x2,x_bnd_2,dx,a,sigma_yield,crack_model)
	sigma_closure_side2_full = inverse_closure2(tippos_side2_full[np.argwhere(tippos_side2_full>crack_middle)]-crack_middle,output_loads_full_2[np.argwhere(tippos_side2_full>crack_middle)],x2,x_bnd_2,dx,a,sigma_yield,crack_model)	#Use this for FEM analysis
	sigma_closure_side1_full = inverse_closure2(-(tippos_side1_full[np.argwhere(tippos_side1_full<crack_middle)]-crack_middle),output_loads_full_1[np.argwhere(tippos_side1_full<crack_middle)],x1,x_bnd_1,dx,a,sigma_yield,crack_model)

	sigma_closure_side2_depth = inverse_closure2(tippos_side2_depth[np.argwhere(tippos_side2_depth>0.0)]-0.0,output_loads2_depth[np.argwhere(tippos_side2_depth>0.0)],x2,x_bnd_2,dx,a,sigma_yield,crack_model)

        #Plots of the left and right closure point
	pl.figure()
	pl.plot(tippos_side1_avg*1.0e6,output_loads1/1.0e6,'x-g')	
	#pl.plot(tippos_side1_full*1.0e3,output_loads_full_1/1.0e6,'o-g')
	pl.plot(tippos_side2_avg*1.0e6,output_loads2/1.0e6,'x-b')
	#pl.plot(tippos_side2_full*1.0e3,output_loads_full_2/1.0e6,'o-b')	
	pl.xlabel(r'$X (\mu m)$')
	pl.ylabel(r'$Applied \ Load \ (MPa)$')
	pl.legend(('YZ Plane Differential Model','YZ Plane Integral Model','XY Plane Differential Model','XY Plane Integral Model'),loc="best", prop={'size':14})
	pl.grid()
	pl.savefig('sample_flanks.png',dpi=600)
	
	if relshift_middleimg_lowerleft_corner_x_ref is not None:
		# Only applies to DIC dgs files generated through dc_process that have additional registration info added!
		TestRegistration(nloads,Xposvecs,u_disps,v_disps,
						 ROI_out_arrays,
						 relshift_middleimg_lowerleft_corner_x_ref=relshift_middleimg_lowerleft_corner_x_ref,
						 relshift_middleimg_lowerleft_corner_x_diff=relshift_middleimg_lowerleft_corner_x_diff,
						 relshift_middleimg_lowerleft_corner_y_ref=relshift_middleimg_lowerleft_corner_y_ref,
						 relshift_middleimg_lowerleft_corner_y_diff=relshift_middleimg_lowerleft_corner_y_diff)
		pass
	
	#Comparison with direct output from FEA with the calculation of the stress distribution using the DIC process
	#for 2D crack us -35
	#for 3D Ideal
	from load_params import load_params
	#Loads Para file which contains the apprioprate file paths.  
	load_params("/home/cgiuffre/crackclosure_abaqus_modeling_3d/crack_closure/crack_closure_params.py",globals())
	results_closure = json.load(open(results_file))
	stress_ypos_thermal_front = pd.read_csv(results_closure["csv_cpress_ypos_thermal"])
	stress_ypos_thermal_symm = pd.read_csv(results_closure["csv_cpress_ypos_thermal_symm"])
	pl.figure()

	#pl.plot(-(x1)*1e3,sigma_closure_side1_avg/1e6,label="Differential Model Inversion: Left Side",linewidth=2)	


	
	pl.plot(stress_ypos_thermal_front["x_position_mm"][0:-1],np.abs(stress_ypos_thermal_front["CPRESS_MPa"])[0:-1],'k-x',
	label="FEA Output: XY Plane",linewidth=2)
	pl.plot(-stress_ypos_thermal_symm["x_position_mm"][0:-1],np.abs(stress_ypos_thermal_symm["CPRESS_MPa"])[0:-1],'k-o',
	label="FEA Output: YZ Plane",linewidth=2)

	
	pl.plot((x2)[0:-1]*1e3,sigma_closure_side2_avg[0:-1]/1e6,'-',label="Differential Model Inversion",linewidth=2)	
	pl.plot((x2)[0:-1]*1e3,sigma_closure_side2_full[0:-1]/1e6,'-',label="Integral Model Inversion",linewidth=2)	
	pl.plot(x2*1e3,sigma_closure_side2_depth/1e6,'-',label="Depth  Inversion",linewidth=2)

	
	#ymaxaxis = max(-stress_ypos_thermal_front["CPRESS_MPa"])
	#pl.axis([0, 5.0, 0, 1000])
	pl.xlabel(r'$X (mm)$')
	pl.ylabel(r'$Closure \ Stress (MPa)$')
	pl.grid()
	pl.legend(loc="best", prop={'size':16})
	pl.savefig('closure_inverse.png',dpi=600)
        	
        #save_closurestress(profile_savefile_side1,x1,np.nan_to_num(sigma_closure_side1_full,copy=False),a)

        #save_closurestress(profile_savefile_side2,x2,np.nan_to_num(sigma_closure_side2_full,copy=False),a)

        #Saves the closure stress distirbution calculated using the full/integral model
        save_closurestress(profile_savefile_side1,x1,sigma_closure_side1_full,a)

        save_closurestress(profile_savefile_side2,x2,sigma_closure_side2_full,a)
	
	#Calculation of the crack depth behavior based on direct outputs from FEA simulation of a surface crack
	#Aspect Ratio Coefficient
        #Applied Load to the FEA simulation
	fea_load =np.array([0.4168250753978888, 25.58840607325236, 50.75975881958008, 75.93054493459066, 101.10072230529785, 126.27031204223633, 151.43932850646974, 176.60781278483074, 201.77576710001628, 226.943239654541, 252.109810353597, 277.2764172566732, 302.44256934611, 327.6083283589681, 352.773636311849, 377.9389294840495, 403.10422296142576, 428.2695157674154, 453.4348095296224, 478.6001022949219])*1e6
	#Aspect Ratio as outputted from the FEA simulation
	aspect_ratio_fea = np.array([0.33, 0.4545454446934474, 0.6428568813265542, 0.7500004656612804, 0.7777775864542216, 0.8421054744984661, 0.8500003814697266, 0.9047617641976567, 0.9090908893868948, 0.9090908893868948, 0.9130436630465396, 0.8800003051757812, 0.8800003051757812, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
	#Approximated depth and surface closure points based on the node of last contact
	depth_fea = np.array([0.4,1.0,1.8,2.5,2.8,3.2,3.2,3.8,4.0,4.0,4.2,4.2,4.2,5.0,5.0,5.0,5.0,5.0,5.0,5.0])
	front_fea = np.array([1.2,2.0,2.8,3.2,3.6,3.8,4.0,4.2,4.4,4.4,4.6,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0])
	
        #The three below sections calculate the intiail coefficient based on the differential, integral, and FEA outputted data.  Since the crack is known to never exceed an aspect ratio of 1 during the opening, only a portion of the model propsosed by Ravichandran (1997) is used and there is no need for both parts
	
	aoc_avg = (-tippos_side1_avg*1.0e6-crack_middle*1e6)/(tippos_side2_avg*1.0e6-crack_middle*1e6)
	aot = (-tippos_side1_avg*1.0e6-crack_middle*1e6)/(50*1e3)

	g = 1+(0.1+0.35*aot**2)
	M1 = 1.13-0.09*aoc_avg
	M2 = -0.54+(0.89/(0.2+aoc_avg))
	M3 = 0.5-(1/(0.65+aoc_avg))+14*(1-aoc_avg)**24
	Q = 1+1.464*aoc_avg**1.65
	F = (M1+M2*aot**2+M3*aot**4)*np.sqrt(np.pi/Q)
	fw = np.sqrt(np.cos((np.pi*(tippos_side2_avg*1.0e6-crack_middle*1e6))/(2*half_width*1e3)*np.sqrt(aot))**-1.0)
    
	xi = 1.6*np.sqrt(8)/np.sqrt(np.pi)*(1-0.3**2)*aoc_avg*F*fw*g
	#pl.figure(1001)
	#pl.plot(output_loads2*1e-6,xi/youngs_modulus*1e-6)
	#Integral Model
	aoc_full = (-tippos_side1_full*1.0e6-crack_middle*1e6)/(tippos_side2_full*1.0e6-crack_middle*1e6)
	aot = (-tippos_side1_full*1.0e6-crack_middle*1e6)/(50*1e3)

	g = 1+(0.1+0.35*aot**2)
	M1 = 1.13-0.09*aoc_full
	M2 = -0.54+(0.89/(0.2+aoc_full))
	M3 = 0.5-(1/(0.65+aoc_full))+14*(1-aoc_full)**24
	Q = 1+1.464*aoc_full**1.65
	F = (M1+M2*aot**2+M3*aot**4)*np.sqrt(np.pi/Q)
	fw = np.sqrt(np.cos((np.pi*(tippos_side2_full*1.0e6-crack_middle*1e6))/(2*half_width*1e3)*np.sqrt(aot))**-1.0)
    
	xi = 1.6*np.sqrt(8)/np.sqrt(np.pi)*(1-0.3**2)*aoc_full*F*fw*g
	#pl.figure(1001)
	#pl.plot(output_loads_full_2*1e-6,xi/youngs_modulus*1e-6)
	
	#FEA Output 
	aoc_fea = np.array(aspect_ratio_fea)
	aot = np.array(depth_fea)/50.0

	g = 1+(0.1+0.35*aot**2)
	M1 = 1.13-0.09*aoc_fea
	M2 = -0.54+(0.89/(0.2+aoc_fea))
	M3 = 0.5-(1/(0.65+aoc_fea))+14*(1-aoc_fea)**24
	Q = 1+1.464*aoc_fea**1.65
	F = (M1+M2*aot**2+M3*aot**4)*np.sqrt(np.pi/Q)
	fw = np.sqrt(np.cos((np.pi*(front_fea*1e3))/(2*half_width*1e3)*np.sqrt(aot))**-1.0)
    
	xi = 1.6*np.sqrt(8)/np.sqrt(np.pi)*(1-0.3**2)*aoc_fea*F*fw*g
	pl.figure(1001)
	pl.plot(xi,fea_load*1e-6)
	pl.legend([r'$\xi$',r'$FEA$'],loc='upper left')
	pl.xlabel(r'$\xi$')
	pl.ylabel(r'$Applied \ Load \ (MPa)$')
	pl.title('')
	#pl.grid()
	pl.savefig('fea_c5.png',dpi=600)
	
	pl.figure(1003)
	#pl.plot(output_loads1/1e6,aoc_avg)
	#pl.plot(output_loads_full_1/1e6,aoc_full)
	pl.plot(depth_fea*1e3,fea_load/1e6)
	pl.legend([r'$x_{d}$',r'$FEA$'],loc='upper left')   
	pl.xlabel(r'$x_{d} (\mu m)$')
	pl.ylabel(r'$Applied \ Load \ (MPa)$')
	pl.title('')
	#pl.grid()
	pl.savefig('fea_xd.png',dpi=600)
	
	pl.figure(1004)
	#pl.plot(output_loads1/1e6,aoc_avg)
	#pl.plot(output_loads_full_1/1e6,aoc_full)
	pl.plot(aspect_ratio_fea,fea_load/1e6)
	pl.legend([r'$x_{d}/x_{t}$',r'$FEA$'],loc='upper left')  
	pl.xlabel(r'$x_{d}/x_{t}$')
	pl.ylabel(r'$Applied \ Load \ (MPa)$')
	pl.title('')
	#pl.grid()
	pl.savefig('fea_car.png',dpi=600)	
	
	pl.figure(2001)
	#pl.plot(output_loads1/1e6,aoc_avg)
	#pl.plot(output_loads_full_1/1e6,aoc_full)
	pl.plot(front_fea*1e3,fea_load/1e6)
	pl.legend([r'$x_{t}$',r'$FEA$'],loc='upper left')  
	pl.xlabel(r'$x (\mu m)$')
	pl.ylabel(r'$Applied \ Load \ (MPa)$')
	pl.title('')
	pl.grid()
	pl.savefig('fea_xt.png',dpi=600)
	
	## Plot diagnostics...
	# Should have at least one plot that evaluates
	# overall performance in fitting the entire data set.
	#
	# ... how to collapse it down to 2D?
	# (Could add additional lines to the InitialModel plots...)
	#pl.figure()
	
	"""
	#This is for 013L
	xaxis = (np.arange(0.0,(225-25)*dic_dx,dic_dx)-crack_middle)*1e6
	yaxis = (np.arange(0.0,(260-150)*dic_dy,dic_dy)-crack_middle)*1e6
	
	idx1_plot = 5
	idx2_plot = 7
	#Offset by -11.25 micron
	
	sigma_ave = (load2[idx1_plot,idx2_plot,0]/1.e6+load1[idx1_plot,idx2_plot,0]/1.e6)/2
	sigma_diff = (load2[idx1_plot,idx2_plot,0]/1.e6-load1[idx1_plot,idx2_plot,0]/1.e6)
	
        #The contourf figures plot the direct DIC results for the first attempt at gathering data form C18-AFVT-013L with the file call being in the large block of short titanium cracks.  This is formatted such that it will only plot the region of the image where DIC displacement fields are present.  Does not work for other applied load increments and other data sets.  The line plot calls work for any applied load increment and for any sample.  It is important to have the sigma_ave and sigma_diff claculations anytime idx1_plot and idx2_plot are changed  The axis are offset such that the origin is around the crack center.  This is currently configured for the above stated sample  but also works for the other C18-AFVT-013L data since the crack is in the same position with the same center measurement  All of the generated plots save at 600 dpi.  
	pl.figure()
	pl.contourf(xaxis+600,yaxis,np.transpose(v_disps[25:225,150:260,idx1_plot,idx2_plot,4])*1e6+11.25)
	pl.xlabel(r'$X (\mu m)$')
	pl.ylabel(r'$Y (\mu m)$')
	pl.title(r'$\sigma_{ave} = %i;\ \Delta \sigma = %i$' % (sigma_ave,sigma_diff))
	cbar = pl.colorbar()
	cbar.set_label(r'$(\mu m)$', rotation=0)
	pl.grid()
	pl.savefig('mostly_closed.png',dpi=600)
	
	
	XPositionsSort_side2=np.argsort(XPositions_side2[idx1_plot,idx2_plot])
	XPositionsSorted_side2=XPositions_side2[idx1_plot,idx2_plot][XPositionsSort_side2]
	InitialModelValuesSorted_side2=InitialModels_side2[idx1_plot,idx2_plot][XPositionsSort_side2]
	pl.figure()
	pl.plot((XPositions_side2[idx1_plot,idx2_plot]-np.min(XPositions_side2[idx1_plot,idx2_plot]))*1e6,CTODValues_side2[idx1_plot,idx2_plot]*1e6,'k.',label="Corrected COD")
	pl.plot((XPositionsSorted_side2-np.min(XPositions_side2[idx1_plot,idx2_plot]))*1e6,InitialModelValuesSorted_side2*1e6,'r-',label="Initial Fit")
	pl.xlabel(r'$X  (\mu m)$')
	pl.ylabel(r'$COD  (\mu m)$')
	pl.title(r'$\sigma_{ave} = %i \ and \ \sigma_{diff} = %i$' % (sigma_ave,sigma_diff))
	pl.grid()
	pl.savefig('lowdiff_model.png',dpi=600)
	
	idx1_plot = 15
	idx2_plot = 17
	#Offset by -8.85 micron
	
	sigma_ave = (load2[idx1_plot,idx2_plot,0]/1.e6+load1[idx1_plot,idx2_plot,0]/1.e6)/2
	sigma_diff = (load2[idx1_plot,idx2_plot,0]/1.e6-load1[idx1_plot,idx2_plot,0]/1.e6)
	
	pl.figure()
	pl.contourf(xaxis+600,yaxis+900,np.transpose(v_disps[25:225,150:260,idx1_plot,idx2_plot,4])*1e6+8.85)
	pl.title(r'$\sigma_{ave} = %i;\ \Delta \sigma = %i$' % (sigma_ave,sigma_diff))
	pl.xlabel(r'$X (\mu m)$')
	pl.ylabel(r'$Y (\mu m)$')
	cbar = pl.colorbar()
	cbar.set_label(r'$(\mu m)$', rotation=0)
	pl.savefig('mostly_open.png',dpi=600)
	pl.grid()
	
	XPositionsSort_side2=np.argsort(XPositions_side2[idx1_plot,idx2_plot])
	XPositionsSorted_side2=XPositions_side2[idx1_plot,idx2_plot][XPositionsSort_side2]
	InitialModelValuesSorted_side2=InitialModels_side2[idx1_plot,idx2_plot][XPositionsSort_side2]
	pl.figure()
	#pl.plot(XPositions_side1[idx1_plot,idx2_plot]*1e3,CTODValues_side1[idx1_plot,idx2_plot]*1e6,'k.',
			#XPositionsSorted_side1*1e3,InitialModelValuesSorted_side1*1e6,'r-')
	pl.plot((XPositions_side2[idx1_plot,idx2_plot]-np.min(XPositions_side2[idx1_plot,idx2_plot]))*1e6,CTODValues_side2[idx1_plot,idx2_plot]*1e6,'k.')
	pl.xlabel(r'$X (\mu m)$')
	pl.ylabel(r'$COD  (\mu m)$')
	pl.title(r'$\sigma_{ave} = %i;\ \Delta \sigma = %i$' % (sigma_ave,sigma_diff))
	pl.grid()
	pl.savefig('mostly_open_profile.png',dpi=600)
	
	XPositionsSort_side2=np.argsort(XPositions_side2[idx1_plot,idx2_plot])
	XPositionsSorted_side2=XPositions_side2[idx1_plot,idx2_plot][XPositionsSort_side2]
	InitialModelValuesSorted_side2=InitialModels_side2[idx1_plot,idx2_plot][XPositionsSort_side2]
	
	pl.figure()
	pl.plot((XPositions_side2[idx1_plot,idx2_plot]-np.min(XPositions_side2[idx1_plot,idx2_plot]))*1e6,CTODValues_side2[idx1_plot,idx2_plot]*1e6,'k.',label="Corrected COD")
	pl.plot((XPositionsSorted_side2-np.min(XPositions_side2[idx1_plot,idx2_plot]))*1e6,InitialModelValuesSorted_side2*1e6,'r-',label="Initial Fit")
	pl.xlabel(r'$X (\mu m)$')
	pl.ylabel(r'$COD (\mu m)$')
	pl.grid()
	pl.title(r'$\sigma_{ave} = %i;\ \Delta \sigma = %i$' % (sigma_ave,sigma_diff))
	pl.legend(loc="best")
	pl.savefig('mostly_open_fitt.png',dpi=600)
	

	
	idx1_plot = 11
	idx2_plot = 21
	#Offset by -40 micron
	sigma_ave = (load2[idx1_plot,idx2_plot,0]/1.e6+load1[idx1_plot,idx2_plot,0]/1.e6)/2
	sigma_diff = (load2[idx1_plot,idx2_plot,0]/1.e6-load1[idx1_plot,idx2_plot,0]/1.e6)
	
	pl.figure()
	pl.contourf(xaxis+600,yaxis+900,np.transpose(v_disps[25:225,150:260,idx1_plot,idx2_plot,4])*1e6+40)
	pl.xlabel(r'$X \ Position (\mu m)$')
	pl.ylabel(r'$Y \ Position (\mu m)$')
	pl.title(r'$\sigma_{ave} = %i;\ \Delta \sigma = %i$' % (sigma_ave,sigma_diff))
	cbar = pl.colorbar()
	cbar.set_label(r'$(\mu m)$', rotation=0)
	pl.grid()
	pl.savefig('highdiff_field.png',dpi=600)
	
	pl.figure()
	pl.plot(yaxis+900,v_disps[175,150:260,idx1_plot,idx2_plot,4]*1e6+40.15,'k+',label="Line A")	
	pl.plot(yaxis+900,v_disps[210,150:260,idx1_plot,idx2_plot,4]*1e6+40.15,'kx',label="Line B")
	pl.plot(yaxis+900,v_disps[220,150:260,idx1_plot,idx2_plot,4]*1e6+40.15,'k.',label="Line C")
	pl.plot(yaxis+900,(sigma_diff/YoungsModulus)*(yaxis+900)*1e6,label="COD Correction")
	pl.legend(loc="lower right")
	pl.xlabel(r'$Y  (\mu m)$')
	pl.ylabel(r'$Vertical \ Displacement (\mu m)$')
	pl.title(r'$\sigma_{ave} = %i;\ \Delta \sigma = %i$' % (sigma_ave,sigma_diff))
	pl.ylim([-1.75, 2])
	pl.xlim([-250, 400])
	pl.grid()
	pl.savefig('highdiff_columns.png',dpi=600)	
	
	XPositionsSort_side2=np.argsort(XPositions_side2[idx1_plot,idx2_plot])
	XPositionsSorted_side2=XPositions_side2[idx1_plot,idx2_plot][XPositionsSort_side2]
	InitialModelValuesSorted_side2=InitialModels_side2[idx1_plot,idx2_plot][XPositionsSort_side2]
	
	pl.figure()
	pl.plot((XPositions_side2[idx1_plot,idx2_plot]-np.min(XPositions_side2[idx1_plot,idx2_plot]))*1e6,CTODValues_side2[idx1_plot,idx2_plot]*1e6,'k.',label="Corrected COD")
	pl.plot((XPositionsSorted_side2-np.min(XPositions_side2[idx1_plot,idx2_plot]))*1e6,InitialModelValuesSorted_side2*1e6,'r-',label="Initial Fit")
	pl.xlabel(r'$X (\mu m)$')
	pl.ylabel(r'$COD (\mu m)$')
	pl.grid()
	pl.title(r'$\sigma_{ave} = %i;\ \Delta \sigma = %i$' % (sigma_ave,sigma_diff))
	pl.legend(loc="best")
	pl.savefig('highdiff_model.png',dpi=600)
	
	pl.figure()
	pl.plot((XPositions_side2[idx1_plot,idx2_plot]-np.min(XPositions_side2[idx1_plot,idx2_plot]))*1e6,CTODValues_side2[idx1_plot,idx2_plot]*1e6,'k.',label="Corrected COD")
	pl.plot((XPositions_side2[idx1_plot,idx2_plot]-np.min(XPositions_side2[idx1_plot,idx2_plot]))*1e6,(CTODValues_side2[idx1_plot,idx2_plot])*1e6+(sigma_diff*1e6/YoungsModulus)*(dic_dy*dic_span)*1e6,'r.',label="Raw Data")
	pl.xlabel(r'$X  (\mu m)$')
	pl.ylabel(r'$COD (\mu m)$')
	pl.grid()
	pl.title(r'$\sigma_{ave} = %i;\ \Delta \sigma = %i$' % (sigma_ave,sigma_diff))
	pl.legend(loc="best")
	pl.savefig('strain_correction.png',dpi=600)
	
	idx1_plot = 10
	idx2_plot = 20
	#Offset by -40 micron
	sigma_ave = (load2[idx1_plot,idx2_plot,0]/1.e6+load1[idx1_plot,idx2_plot,0]/1.e6)/2
	sigma_diff = (load2[idx1_plot,idx2_plot,0]/1.e6-load1[idx1_plot,idx2_plot,0]/1.e6)
	
	XPositionsSort_side2=np.argsort(XPositions_side2[idx1_plot,idx2_plot])
	XPositionsSorted_side2=XPositions_side2[idx1_plot,idx2_plot][XPositionsSort_side2]
	InitialModelValuesSorted_side2=InitialModels_side2[idx1_plot,idx2_plot][XPositionsSort_side2]
	
	pl.figure()
	pl.plot((XPositions_side2[idx1_plot,idx2_plot]-np.min(XPositions_side2[idx1_plot,idx2_plot]))*1e6,CTODValues_side2[idx1_plot,idx2_plot]*1e6,'k.',label="Corrected COD")
	pl.plot((XPositionsSorted_side2-np.min(XPositions_side2[idx1_plot,idx2_plot]))*1e6,InitialModelValuesSorted_side2*1e6,'r-',label="Differential Fit")
	pl.xlabel(r'$X (\mu m)$')
	pl.ylabel(r'$COD (\mu m)$')
	pl.grid()
	pl.title(r'$\sigma_{ave} = %i;\ \Delta \sigma = %i$' % (sigma_ave,sigma_diff))
	pl.legend(loc="best")
	pl.savefig('highdiff_model.png',dpi=600)
	
	idx1_plot = 12
	idx2_plot = 18
	
	#Offset by -40 micron
	sigma_ave = (load2[idx1_plot,idx2_plot,0]/1.e6+load1[idx1_plot,idx2_plot,0]/1.e6)/2
	sigma_diff = (load2[idx1_plot,idx2_plot,0]/1.e6-load1[idx1_plot,idx2_plot,0]/1.e6)
	
	XPositionsSort_side2=np.argsort(XPositions_side2[idx1_plot,idx2_plot])
	XPositionsSorted_side2=XPositions_side2[idx1_plot,idx2_plot][XPositionsSort_side2]
	InitialModelValuesSorted_side2=InitialModels_side2[idx1_plot,idx2_plot][XPositionsSort_side2]
	
	pl.figure()
	pl.plot((XPositions_side2[idx1_plot,idx2_plot]-np.min(XPositions_side2[idx1_plot,idx2_plot]))*1e6,CTODValues_side2[idx1_plot,idx2_plot]*1e6,'k.',label="Corrected COD")
	pl.plot((XPositionsSorted_side2-np.min(XPositions_side2[idx1_plot,idx2_plot]))*1e6,InitialModelValuesSorted_side2*1e6,'r-',label="Differential Fit")
	pl.xlabel(r'$X (\mu m)$')
	pl.ylabel(r'$COD (\mu m)$')
	pl.grid()
	pl.title(r'$\sigma_{ave} = %i;\ \Delta \sigma = %i$' % (sigma_ave,sigma_diff))
	pl.legend(loc="best")
	pl.savefig('meddiff_model.png',dpi=600)	
	
	idx1_plot = 14
	idx2_plot = 16
	#Offset by -40 micron
	sigma_ave = (load2[idx1_plot,idx2_plot,0]/1.e6+load1[idx1_plot,idx2_plot,0]/1.e6)/2
	sigma_diff = (load2[idx1_plot,idx2_plot,0]/1.e6-load1[idx1_plot,idx2_plot,0]/1.e6)
	
	XPositionsSort_side2=np.argsort(XPositions_side2[idx1_plot,idx2_plot])
	XPositionsSorted_side2=XPositions_side2[idx1_plot,idx2_plot][XPositionsSort_side2]
	InitialModelValuesSorted_side2=InitialModels_side2[idx1_plot,idx2_plot][XPositionsSort_side2]
	
	pl.figure()
	pl.plot((XPositions_side2[idx1_plot,idx2_plot]-np.min(XPositions_side2[idx1_plot,idx2_plot]))*1e6,CTODValues_side2[idx1_plot,idx2_plot]*1e6,'k.',label="Corrected COD")
	pl.plot((XPositionsSorted_side2-np.min(XPositions_side2[idx1_plot,idx2_plot]))*1e6,InitialModelValuesSorted_side2*1e6,'r-',label="Differential Fit")
	pl.xlabel(r'$X (\mu m)$')
	pl.ylabel(r'$COD (\mu m)$')
	pl.grid()
	pl.title(r'$\sigma_{ave} = %i;\ \Delta \sigma = %i$' % (sigma_ave,sigma_diff))
	pl.legend(loc="best")
	pl.savefig('lowdiff_model.png',dpi=600)	
    
    
	idx1_plot = 15
	idx2_plot = 17
	#Offset by -40 micron
	sigma_ave = (load2[idx1_plot,idx2_plot,0]/1.e6+load1[idx1_plot,idx2_plot,0]/1.e6)/2
	sigma_diff = (load2[idx1_plot,idx2_plot,0]/1.e6-load1[idx1_plot,idx2_plot,0]/1.e6)
	XPositionsSort_side2=np.argsort(XPositions_side2[idx1_plot,idx2_plot])
	XPositionsSorted_side2=XPositions_side2[idx1_plot,idx2_plot][XPositionsSort_side2]
	InitialModelValuesSorted_side2=InitialModels_side2[idx1_plot,idx2_plot][XPositionsSort_side2]
	
	pl.figure()
	pl.plot((XPositions_side2[idx1_plot,idx2_plot]-np.min(XPositions_side2[idx1_plot,idx2_plot]))*1e6,CTODValues_side2[idx1_plot,idx2_plot]*1e6,'k.',label="Corrected COD")
	pl.plot((XPositionsSorted_side2-np.min(XPositions_side2[idx1_plot,idx2_plot]))*1e6,InitialModelValuesSorted_side2*1e6,'r-',label="Differential Fit")

	idx1_plot = 11
	idx2_plot = 21
	#Offset by -40 micron
	sigma_ave = (load2[idx1_plot,idx2_plot,0]/1.e6+load1[idx1_plot,idx2_plot,0]/1.e6)/2
	sigma_diff = (load2[idx1_plot,idx2_plot,0]/1.e6-load1[idx1_plot,idx2_plot,0]/1.e6)
	XPositionsSort_side2=np.argsort(XPositions_side2[idx1_plot,idx2_plot])
	XPositionsSorted_side2=XPositions_side2[idx1_plot,idx2_plot][XPositionsSort_side2]
	InitialModelValuesSorted_side2=InitialModels_side2[idx1_plot,idx2_plot][XPositionsSort_side2]
	

	pl.plot((XPositions_side2[idx1_plot,idx2_plot]-np.min(XPositions_side2[idx1_plot,idx2_plot]))*1e6,CTODValues_side2[idx1_plot,idx2_plot]*1e6,'k.',label="Corrected COD")
	pl.plot((XPositionsSorted_side2-np.min(XPositions_side2[idx1_plot,idx2_plot]))*1e6,InitialModelValuesSorted_side2*1e6,'r-',label="Differential Fit")
    
    """
   
	
	#Plotting specific to FEA dat with a new contourf plot and new opening plots.  This is formatted for the 2D through crack case since this plots both sides of the crack which are not present in the 3D simulations due to symmetry.  However, it would be interesting to see the difference between the displacement fields for the surface crack case.
	"""
	xaxis = (np.arange(0.0,(151-0)*dic_dx,dic_dx)-crack_middle)*1e6
	yaxis = (np.arange(0.0,(151-0)*dic_dy,dic_dy)-crack_middle)*1e6

	idx1_plot=5
	idx2_plot=7
	
	sigma_ave = (load2[idx1_plot,idx2_plot,0]/1.e6+load1[idx1_plot,idx2_plot,0]/1.e6)/2
	sigma_diff = (load2[idx1_plot,idx2_plot,0]/1.e6-load1[idx1_plot,idx2_plot,0]/1.e6)	
	
	pl.figure()
	pl.contourf(xaxis,yaxis,np.transpose(v_disps[:,:,idx1_plot,idx2_plot,0])*1e6-73)	
	pl.xlabel(r'$X (\mu m)$')
	pl.ylabel(r'$Y (\mu m)$')
	pl.title(r'$\sigma_{ave} = %i;\ \Delta \sigma = %i$' % (sigma_ave,sigma_diff))
	cbar = pl.colorbar()
	cbar.set_label(r'$(\mu m)$', rotation=0)
	pl.grid()
	pl.savefig('dic_field.png',dpi=600)
	XPositionsSort_side2=np.argsort(XPositions_side2[idx1_plot,idx2_plot])
	XPositionsSorted_side2=XPositions_side2[idx1_plot,idx2_plot][XPositionsSort_side2]
	InitialModelValuesSorted_side2=InitialModels_side2[idx1_plot,idx2_plot][XPositionsSort_side2]
	pl.savefig('fea_contour.png',dpi=600)	
	
	pl.figure()
	pl.plot((XPositions_side2[idx1_plot,idx2_plot]-np.min(XPositions_side2[idx1_plot,idx2_plot]))*1e6,CTODValues_side2[idx1_plot,idx2_plot]*1e6,'k.',label="Corrected COD")
	pl.plot((XPositionsSorted_side2-np.min(XPositions_side2[idx1_plot,idx2_plot]))*1e6,InitialModelValuesSorted_side2*1e6,'r-',label="Initial Fit")
	pl.xlabel(r'$X (\mu m)$')
	pl.ylabel(r'$COD (\mu m)$')
	pl.grid()
	pl.title(r'$\sigma_{ave} = %i \ and \ \Delta \sigma = %i$' % (sigma_ave,sigma_diff))
	pl.legend(loc="best")
	"""
	
	idx1_plot=5
	idx2_plot=7
	
	sigma_ave = (load2[idx1_plot,idx2_plot,0]/1.e6+load1[idx1_plot,idx2_plot,0]/1.e6)/2
	sigma_diff = (load2[idx1_plot,idx2_plot,0]/1.e6-load1[idx1_plot,idx2_plot,0]/1.e6)	
	
	XPositionsSort_side2=np.argsort(XPositions_side2[idx1_plot,idx2_plot])
	XPositionsSorted_side2=XPositions_side2[idx1_plot,idx2_plot][XPositionsSort_side2]
	InitialModelValuesSorted_side2=InitialModels_side2[idx1_plot,idx2_plot][XPositionsSort_side2]
	
	pl.figure()
	pl.plot((XPositions_side2[idx1_plot,idx2_plot]-np.min(XPositions_side2[idx1_plot,idx2_plot]))*1e6,CTODValues_side2[idx1_plot,idx2_plot]*1e6,'k.')
	pl.plot((XPositionsSorted_side2-np.min(XPositions_side2[idx1_plot,idx2_plot]))*1e6,InitialModelValuesSorted_side2*1e6,'r-',label=r'$\Delta \sigma = 50 MPa')
	pl.xlabel(r'$X (\mu m)$')
	pl.ylabel(r'$COD (\mu m)$')
	pl.grid()
	#pl.title(r'$\sigma_{ave} = %i \ and \ \Delta \sigma = %i$' % (sigma_ave,sigma_diff))
	#pl.legend(loc="best")
	
	idx1_plot=1
	idx2_plot=15
	
	sigma_ave = (load2[idx1_plot,idx2_plot,0]/1.e6+load1[idx1_plot,idx2_plot,0]/1.e6)/2
	sigma_diff = (load2[idx1_plot,idx2_plot,0]/1.e6-load1[idx1_plot,idx2_plot,0]/1.e6)	
	
	XPositionsSort_side2=np.argsort(XPositions_side2[idx1_plot,idx2_plot])
	XPositionsSorted_side2=XPositions_side2[idx1_plot,idx2_plot][XPositionsSort_side2]
	InitialModelValuesSorted_side2=InitialModels_side2[idx1_plot,idx2_plot][XPositionsSort_side2]
	

	pl.plot((XPositions_side2[idx1_plot,idx2_plot]-np.min(XPositions_side2[idx1_plot,idx2_plot]))*1e6,CTODValues_side2[idx1_plot,idx2_plot]*1e6,'k.')
	pl.plot((XPositionsSorted_side2-np.min(XPositions_side2[idx1_plot,idx2_plot]))*1e6,InitialModelValuesSorted_side2*1e6,'b-',label=r'$\Delta \sigma = 350 MPa')
	pl.legend(loc="best")
	
	pl.show()	

	pass
