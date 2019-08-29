import sys
import ast
import copy
import posixpath
import multiprocessing
import numpy as np

from matplotlib import pyplot as pl

import limatix.timestamp
from closure_measurements import process_dic

from closure_measurements.process_dic import load_dgs,Calc_CTODs,CalcInitialModel,EvalEffectiveTip,InitializeFullModel

from limatix import dc_value
from limatix import xmldoc
from limatix.dc_value import numericunitsvalue as numericunitsv
from limatix.dc_value import hrefvalue as hrefv
from limatix.dc_value import xmltreevalue as xmltreev

def run(_xmldoc,_element,
        _dest_href,
        _inputfilename,
        dc_measnum_int,
        dc_scan_outdic_href,
        dc_coordinatetransform,
        dc_spcYoungsModulus_numericunits,
        dc_specimen_str,
        dc_dic_span_int=20,
        dc_dic_smoothing_window_int=3,
        # NOTE: Default tip tolerance is smaller here than in command line version
        # because typically this is used with data that has gone through a registration
        # process and therefore has tighter tolerances
        dc_dic_tip_tolerance_numericunits=numericunitsv(50e-6,"m"),
        dc_min_dic_points_per_meter_numericunits=numericunitsv(40000,"1/m"),
        dc_symmetric_cod_bool=True,
        debug_bool=False):
    
    # Non-adjustable parameters
    nominal_length=2e-3 # nominal crack length, for nondimensional normalization
    #nominal_modulus=100.0e9 # nominal modulus
    nominal_stress=50e6 # nominal stress
    
    num_output_loads = 15
    
    # Coordinates generated by
    # this processtrak step shall be in stitched image
    # coordinates (see ExtractAndSitchOpticalImageJ.py _aligned_measnum.dgs and RegisterOpticalData.py)

    min_dic_points_per_meter=dc_min_dic_points_per_meter_numericunits.value("1/m")

    # To get fiducial coordinates, add xshift and yshift.
    
    YoungsModulus = dc_spcYoungsModulus_numericunits.value("Pa")
    
    xshift = numericunitsv.fromxml(_xmldoc, _xmldoc.xpathsinglecontext(dc_coordinatetransform, 'dc:translation/dc:xtranslation')).value('m')
    yshift = numericunitsv.fromxml(_xmldoc, _xmldoc.xpathsinglecontext(dc_coordinatetransform, 'dc:translation/dc:ytranslation')).value('m')
    
    #keypoints = _opxmldoc.xpathcontext(dc_crackpath, 'dc:segment/dc:keypoint')

    tip_tolerance = dc_dic_tip_tolerance_numericunits.value("m")
    
    (dic_dx,dic_dy,
     dic_nx,dic_ny,
     XRangeSize,
     nloads,
     Xinivec,Xposvecs,
     load1,load2,u_disps,v_disps,
     ROI_out_arrays,
     CrackCenterX,TipCoords1,TipCoords2,
     ROI_dic_yminidx,ROI_dic_ymaxidx,
     relshift_middleimg_lowerleft_corner_x_ref,
     relshift_middleimg_lowerleft_corner_x_diff,
     relshift_middleimg_lowerleft_corner_y_ref,
     relshift_middleimg_lowerleft_corner_y_diff) = load_dgs(dc_scan_outdic_href.getpath())



    CTODs = Calc_CTODs(dic_nx,nloads,XRangeSize,Xposvecs,v_disps,ROI_out_arrays,ROI_dic_yminidx,ROI_dic_ymaxidx,dc_dic_span_int,dc_dic_smoothing_window_int)

    (InitialModels_side1,
     InitialCoeffs_side1,
     Error_side1,
     npoints_side1,
     XPositions_side1,
     CTODValues_side1) = CalcInitialModel(nloads,CTODs,
                                          load1,load2,
                                          Xposvecs,CrackCenterX,
                                          dic_dy,dc_dic_span_int,
                                          dc_symmetric_cod_bool,1,YoungsModulus,
                                          relshift_middleimg_lowerleft_corner_x_ref=relshift_middleimg_lowerleft_corner_x_ref,
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
                                          dic_dy,dc_dic_span_int,
                                          dc_symmetric_cod_bool,2,YoungsModulus,
                                          relshift_middleimg_lowerleft_corner_x_ref=relshift_middleimg_lowerleft_corner_x_ref,
                                          nominal_length=nominal_length,nominal_stress=nominal_stress,
                                          doplots=False)


    (minload_side1,maxload_side1,seed_param_side1,fm_plots_side1) = InitializeFullModel(load1,load2,TipCoords1,TipCoords2,InitialCoeffs_side1,Error_side1,npoints_side1,XPositions_side1,CTODValues_side1,InitialModels_side1,CrackCenterX,tip_tolerance,min_dic_points_per_meter,dc_symmetric_cod_bool,side=1,doplots=True)
    (fitplot_side1,pickableplot_side1,c5plot_side1)=fm_plots_side1

    (minload_side2,maxload_side2,seed_param_side2,fm_plots_side2) = InitializeFullModel(load1,load2,TipCoords1,TipCoords2,InitialCoeffs_side2,Error_side2,npoints_side2,XPositions_side2,CTODValues_side2,InitialModels_side2,CrackCenterX,tip_tolerance,min_dic_points_per_meter,dc_symmetric_cod_bool,side=2,doplots=True)
    (fitplot_side2,pickableplot_side2,c5plot_side2)=fm_plots_side2


    (output_loads,tippos_side1,tippos_side2) =  process_dic.calculate_closureprofile(load1,num_output_loads,seed_param_side1,seed_param_side2,TipCoords1,TipCoords2)
    
    #closureprofile_side1 = xmldoc.xmldoc.newdoc("dc:closureprofile",nsmap={"dc":"http://limatix.org/datacollect"},contexthref=_dest_href)
    #for loadcnt in range(num_output_loads):
    #    newel = closureprofile_side1.addsimpleelement(closureprofile_side1.getroot(),"dc:tippos",(tippos_side1[loadcnt],"m"))
    #    closureprofile_side1.setattr(newel,"load_pascals",str(output_loads[loadcnt]))
    #    pass
    #closureprofile_side1_tree = xmltreev(closureprofile_side1)
    #    
    #closureprofile_side2 = xmldoc.xmldoc.newdoc("dc:closureprofile",nsmap={"dc":"http://limatix.org/datacollect"},contexthref=_dest_href)
    #for loadcnt in range(num_output_loads):
    #    newel = closureprofile_side2.addsimpleelement(closureprofile_side2.getroot(),"dc:tippos",(tippos_side2[loadcnt],"m"))
    #    closureprofile_side2.setattr(newel,"load_pascals",str(output_loads[loadcnt]))
    #    pass
    #closureprofile_side2_tree = xmltreev(closureprofile_side2)

    outdic_basename = posixpath.splitext(dc_scan_outdic_href.get_bare_quoted_filename())[0]
    if posixpath.splitext(outdic_basename)[1]==".dgs":  # Would happen if what we just split off was a .bz2, .gz, etc.
        outdic_basename = posixpath.splitext(outdic_basename)[0]
        pass

    closureprofile_href = hrefv(outdic_basename+"_closureprofile.csv",contexthref=_dest_href)

    process_dic.save_closureprofile(closureprofile_href.getpath(),output_loads,tippos_side1,tippos_side2)

    
    fitplot_side1_href = hrefv(outdic_basename+"_tipfit_side1.png",contexthref=_dest_href)
    pl.figure(fitplot_side1.number)
    pl.savefig(fitplot_side1_href.getpath(),dpi=300,transparent=True)

    fitplot_side2_href = hrefv(outdic_basename+"_tipfit_side2.png",contexthref=_dest_href)
    pl.figure(fitplot_side2.number)
    pl.savefig(fitplot_side2_href.getpath(),dpi=300,transparent=True)

    closureprofile_plot_href = hrefv(outdic_basename+"_closureprofile.png",contexthref=_dest_href)
    pl.figure()
    pl.plot(tippos_side1*1e3,output_loads/1e6,'-',
            tippos_side2*1e3,output_loads/1e6,'-')
    pl.grid()
    pl.xlabel('Tip position (relative to stitched image, mm)')
    pl.ylabel('Applied load (MPa)')
    pl.legend(('Side 1','Side 2'),loc="best")
    pl.title(dc_specimen_str)
    pl.savefig(closureprofile_plot_href.getpath(),dpi=300,transparent=True)

                         
    pl.close('all')
    return [
        #(("dc:closureprofile",{"side": "1"}),closureprofile_side1_tree),
        #(("dc:closureprofile",{"side": "2"}),closureprofile_side2_tree),
        ("dc:closureprofile",closureprofile_href),
        ("dc:closureprofileplot",closureprofile_plot_href),
        (("dc:dic_tip_fit",{"side": "1"}),fitplot_side1_href),
        (("dc:dic_tip_fit",{"side": "2"}),fitplot_side2_href),
    ]
