import sys
import ast
import posixpath
import multiprocessing
import numpy as np

import limatix.timestamp
from closure_measurements import perform_dic

from closure_measurements.perform_dic import execute_dic_loaded_data

from limatix import dc_value
from limatix.dc_value import numericunitsvalue as numericunitsv
from limatix.dc_value import hrefvalue as hrefv

def run(_xmldoc,_element,
        _dest_href,
        _inputfilename,
        dc_measnum_int,
        #dc_scan_outdgd_str,
        dc_scan_outdgd_href,
        dc_crackpath,
        dc_coordinatetransform,
        dic_scalefactor_int=5,
        dic_radius_int=20,
        dic_span_int=20,
        n_threads=multiprocessing.cpu_count(),
        debug_bool=False):

    # Coordinates in _dic.dgs files genereted by
    # this processtrak step shall be in stitched optical image
    # coordinates (see ExtractAndSitchOpticalImageJ.py _aligned_measnum.dgs and RegisterOpticalData.py)

    # To generate fiducial coordinates from coordinates in the _dic.dgs files:
    #   * add (xshift,yshift)
    

    xshift = numericunitsv.fromxml(_xmldoc, _xmldoc.xpathsinglecontext(dc_coordinatetransform, 'dc:translation/dc:xtranslation')).value('m')
    yshift = numericunitsv.fromxml(_xmldoc, _xmldoc.xpathsinglecontext(dc_coordinatetransform, 'dc:translation/dc:ytranslation')).value('m')
    
    #keypoints = _opxmldoc.xpathcontext(dc_crackpath, 'dc:segment/dc:keypoint')

    crackstartx = numericunitsv.fromxml(_xmldoc,_xmldoc.xpathsinglecontext(dc_crackpath, 'dc:segment/dc:keypoint[1]/dc:xcoordinate')).value('m') - xshift
    crackstarty = numericunitsv.fromxml(_xmldoc,_xmldoc.xpathsinglecontext(dc_crackpath, 'dc:segment/dc:keypoint[1]/dc:ycoordinate')).value('m') - yshift

    (crackally,crackally_units) = _xmldoc.xpathcontextnumpy(dc_crackpath, 'dc:segment/dc:keypoint/dc:ycoordinate')
    assert(limatix.lm_units.compareunits(crackally_units,limatix.lm_units.parseunits('m'))==1.0) # units of crack coordinates must be meters
    crackmaxy = np.max(crackally) - yshift
    crackminy = np.min(crackally) - yshift
    
    crackendx = numericunitsv.fromxml(_xmldoc,_xmldoc.xpathsinglecontext(dc_crackpath, 'dc:segment/dc:keypoint[last()]/dc:xcoordinate')).value('m') - xshift
    crackendy = numericunitsv.fromxml(_xmldoc,_xmldoc.xpathsinglecontext(dc_crackpath, 'dc:segment/dc:keypoint[last()]/dc:ycoordinate')).value('m') - yshift

    assert(crackstartx < crackendx)
    TipCoords1=(crackstartx,crackstarty)
    TipCoords2=(crackendx,crackendy)

    
    #dc_scan_outdgd_href = hrefv(dc_scan_outdgd_str,contexthref=_dest_href)

    (Images,x0,y0,dx,dy,nx,ny,nimages,nloads,ybase,YMotionPosns,StressPosns,ActualStressPosns,LowerLeft_XCoordinates_MotionController) = perform_dic.load_dgd(dc_scan_outdgd_href.getpath())    

    # LowerLeft_XCoordinates_MotionController are the
    # X coordinates of the lower left of each image, as
    # measured by the motion controller.
    # The middle entry is defined to be 0 (reference motion position)
    # and values decrease for increasing motion controller "Y" position.

    # In this case, since we have the registration data from
    # the ExtractAndStitchOpticalImageJ alignment,
    # we can evaluate more accurate values from the stitching positions.
    #
    # ... but this means our LowerLeft_XCoordinates, instead of
    # just being a function of position, will be a function of
    # position and load (alignment slightly different at each load)
    
    
    numposns = YMotionPosns.shape[0]

    YRange=(crackminy - 2*(2*dic_radius_int + dic_span_int*dic_scalefactor_int)*dy,
            crackmaxy + 2*(2*dic_radius_int + dic_span_int*dic_scalefactor_int)*dy)
    
    shift_middleimg_lowerleft_corner_x=[]
    shift_middleimg_lowerleft_corner_y=[]

    registrationshift_x = np.zeros((numposns,ActualStressPosns.shape[1]),dtype='d')
    registrationshift_y = np.zeros((numposns,ActualStressPosns.shape[1]),dtype='d')
    # Iterate throuh the various stress positions
    for stresscnt in range(ActualStressPosns.shape[1]):

        # Load in the registration shift for this stress level.
        shift_middleimg_lowerleft_corner_x.append(numericunitsv.fromxml(_xmldoc, _xmldoc.xpathsinglecontext(_element, 'dc:shift_middleimg_lowerleft_corner_x[@stepnum_stress="%d"]' % (stresscnt))).value('m'))
        shift_middleimg_lowerleft_corner_y.append(numericunitsv.fromxml(_xmldoc, _xmldoc.xpathsinglecontext(_element, 'dc:shift_middleimg_lowerleft_corner_y[@stepnum_stress="%d"]' % (stresscnt))).value('m'))

        for posncnt in range(numposns):
            registrationshift_x[posncnt,stresscnt] = numericunitsv.fromxml(_xmldoc, _xmldoc.xpathsinglecontext(_element, 'dc:registrationshift_x[@stepnum_posn="%d" and @stepnum_stress="%d"]' % (posncnt,stresscnt))).value('pixels')
            registrationshift_y[posncnt,stresscnt] = numericunitsv.fromxml(_xmldoc, _xmldoc.xpathsinglecontext(_element, 'dc:registrationshift_y[@stepnum_posn="%d" and @stepnum_stress="%d"]' % (posncnt,stresscnt))).value('pixels')
            pass

        pass

    # build LowerLeft_XCoordinates from alignment data
    # indexed by (posncnt,stresscnt)... positions relative to  image @ middle position, first stress. 
    LowerLeft_XCoordinates = (registrationshift_x-registrationshift_x[numposns//2,:])*dx  # subtract out shift from middle image at this stress level
    
    relshift_middleimg_lowerleft_corner_x = np.array(shift_middleimg_lowerleft_corner_x,dtype='d')-shift_middleimg_lowerleft_corner_x[0]
    relshift_middleimg_lowerleft_corner_y = np.array(shift_middleimg_lowerleft_corner_y,dtype='d')-shift_middleimg_lowerleft_corner_y[0]

    
    # Now (shift_middleimg_lowerleft_corner_x[0],shift_middleimg_lowerleft_corner_y[0]) are the baseline shift (first load level)
    # and relshift_middleimg_lowerleft_corner_x and relshift_middleimg_lowerleft_corner_y are relative to them.
    
    # Add the shifts into LowerLeft_XCoordinates and ybase
    LowerLeft_XCoordinates += shift_middleimg_lowerleft_corner_x[0]
    ybase += shift_middleimg_lowerleft_corner_y[0]  

    dgs_outfilehref = hrefv(posixpath.splitext(dc_scan_outdgd_href.get_bare_quoted_filename())[0]+"_dic.dgs",contexthref=dc_scan_outdgd_href.leafless())

    extra_wfmdict={}

    

    processpool = multiprocessing.Pool(multiprocessing.cpu_count()/2+1)


    
    
    (outwfmdict,
     outmetadata,
     u_disps,
     v_disps,
     ROI_out_arrays,
     Xposvecs,
     Xinivec,
     CrackCenterX,
     dic_dx,dic_dy)=execute_dic_loaded_data(Images,dx,dy,ybase,ActualStressPosns,LowerLeft_XCoordinates,
                                            dgs_outfilehref.getpath(),dic_scalefactor_int,dic_radius_int,TipCoords1,TipCoords2,YRange,
                                            extra_wfmdict=extra_wfmdict,
                                            relshift_middleimg_lowerleft_corner_x=relshift_middleimg_lowerleft_corner_x,
                                            relshift_middleimg_lowerleft_corner_y=relshift_middleimg_lowerleft_corner_y,
                                            n_threads=n_threads,processpool=processpool,debug=debug_bool)
    
    
    return {
        "dc:scan_outdic": dgs_outfilehref,
        }
