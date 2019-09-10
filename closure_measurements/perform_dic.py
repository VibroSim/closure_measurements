import sys
import collections
import numpy as np


import scipy
import scipy.optimize
import multiprocessing

import dataguzzler as dg
import dg_file as dgf
import dg_metadata as dgm
import dg_dgdread

#import pyximport
#pyximport.install()

from . import correlate
from . import initial_fit


#dgdfilename = sys.argv[1]

#TipCoords1=(.313e-3,3.49e-3) # should have smaller value of x
#TipCoords2=(.335e-3,7.20e-3) # Should have larger value of x
#XRange=(.15e-3,.8e-3)


# X Positions used in DIC (as seen in closure_measurement_coords script)
# come from:
#    * Ignoring (setting to 0) IniValX and IniValY from .dgd file (use_x0,use_y0)
#    * Treating lower-left corner of middle image as the origin.
#    * X increases to the right, Y increases up
#    * Motion stage "Y" position relative to the first such position
#      corresponds to decreasing X.
#    * pixels correspond to increasing x by dx.

def load_dgd(dgdfilename):
    (metadatadict,wfmmetadatadict,expandedwfmdict)=dg_dgdread.dg_dgdread(dgdfilename)


    #ActualYPos=metadatadict["ActualYPos"]
    #YPos=metadatadict["YPos"]
    #Stress=metadatadict["Stress"]
    #RasterStress=metadatadict["RasterStress"]
    #ActualStress=metadatadict["ActualStress"]
    #RasterY=metadatadict["RasterY"]
    
    # Fudge metadata for reshapewfms() because
    # acquisition doesn't save it quite the right way
    metadatadict["StressPos"]=metadatadict["Stress"]
    metadatadict["ActualStressPos"]=metadatadict["ActualStress"]
    
    #Images=expandedwfmdict["GEV"].data
    
    (SortedAxisNames,
     Posns,
     ActualPosns,
     PosnUnits,
     ReshapedWfmMetadata,
     ReshapedWfmData)=dg_dgdread.reshapewfms(metadatadict,wfmmetadatadict,expandedwfmdict)

    dx=ReshapedWfmMetadata["GEV"]["Step1"][0,0] # Not sure if this should be Step1 or Step2...
    dy=ReshapedWfmMetadata["GEV"]["Step2"][0,0] # Not sure if this should be Step1 or Step2...
    assert((ReshapedWfmMetadata["GEV"]["Step1"]==dx).all())
    assert((ReshapedWfmMetadata["GEV"]["Step2"]==dy).all())


    
    y0=ReshapedWfmMetadata["GEV"]["IniVal2"][0,0]
    x0=ReshapedWfmMetadata["GEV"]["IniVal1"][0,0]

    assert(PosnUnits["Y"]=="mm")  # Y axis of motion stages moves images in -X direction
    
    Images=ReshapedWfmData["GEV"].data
    # Images is (in one case) 1200x1920x19x2
    # axes are X,Y,motion stage Y,...
    # pl.imshow(Images[:,:,2,0].T,vmin=0.3,vmax=0.5)
    
    nx=Images.shape[0]
    ny=Images.shape[1]
    nimages=Images.shape[2]
    nloads=Images.shape[3]
    use_y0=0.0

    ybase=use_y0+np.arange(ny,dtype='d')*dy

    YMotionPosns=Posns["Y"]
    StressPosns=Posns["Stress"]
    ActualStressPosns = ActualPosns["Stress"]

    # Check for reversed ActualStress 
    # First row of parameter (variable a) is StressPosns
    # Second row of parameter (variable b) is ActualStressPosns
    correlation=np.corrcoef(np.array((StressPosns,np.mean(ActualStressPosns,axis=0)),dtype='d'))[0,1]
    if correlation < 0.0:
        sys.stderr.write("Found reversed ActualStress... Correcting!\n")
        ActualStressPosns=-ActualStressPosns
        pass


    # We actually ignore the IniVal1/2 from the .dgd file
    use_x0=0
    use_y0=0 

    ## *** NOTE: Y axis as defined by motion stages and Y axis from images
    ## are flipped in the recorded data. So here we flip the Y axis from the motion stages
    XPosn_relmiddle = -(YMotionPosns-YMotionPosns[nimages//2])*1e-3 # 1e-3 converts motion stage mm into meters
    LowerLeft_XCoordinates = use_x0 + XPosn_relmiddle
    LowerLeft_YCoordinates = use_y0 + YPosn_relmiddle
    
    return (Images,x0,y0,dx,dy,nx,ny,nimages,nloads,ybase,YMotionPosns,StressPosns,ActualStressPosns,LowerLeft_XCoordinates,LowerLeft_YCoordinates)


def dic_plot_click_handler(event):
    print("Selected point (%g,%g)" % (event.xdata/1.e3,event.ydata/1.e3))
    pass

def dic_raw_plots(dgdfilename):
    from matplotlib import pyplot as pl
    (Images,x0,y0,dx,dy,nx,ny,nimages,nloads,ybase,YMotionPosns,StressPosns,ActualStressPosns,LowerLeft_XCoordinates,LowerLeft_YCoordinates)=load_dgd(dgdfilename)
    
    maxstress_idx=np.argmax(np.abs(StressPosns))
    for XMotionidx in range(LowerLeft_XCoordinates.shape[0]):
        LowerLeft_XCoordinate = LowerLeft_XCoordinates[XMotionidx]
        
        Xposvec=LowerLeft_XCoordinate + np.arange(nx,dtype='d')*dx
        #Yposvec=YPosn - y0-np.arange(ny,dtype='d')*dy
        use_y0=ybase[0]
        extent=np.array((Xposvec[0]-dx/2.0,Xposvec[-1]+dx/2.0,use_y0-dy/2.0,use_y0+ny*dy-dy/2.0,))*1e3
        fig=pl.figure()
        pl.imshow(Images[:,:,XMotionidx,maxstress_idx].T,origin='lower',extent=extent)
        pl.xlabel('X position')
        pl.ylabel('Y position')
        pl.title("XMotionidx=%d; LowerLeft_XCoordinate=%f mm" % (XMotionidx,LowerLeft_XCoordinate*1e3))

        fig.canvas.mpl_connect('button_press_event',dic_plot_click_handler)
        pass
    pass


def execute_one_dic(params):
    (idx2,input1,input2,ROI_buf,dic_scalefactor,dic_radius,n_threads,debug)=params
    
    (v_array,u_array,ROI_out_array) = correlate.correlate(input1,input2,ROI_buf,dic_scalefactor,dic_radius,n_threads=n_threads,debug=debug)

    return (idx2,v_array,u_array,ROI_out_array)


def execute_dic_loaded_data(Images,dx,dy,ybase,ActualStressPosns,LowerLeft_XCoordinates,LowerLeft_YCoordinates
                            dgs_outfilename,dic_scalefactor,dic_radius,TipCoords1,TipCoords2,YRange,extra_wfmdict={},relshift_middleimg_lowerleft_corner_x=None,relshift_middleimg_lowerleft_corner_y=None,motioncontroller_tiptolerance=0.0,n_threads=multiprocessing.cpu_count(),processpool=None,debug=True):
    """ Perform DIC on data already loaded into memory """
    

    #dgs_outfilename=os.path.splitext(dgdfilename)[0]+"_dic.dgs"

    nx=Images.shape[0]
    ny=Images.shape[1]
    nimages=Images.shape[2]
    nloads=Images.shape[3]

    print("Perform_dic: Got %d %dx%d images at %d loads" % (Images.shape[2],nx,ny,nloads))
    


    CrackCenterX=(TipCoords1[0]+TipCoords2[0])/2.0

    

    #sys.modules["__main__"].__dict__.update(globals())
    #sys.modules["__main__"].__dict__.update(locals())
    #raise ValueError("Break")

    # XRange selects images  where the right hand edge of each image
    # must be to the right of the left tip, and the left hand edge of
    # each image must be to the left of the right tip
    if len(LowerLeft_XCoordinates.shape) > 1:
        XRange=(np.mean(LowerLeft_XCoordinates,axis=1)+nx*dx+motioncontroller_tiptolerance > TipCoords1[0]) & (np.mean(LowerLeft_XCoordinates,axis=1)-motioncontroller_tiptolerance < TipCoords2[0])
        pass
    else:
        XRange=(LowerLeft_XCoordinates+nx*dx+motioncontroller_tiptolerance > TipCoords1[0]) & (LowerLeft_XCoordinates-motioncontroller_tiptolerance < TipCoords2[0])
        pass
    XRangeSize=np.count_nonzero(XRange)

    dic_ny = ny//dic_scalefactor
    dic_nx = nx//dic_scalefactor

    dic_dx = dx*dic_scalefactor
    dic_dy = dy*dic_scalefactor
    
    u_disps=np.ones((dic_nx,dic_ny,nloads,nloads,XRangeSize),dtype='f',order='F') # matrix of Python Objects to store u displacements
    u_disps[...]=np.nan
    v_disps=np.ones((dic_nx,dic_ny,nloads,nloads,XRangeSize),dtype='f',order='F') # matrix of Python Objects to store v displacements
    v_disps[...]=np.nan
    ROI_out_arrays=np.ones((dic_nx,dic_ny,nloads,nloads,XRangeSize),dtype='f',order='F')
    ROI_out_arrays[...]=np.nan

    load1=np.zeros((nloads,nloads,XRangeSize),dtype='f',order='F')
    load1[...]=np.nan
    load2=np.zeros((nloads,nloads,XRangeSize),dtype='f',order='F')
    load2[...]=np.nan
    
    XRange_idxs = np.where(XRange)[0]
    Xposvecs=np.ones((dic_nx,XRangeSize,nloads),dtype='f',order='F')
    Xposvecs[...]=np.nan

    Xinivec = np.ones((XRangeSize,nloads),dtype='f')
    Xinivec[...]=np.nan


    if relshift_middleimg_lowerleft_corner_x is not None:
        relxmtx_ref=np.zeros((nloads,nloads),dtype='f',order='F')
        relxmtx_diff=np.zeros((nloads,nloads),dtype='f',order='F')
        pass
    
    if relshift_middleimg_lowerleft_corner_y is not None:
        relymtx_ref=np.zeros((nloads,nloads),dtype='f',order='F')
        relymtx_diff=np.zeros((nloads,nloads),dtype='f',order='F')
        pass

    #import pdb
    #pdb.set_trace()
    
    for XCnt in range(XRange_idxs.shape[0]):
        #if YCnt != 1:
        #    continue
        Xidx = XRange_idxs[XCnt]
        for idx1 in range(nloads):

            if len(LowerLeft_XCoordinates.shape)==1:  # just single dimension -- no load dependence
                LowerLeft_XCoordinate = LowerLeft_XCoordinates[Xidx]
                pass
            else:
                LowerLeft_XCoordinate = LowerLeft_XCoordinates[Xidx,idx1]
                pass

            if len(LowerLeft_YCoordinates.shape)==1:  # just single dimension -- no load dependence
                LowerLeft_YCoordinate = LowerLeft_YCoordinates[Xidx]
                pass
            else:
                LowerLeft_YCoordinate = LowerLeft_YCoordinates[Xidx,idx1]
                pass

            Xinivec[XCnt,idx1]=LowerLeft_XCoordinate
            Xposvec=Xinivec[XCnt,idx1] + np.arange(nx//dic_scalefactor,dtype='d')*dx*dic_scalefactor
            Xposvecs[:,XCnt,idx1]=Xposvec

            #if idx1 != 0:
            #    continue

            # build up correlate_params so we can use map() (potentially with multiprocessing)
            print("X=%d/%d; idx1=%d/%d" % (XCnt,XRange_idxs.shape[0],idx1,nloads))

            correlate_params=[]
            
            for idx2 in range(idx1+1,nloads):
                #if idx2 != nloads-1:
                #    continue
                load1[idx1,idx2,XCnt]=ActualStressPosns[Xidx,idx1]
                load2[idx1,idx2,XCnt]=ActualStressPosns[Xidx,idx2]
                load1[idx2,idx1,XCnt]=ActualStressPosns[Xidx,idx2]
                load2[idx2,idx1,XCnt]=ActualStressPosns[Xidx,idx1]

                if relshift_middleimg_lowerleft_corner_x is not None:
                    relxmtx_ref[idx1,idx2]=relshift_middleimg_lowerleft_corner_x[idx1]
                    relxmtx_diff[idx1,idx2]=relshift_middleimg_lowerleft_corner_x[idx2]
                    #relxmtx_ref[idx2,idx1]=relshift_middleimg_lowerleft_corner_x[idx1]
                    pass
                
                if relshift_middleimg_lowerleft_corner_y is not None:
                    relymtx_ref[idx1,idx2]=relshift_middleimg_lowerleft_corner_y[idx1]
                    relymtx_diff[idx1,idx2]=relshift_middleimg_lowerleft_corner_y[idx2]
                    #relymtx_ref[idx2,idx1]=relshift_middleimg_lowerleft_corner_y[idx1]
                    pass


                ROI=np.zeros((nx,ny),dtype=np.uint8,order='F')
                # Should use better process to locate crack and identify ROI
                #ROI[450:840,:]=1

                # Y coordinates are generally relative to the reference stress level picked 
                # when identifying the crack tips
    
                ROI_yminidx=np.where(ybase + (shift_middleimg_lowerleft_corner_y[idx1]+shift_middleimg_lowerleft_corner_y[idx2])/2.0 > YRange[0])[0][0]
                
                ROI_ymaxidx=np.where(ybase + (shift_middleimg_lowerleft_corner_y[idx1]+shift_middleimg_lowerleft_corner_y[idx2])/2.0 < YRange[1])[0][-1]
                
                ROI[:,ROI_yminidx:ROI_ymaxidx]=1
                

                # DIC represents idx2 state minus idx1 state
                #
                # Relxmtx_ref values are positions of idx1 state
                # Relxmtx_diff values are positions of idx2 state

                # Expect DIC displacements to approximately match relxmtx_diff - relxmtx_ref
                
                input1=np.asfortranarray(Images[:,:,Xidx,idx1].T.astype(np.float64))
                input2=np.asfortranarray(Images[:,:,Xidx,idx2].T.astype(np.float64))
                ROI_buf=np.asfortranarray(ROI[:,:].T.astype(np.uint8))

                correlate_params.append((idx2,input1,input2,ROI_buf,dic_scalefactor,dic_radius,n_threads,debug))
                pass

            if processpool is None:
                correlate_results=map(execute_one_dic,correlate_params)
                pass
            else:
                correlate_results=processpool.map(execute_one_dic,correlate_params)
                pass
            for (idx2,v_array,u_array,ROI_out_array) in correlate_results:
                
                u_disps[:,:,idx1,idx2,XCnt]=u_array.T*dx
                v_disps[:,:,idx1,idx2,XCnt]=v_array.T*dy
                
                #u_disps[:,:,idx2,idx1,XCnt]=-u_array.T*dx
                #v_disps[:,:,idx2,idx1,XCnt]=-v_array.T*dy

                ROI_out_arrays[:,:,idx1,idx2,XCnt] = ROI_out_array.T
                ROI_out_arrays[:,:,idx2,idx1,XCnt] = ROI_out_array.T
                pass
            
            pass
        #break
        pass
    
    # Write output to a .dgs file


    outwfmdict=collections.OrderedDict()

    for XCnt in range(XRange_idxs.shape[0]):
        
        outwfmdict["u_disps%.3d" % (XCnt)]=dg.wfminfo()
        outwfmdict["u_disps%.3d" % (XCnt)].Name="u_disps%.3d" % (XCnt)
        outwfmdict["u_disps%.3d" % (XCnt)].data=u_disps[...,XCnt]
        outwfmdict["u_disps%.3d" % (XCnt)].dimlen=np.array(u_disps.shape[:-1])
        outwfmdict["u_disps%.3d" % (XCnt)].ndim=4
        dgm.AddMetaDatumWI(outwfmdict["u_disps%.3d" % (XCnt)],dgm.CreateMetaDatumStr("Coord1","X Position"))
        dgm.AddMetaDatumWI(outwfmdict["u_disps%.3d" % (XCnt)],dgm.CreateMetaDatumStr("Units1","meters"))
        dgm.AddMetaDatumWI(outwfmdict["u_disps%.3d" % (XCnt)],dgm.CreateMetaDatumDbl("IniVal1",np.mean(Xinivec[XCnt,:])))
        dgm.AddMetaDatumWI(outwfmdict["u_disps%.3d" % (XCnt)],dgm.CreateMetaDatumDbl("Step1",dx*dic_scalefactor))
        dgm.AddMetaDatumWI(outwfmdict["u_disps%.3d" % (XCnt)],dgm.CreateMetaDatumStr("Coord2","Y Position"))
        dgm.AddMetaDatumWI(outwfmdict["u_disps%.3d" % (XCnt)],dgm.CreateMetaDatumStr("Units2","meters"))
        #dgm.AddMetaDatumWI(outwfmdict["u_disps%.3d" % (XCnt)],dgm.CreateMetaDatumDbl("IniVal2",y0))
        dgm.AddMetaDatumWI(outwfmdict["u_disps%.3d" % (XCnt)],dgm.CreateMetaDatumDbl("IniVal2",ybase[0]))
        dgm.AddMetaDatumWI(outwfmdict["u_disps%.3d" % (XCnt)],dgm.CreateMetaDatumDbl("Step2",dy*dic_scalefactor))
        dgm.AddMetaDatumWI(outwfmdict["u_disps%.3d" % (XCnt)],dgm.CreateMetaDatumStr("Coord3","Stress Level DIC input 1"))
        dgm.AddMetaDatumWI(outwfmdict["u_disps%.3d" % (XCnt)],dgm.CreateMetaDatumStr("Coord4","Stress Level DIC input 2"))
        dgm.AddMetaDatumWI(outwfmdict["u_disps%.3d" % (XCnt)],dgm.CreateMetaDatumStr("Coord5","Image index"))
        dgm.AddMetaDatumWI(outwfmdict["u_disps%.3d" % (XCnt)],dgm.CreateMetaDatumStr("Units5","Unitless"))
        dgm.AddMetaDatumWI(outwfmdict["u_disps%.3d" % (XCnt)],dgm.CreateMetaDatumStr("AmplUnits","meters"))
        dgm.AddMetaDatumWI(outwfmdict["u_disps%.3d" % (XCnt)],dgm.CreateMetaDatumStr("X Displacement","meters"))
        
        
        outwfmdict["v_disps%.3d" % (XCnt)]=dg.wfminfo()
        outwfmdict["v_disps%.3d" % (XCnt)].Name="v_disps%.3d" % (XCnt)
        outwfmdict["v_disps%.3d" % (XCnt)].data=v_disps[...,XCnt]
        outwfmdict["v_disps%.3d" % (XCnt)].dimlen=np.array(v_disps.shape[:-1])
        outwfmdict["v_disps%.3d" % (XCnt)].ndim=4
        dgm.AddMetaDatumWI(outwfmdict["v_disps%.3d" % (XCnt)],dgm.CreateMetaDatumStr("Coord1","X Position"))
        dgm.AddMetaDatumWI(outwfmdict["v_disps%.3d" % (XCnt)],dgm.CreateMetaDatumStr("Units1","meters"))
        dgm.AddMetaDatumWI(outwfmdict["v_disps%.3d" % (XCnt)],dgm.CreateMetaDatumDbl("IniVal1",np.mean(Xinivec[XCnt,:])))
        dgm.AddMetaDatumWI(outwfmdict["v_disps%.3d" % (XCnt)],dgm.CreateMetaDatumDbl("Step1",dx*dic_scalefactor))
        dgm.AddMetaDatumWI(outwfmdict["v_disps%.3d" % (XCnt)],dgm.CreateMetaDatumStr("Coord2","Y Position"))
        dgm.AddMetaDatumWI(outwfmdict["v_disps%.3d" % (XCnt)],dgm.CreateMetaDatumStr("Units2","meters"))
        #dgm.AddMetaDatumWI(outwfmdict["v_disps%.3d" % (XCnt)],dgm.CreateMetaDatumDbl("IniVal2",y0))
        dgm.AddMetaDatumWI(outwfmdict["v_disps%.3d" % (XCnt)],dgm.CreateMetaDatumDbl("IniVal2",ybase[0]))
        dgm.AddMetaDatumWI(outwfmdict["v_disps%.3d" % (XCnt)],dgm.CreateMetaDatumDbl("Step2",dy*dic_scalefactor))
        dgm.AddMetaDatumWI(outwfmdict["v_disps%.3d" % (XCnt)],dgm.CreateMetaDatumStr("Coord3","Stress Level DIC input 1"))
        dgm.AddMetaDatumWI(outwfmdict["v_disps%.3d" % (XCnt)],dgm.CreateMetaDatumStr("Coord4","Stress Level DIC input 2"))
        dgm.AddMetaDatumWI(outwfmdict["v_disps%.3d" % (XCnt)],dgm.CreateMetaDatumStr("Coord5","Image index"))
        dgm.AddMetaDatumWI(outwfmdict["v_disps%.3d" % (XCnt)],dgm.CreateMetaDatumStr("Units5","Unitless"))
        dgm.AddMetaDatumWI(outwfmdict["v_disps%.3d" % (XCnt)],dgm.CreateMetaDatumStr("AmplUnits","meters"))
        dgm.AddMetaDatumWI(outwfmdict["v_disps%.3d" % (XCnt)],dgm.CreateMetaDatumStr("Y Displacement","meters"))
        

        outwfmdict["ROI_out%.3d" % (XCnt)]=dg.wfminfo()
        outwfmdict["ROI_out%.3d" % (XCnt)].Name="ROI_out%.3d" % (XCnt)
        outwfmdict["ROI_out%.3d" % (XCnt)].data=ROI_out_arrays[...,XCnt]
        outwfmdict["ROI_out%.3d" % (XCnt)].dimlen=np.array(ROI_out_arrays.shape[:-1])
        outwfmdict["ROI_out%.3d" % (XCnt)].ndim=4
        dgm.AddMetaDatumWI(outwfmdict["ROI_out%.3d" % (XCnt)],dgm.CreateMetaDatumStr("Coord1","X Position"))
        dgm.AddMetaDatumWI(outwfmdict["ROI_out%.3d" % (XCnt)],dgm.CreateMetaDatumStr("Units1","meters"))
        dgm.AddMetaDatumWI(outwfmdict["ROI_out%.3d" % (XCnt)],dgm.CreateMetaDatumDbl("IniVal1",np.mean(Xinivec[XCnt,:])))
        dgm.AddMetaDatumWI(outwfmdict["ROI_out%.3d" % (XCnt)],dgm.CreateMetaDatumDbl("Step1",dx*dic_scalefactor))
        dgm.AddMetaDatumWI(outwfmdict["ROI_out%.3d" % (XCnt)],dgm.CreateMetaDatumStr("Coord2","Y Position"))
        dgm.AddMetaDatumWI(outwfmdict["ROI_out%.3d" % (XCnt)],dgm.CreateMetaDatumStr("Units2","meters"))
        #dgm.AddMetaDatumWI(outwfmdict["ROI_out%.3d" % (XCnt)],dgm.CreateMetaDatumDbl("IniVal2",y0))
        dgm.AddMetaDatumWI(outwfmdict["ROI_out%.3d" % (XCnt)],dgm.CreateMetaDatumDbl("IniVal2",ybase[0]))
        dgm.AddMetaDatumWI(outwfmdict["ROI_out%.3d" % (XCnt)],dgm.CreateMetaDatumDbl("Step2",dy*dic_scalefactor))
        dgm.AddMetaDatumWI(outwfmdict["ROI_out%.3d" % (XCnt)],dgm.CreateMetaDatumStr("Coord3","Stress Level DIC input 1"))
        dgm.AddMetaDatumWI(outwfmdict["ROI_out%.3d" % (XCnt)],dgm.CreateMetaDatumStr("Coord4","Stress Level DIC input 2"))
        dgm.AddMetaDatumWI(outwfmdict["ROI_out%.3d" % (XCnt)],dgm.CreateMetaDatumStr("Coord5","Image index"))
        dgm.AddMetaDatumWI(outwfmdict["ROI_out%.3d" % (XCnt)],dgm.CreateMetaDatumStr("Units5","Unitless"))

        pass
    
    outwfmdict["Xposvecs"]=dg.wfminfo()
    outwfmdict["Xposvecs"].Name="Xposvecs"
    outwfmdict["Xposvecs"].data=Xposvecs
    outwfmdict["Xposvecs"].dimlen=np.array(Xposvecs.shape)
    outwfmdict["Xposvecs"].ndim=3
    dgm.AddMetaDatumWI(outwfmdict["Xposvecs"],dgm.CreateMetaDatumStr("Coord1","X Position"))
    dgm.AddMetaDatumWI(outwfmdict["Xposvecs"],dgm.CreateMetaDatumStr("Units1","meters"))
    #dgm.AddMetaDatumWI(outwfmdict["Xposvecs"],dgm.CreateMetaDatumDbl("IniVal1",x0))
    dgm.AddMetaDatumWI(outwfmdict["Xposvecs"],dgm.CreateMetaDatumDbl("IniVal1",0.0))
    dgm.AddMetaDatumWI(outwfmdict["Xposvecs"],dgm.CreateMetaDatumDbl("Step1",dx*dic_scalefactor))
    dgm.AddMetaDatumWI(outwfmdict["Xposvecs"],dgm.CreateMetaDatumStr("Coord2","Image index"))
    dgm.AddMetaDatumWI(outwfmdict["Xposvecs"],dgm.CreateMetaDatumStr("Units2","Unitless"))
    dgm.AddMetaDatumWI(outwfmdict["Xposvecs"],dgm.CreateMetaDatumStr("Coord3","Load index"))
    dgm.AddMetaDatumWI(outwfmdict["Xposvecs"],dgm.CreateMetaDatumStr("Units3","Unitless"))
    dgm.AddMetaDatumWI(outwfmdict["Xposvecs"],dgm.CreateMetaDatumStr("AmplCoord","X shiftedscaled"))
    dgm.AddMetaDatumWI(outwfmdict["Xposvecs"],dgm.CreateMetaDatumStr("AmplUnits","meters"))


    outwfmdict["Xinivec"]=dg.wfminfo()
    outwfmdict["Xinivec"].Name="Xinivec"
    outwfmdict["Xinivec"].data=Xinivec
    outwfmdict["Xinivec"].dimlen=np.array(Xinivec.shape)
    outwfmdict["Xinivec"].ndim=2
    dgm.AddMetaDatumWI(outwfmdict["Xinivec"],dgm.CreateMetaDatumStr("Coord1","Image index"))
    dgm.AddMetaDatumWI(outwfmdict["Xinivec"],dgm.CreateMetaDatumStr("Units1","Unitless"))
    dgm.AddMetaDatumWI(outwfmdict["Xinivec"],dgm.CreateMetaDatumStr("Coord2","Load index"))
    dgm.AddMetaDatumWI(outwfmdict["Xinivec"],dgm.CreateMetaDatumStr("Units2","Unitless"))
    dgm.AddMetaDatumWI(outwfmdict["Xinivec"],dgm.CreateMetaDatumStr("AmplCoord","Y initial"))
    dgm.AddMetaDatumWI(outwfmdict["Xinivec"],dgm.CreateMetaDatumStr("AmplUnits","meters"))


    outwfmdict["load1"]=dg.wfminfo()
    outwfmdict["load1"].Name="load1"
    outwfmdict["load1"].data=load1
    outwfmdict["load1"].dimlen=np.array(load1.shape)
    outwfmdict["load1"].ndim=3
    dgm.AddMetaDatumWI(outwfmdict["load1"],dgm.CreateMetaDatumStr("Coord1","First load index"))
    dgm.AddMetaDatumWI(outwfmdict["load1"],dgm.CreateMetaDatumStr("Units1","Unitless"))
    dgm.AddMetaDatumWI(outwfmdict["load1"],dgm.CreateMetaDatumStr("Coord2","Second load index"))
    dgm.AddMetaDatumWI(outwfmdict["load1"],dgm.CreateMetaDatumStr("Units2","Unitless"))
    dgm.AddMetaDatumWI(outwfmdict["load1"],dgm.CreateMetaDatumStr("Coord3","Image index"))
    dgm.AddMetaDatumWI(outwfmdict["load1"],dgm.CreateMetaDatumStr("Units3","Unitless"))
    dgm.AddMetaDatumWI(outwfmdict["load1"],dgm.CreateMetaDatumStr("AmplCoord","Stress"))
    dgm.AddMetaDatumWI(outwfmdict["load1"],dgm.CreateMetaDatumStr("AmplUnits","Pascals"))

    outwfmdict["load2"]=dg.wfminfo()
    outwfmdict["load2"].Name="load2"
    outwfmdict["load2"].data=load2
    outwfmdict["load2"].dimlen=np.array(load2.shape)
    outwfmdict["load2"].ndim=3
    dgm.AddMetaDatumWI(outwfmdict["load2"],dgm.CreateMetaDatumStr("Coord1","First load index"))
    dgm.AddMetaDatumWI(outwfmdict["load2"],dgm.CreateMetaDatumStr("Units1","Unitless"))
    dgm.AddMetaDatumWI(outwfmdict["load2"],dgm.CreateMetaDatumStr("Coord2","Second load index"))
    dgm.AddMetaDatumWI(outwfmdict["load2"],dgm.CreateMetaDatumStr("Units2","Unitless"))
    dgm.AddMetaDatumWI(outwfmdict["load2"],dgm.CreateMetaDatumStr("Coord3","Image index"))
    dgm.AddMetaDatumWI(outwfmdict["load2"],dgm.CreateMetaDatumStr("Units3","Unitless"))
    dgm.AddMetaDatumWI(outwfmdict["load2"],dgm.CreateMetaDatumStr("AmplCoord","Stress"))
    dgm.AddMetaDatumWI(outwfmdict["load2"],dgm.CreateMetaDatumStr("AmplUnits","Pascals"))

    
    outmetadata={}
    dgm.AddMetaDatumL(outmetadata,dgm.CreateMetaDatumDbl("CrackCenterX",CrackCenterX))
    dgm.AddMetaDatumL(outmetadata,dgm.CreateMetaDatumDbl("TipCoords1X",TipCoords1[0]))
    dgm.AddMetaDatumL(outmetadata,dgm.CreateMetaDatumDbl("TipCoords1Y",TipCoords1[1]))
    dgm.AddMetaDatumL(outmetadata,dgm.CreateMetaDatumDbl("TipCoords2X",TipCoords2[0]))
    dgm.AddMetaDatumL(outmetadata,dgm.CreateMetaDatumDbl("TipCoords2Y",TipCoords2[1]))

    dgm.AddMetaDatumL(outmetadata,dgm.CreateMetaDatumInt("NumImages",XRange_idxs.shape[0]))

    dgm.AddMetaDatumL(outmetadata,dgm.CreateMetaDatumInt("ROI_yminidx",ROI_yminidx))
    dgm.AddMetaDatumL(outmetadata,dgm.CreateMetaDatumInt("ROI_ymaxidx",ROI_ymaxidx))

    dgm.AddMetaDatumL(outmetadata,dgm.CreateMetaDatumInt("ROI_dic_yminidx",ROI_yminidx//dic_scalefactor))
    dgm.AddMetaDatumL(outmetadata,dgm.CreateMetaDatumInt("ROI_dic_ymaxidx",ROI_ymaxidx//dic_scalefactor))


    if relshift_middleimg_lowerleft_corner_x is not None:
       relxwfm_ref = dg.wfminfo()
       relxwfm_ref.Name = "relshift_middleimg_lowerleft_corner_x_ref"
       relxwfm_ref.data = relxmtx_ref
       relxwfm_ref.dimlen = np.array(relxmtx_ref.shape)
       relxwfm_ref.ndim=2
       dgm.AddMetaDatumWI(relxwfm_ref,dgm.CreateMetaDatumStr("Coord1","First load index"))
       dgm.AddMetaDatumWI(relxwfm_ref,dgm.CreateMetaDatumStr("Units1","Unitless"))
       dgm.AddMetaDatumWI(relxwfm_ref,dgm.CreateMetaDatumStr("Coord2","Second load index"))
       dgm.AddMetaDatumWI(relxwfm_ref,dgm.CreateMetaDatumStr("Units2","Unitless"))
       dgm.AddMetaDatumWI(relxwfm_ref,dgm.CreateMetaDatumStr("AmplUnits","meters"))
       dgm.AddMetaDatumWI(relxwfm_ref,dgm.CreateMetaDatumStr("AmplCoord","X shift"))
       outwfmdict["relshift_middleimg_lowerleft_corner_x_ref"]=relxwfm_ref

       relxwfm_diff = dg.wfminfo()
       relxwfm_diff.Name = "relshift_middleimg_lowerleft_corner_x_diff"
       relxwfm_diff.data = relxmtx_diff
       relxwfm_diff.dimlen = np.array(relxmtx_diff.shape)
       relxwfm_diff.ndim=2
       dgm.AddMetaDatumWI(relxwfm_diff,dgm.CreateMetaDatumStr("Coord1","First load index"))
       dgm.AddMetaDatumWI(relxwfm_diff,dgm.CreateMetaDatumStr("Units1","Unitless"))
       dgm.AddMetaDatumWI(relxwfm_diff,dgm.CreateMetaDatumStr("Coord2","Second load index"))
       dgm.AddMetaDatumWI(relxwfm_diff,dgm.CreateMetaDatumStr("Units2","Unitless"))
       dgm.AddMetaDatumWI(relxwfm_diff,dgm.CreateMetaDatumStr("AmplUnits","meters"))
       dgm.AddMetaDatumWI(relxwfm_diff,dgm.CreateMetaDatumStr("AmplCoord","X shift"))
       outwfmdict["relshift_middleimg_lowerleft_corner_x_diff"]=relxwfm_diff

       pass

    if relshift_middleimg_lowerleft_corner_y is not None:
        relywfm_ref = dg.wfminfo()
        relywfm_ref.Name = "relshift_middleimg_lowerleft_corner_y_ref"
        relywfm_ref.data = relymtx_ref
        relywfm_ref.dimlen = np.array(relymtx_ref.shape)
        relywfm_ref.ndim=2
        dgm.AddMetaDatumWI(relywfm_ref,dgm.CreateMetaDatumStr("Coord1","First load index"))
        dgm.AddMetaDatumWI(relywfm_ref,dgm.CreateMetaDatumStr("Units1","Unitless"))
        dgm.AddMetaDatumWI(relywfm_ref,dgm.CreateMetaDatumStr("Coord2","Second load index"))
        dgm.AddMetaDatumWI(relywfm_ref,dgm.CreateMetaDatumStr("Units2","Unitless"))
        dgm.AddMetaDatumWI(relywfm_ref,dgm.CreateMetaDatumStr("AmplUnits","meters"))
        dgm.AddMetaDatumWI(relywfm_ref,dgm.CreateMetaDatumStr("AmplCoord","Y shift"))
        outwfmdict["relshift_middleimg_lowerleft_corner_y_ref"]=relywfm_ref

        relywfm_diff = dg.wfminfo()
        relywfm_diff.Name = "relshift_middleimg_lowerleft_corner_y_diff"
        relywfm_diff.data = relymtx_diff
        relywfm_diff.dimlen = np.array(relymtx_diff.shape)
        relywfm_diff.ndim=2
        dgm.AddMetaDatumWI(relywfm_diff,dgm.CreateMetaDatumStr("Coord1","First load index"))
        dgm.AddMetaDatumWI(relywfm_diff,dgm.CreateMetaDatumStr("Units1","Unitless"))
        dgm.AddMetaDatumWI(relywfm_diff,dgm.CreateMetaDatumStr("Coord2","Second load index"))
        dgm.AddMetaDatumWI(relywfm_diff,dgm.CreateMetaDatumStr("Units2","Unitless"))
        dgm.AddMetaDatumWI(relywfm_diff,dgm.CreateMetaDatumStr("AmplUnits","meters"))
        dgm.AddMetaDatumWI(relywfm_diff,dgm.CreateMetaDatumStr("AmplCoord","Y shift"))
        outwfmdict["relshift_middleimg_lowerleft_corner_y_diff"]=relywfm_diff
        pass

    
    # Add extra waveforms provided by caller
    for key in extra_wfmdict:
        outwfmdict[key]=extra_wfmdict[key]
        pass
    
    
    dgfh=dgf.creat(dgs_outfilename)
    dgf.writesnapshot(dgfh,outmetadata,outwfmdict)
    dgf.close(dgfh)

    return (outwfmdict,outmetadata,u_disps,v_disps,ROI_out_arrays,Xposvecs,Xinivec,CrackCenterX,dic_dx,dic_dy)
    


    
    


def execute_dic(dgdfilename,dgs_outfilename,dic_scalefactor,dic_radius,TipCoords1,TipCoords2,YRange,extra_wfmdict={},n_threads=multiprocessing.cpu_count(),processpool=None,debug=True):
    """Perform DIC on optical microscopy .dgd file. 
     dic_scalefactor and dic_radius parameters to ncorr, given in pixels
     TipCoords1 is an (x,y) tuple indicating the coordinates of the tip with a 
     lower value of y, in meters
     TipCoords2 is an (x,y) tuple indicating the coordinates of the tip with a
     larger value of y, in meters
     YRange is a (y1,y2) tuple indicating the lower and upper bounds of the region of 
     of interest in y, in meters

"""
    
    (Images,x0,y0,dx,dy,nx,ny,nimages,nloads,ybase,YMotionPosns,StressPosns,ActualStressPosns,LowerLeft_XCoordinates,LowerLeft_YCoordinates)=load_dgd(dgdfilename)

    return execute_dic_loaded_data(Images,dx,dy,ybase,ActualStressPosns,LowerLeft_XCoordinates,LowerLeft_YCoordinates,
                                   dgs_outfilename,dic_scalefactor,dic_radius,TipCoords1,TipCoords2,YRange,extra_wfmdict=extra_wfmdict,n_threads=n_threads,processpool=processpool,debug=debug)
