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
from . import dic_ctod
from . import initial_fit


#dgdfilename = sys.argv[1]

#TipCoords1=(.313e-3,3.49e-3) # should have smaller value of x
#TipCoords2=(.335e-3,7.20e-3) # Should have larger value of x
#XRange=(.15e-3,.8e-3)

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
    assert((ReshapedWfmMetadata["GEV"]["Step1"]==dy).all())
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

    
    return (Images,x0,y0,dx,dy,nx,ny,nimages,nloads,ybase,YMotionPosns,StressPosns,ActualStressPosns)


def dic_plot_click_handler(event):
    print("Selected point (%g,%g)" % (event.xdata/1.e3,event.ydata/1.e3))
    pass

def dic_raw_plots(dgdfilename):
    from matplotlib import pyplot as pl
    (Images,x0,y0,dx,dy,nx,ny,nimages,nloads,ybase,YMotionPosns,StressPosns,ActualStressPosns)=load_dgd(dgdfilename)
    
    maxstress_idx=np.argmax(np.abs(StressPosns))
    for YMotionidx in range(YMotionPosns.shape[0]):
        YMotionPosn = (YMotionPosns[YMotionidx]-YMotionPosns[0])*1e-3 # Posns is stored in mm
        
        ## *** NOTE: Y axis as defined by motion stages and X axis from images 
        ## are flipped in the recorded data. So here we flip the Y axis from the motion stages so it aligns with the X axis from the images
        ##use_y0 = y0 - (ny//dic_scalefactor-1)*dy*dic_scalefactor
        #use_y0 = y0 - (ny-1)*dy
        #use_y0=y0
        use_x0=0.0
        use_y0=0.0
        Xposvec=-YMotionPosn + use_x0 + np.arange(nx,dtype='d')*dx
        #Yposvec=YPosn - y0-np.arange(ny,dtype='d')*dy
        extent=np.array((Xposvec[0]-dx/2.0,Xposvec[-1]+dx/2.0,use_y0-dy/2.0,use_y0+ny*dy-dy/2.0,))*1e3
        fig=pl.figure()
        pl.imshow(Images[:,:,YMotionidx,maxstress_idx].T,origin='lower',extent=extent)
        pl.xlabel('X position')
        pl.ylabel('Y position')
        pl.title("YMotionidx=%d; YMotionPosn=%f mm" % (YMotionidx,YMotionPosn*1e3))

        fig.canvas.mpl_connect('button_press_event',dic_plot_click_handler)
        pass
    pass


def execute_one_dic(params):
    (idx2,input1,input2,ROI_buf,dic_scalefactor,dic_radius,n_threads,debug)=params
    
    (v_array,u_array,ROI_out_array) = correlate.correlate(input1,input2,ROI_buf,dic_scalefactor,dic_radius,n_threads=n_threads,debug=debug)

    return (idx2,v_array,u_array,ROI_out_array)

def execute_dic(dgdfilename,dgs_outfilename,dic_scalefactor,dic_radius,TipCoords1,TipCoords2,YRange,n_threads=multiprocessing.cpu_count(),processpool=None,debug=True):
    """Perform DIC on optical microscopy .dgd file. 
     dic_scalefactor and dic_radius parameters to ncorr, given in pixels
     TipCoords1 is an (x,y) tuple indicating the coordinates of the tip with a 
     lower value of y, in meters
     TipCoords2 is an (x,y) tuple indicating the coordinates of the tip with a
     larger value of y, in meters
     YRange is a (y1,y2) tuple indicating the lower and upper bounds of the region of 
     of interest in y, in meters

"""
    

    
    (Images,x0,y0,dx,dy,nx,ny,nimages,nloads,ybase,YMotionPosns,StressPosns,ActualStressPosns)=load_dgd(dgdfilename)

    #dgs_outfilename=os.path.splitext(dgdfilename)[0]+"_dic.dgs"
    

    print("Perform_dic: Got %d %dx%d images at %d loads" % (Images.shape[2],nx,ny,nloads))
    


    CrackCenterX=(TipCoords1[0]+TipCoords2[0])/2.0

    
    ROI=np.zeros((nx,ny),dtype=np.uint8,order='F')
    # Should use better process to locate crack and identify ROI
    #ROI[450:840,:]=1
    
    ROI_yminidx=np.where(ybase > YRange[0])[0][0]
    
    ROI_ymaxidx=np.where(ybase < YRange[1])[0][-1]

    ROI[:,ROI_yminidx:ROI_ymaxidx]=1

    #sys.modules["__main__"].__dict__.update(globals())
    #sys.modules["__main__"].__dict__.update(locals())
    #raise ValueError("Break")
    XRange=(-(YMotionPosns-YMotionPosns[0])*1e-3+nx*dx > TipCoords1[0]) & (-(YMotionPosns-YMotionPosns[0])*1e-3 < TipCoords2[0])
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
    Xposvecs=np.ones((dic_nx,XRangeSize),dtype='f',order='F')
    Xposvecs[...]=np.nan

    Xinivec = np.ones(XRangeSize,dtype='f')
    Xinivec[...]=np.nan

    #import pdb
    #pdb.set_trace()
    
    for XCnt in range(XRange_idxs.shape[0]):
        #if YCnt != 1:
        #    continue
        Xidx = XRange_idxs[XCnt]
        YMotionPosn = -(YMotionPosns[Xidx]-YMotionPosns[0])*1e-3 # Posns is stored in mm
        ## *** NOTE: Y axis as defined by motion stages and Y axis from images
        ## are flipped in the recorded data. So here we flip the Y axis from the motion stages

        #use_y0 = y0 - (ny-1)*dy
        #use_y0 = y0
        use_x0 = 0.0
        Xinivec[XCnt]=use_x0+YMotionPosn
        Xposvec=Xinivec[XCnt] + np.arange(nx//dic_scalefactor,dtype='d')*dx*dic_scalefactor
        Xposvecs[:,XCnt]=Xposvec
        for idx1 in range(nloads):
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
                
                u_disps[:,:,idx2,idx1,XCnt]=-u_array.T*dx
                v_disps[:,:,idx2,idx1,XCnt]=-v_array.T*dy

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
        dgm.AddMetaDatumWI(outwfmdict["u_disps%.3d" % (XCnt)],dgm.CreateMetaDatumDbl("IniVal1",Xinivec[XCnt]))
        dgm.AddMetaDatumWI(outwfmdict["u_disps%.3d" % (XCnt)],dgm.CreateMetaDatumDbl("Step1",dx*dic_scalefactor))
        dgm.AddMetaDatumWI(outwfmdict["u_disps%.3d" % (XCnt)],dgm.CreateMetaDatumStr("Coord2","Y Position"))
        dgm.AddMetaDatumWI(outwfmdict["u_disps%.3d" % (XCnt)],dgm.CreateMetaDatumStr("Units2","meters"))
        #dgm.AddMetaDatumWI(outwfmdict["u_disps%.3d" % (XCnt)],dgm.CreateMetaDatumDbl("IniVal2",y0))
        dgm.AddMetaDatumWI(outwfmdict["u_disps%.3d" % (XCnt)],dgm.CreateMetaDatumDbl("IniVal2",0.0))
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
        dgm.AddMetaDatumWI(outwfmdict["v_disps%.3d" % (XCnt)],dgm.CreateMetaDatumDbl("IniVal1",Xinivec[XCnt]))
        dgm.AddMetaDatumWI(outwfmdict["v_disps%.3d" % (XCnt)],dgm.CreateMetaDatumDbl("Step1",dx*dic_scalefactor))
        dgm.AddMetaDatumWI(outwfmdict["v_disps%.3d" % (XCnt)],dgm.CreateMetaDatumStr("Coord2","Y Position"))
        dgm.AddMetaDatumWI(outwfmdict["v_disps%.3d" % (XCnt)],dgm.CreateMetaDatumStr("Units2","meters"))
        #dgm.AddMetaDatumWI(outwfmdict["v_disps%.3d" % (XCnt)],dgm.CreateMetaDatumDbl("IniVal2",y0))
        dgm.AddMetaDatumWI(outwfmdict["v_disps%.3d" % (XCnt)],dgm.CreateMetaDatumDbl("IniVal2",0.0))
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
        dgm.AddMetaDatumWI(outwfmdict["ROI_out%.3d" % (XCnt)],dgm.CreateMetaDatumDbl("IniVal1",Xinivec[XCnt]))
        dgm.AddMetaDatumWI(outwfmdict["ROI_out%.3d" % (XCnt)],dgm.CreateMetaDatumDbl("Step1",dx*dic_scalefactor))
        dgm.AddMetaDatumWI(outwfmdict["ROI_out%.3d" % (XCnt)],dgm.CreateMetaDatumStr("Coord2","Y Position"))
        dgm.AddMetaDatumWI(outwfmdict["ROI_out%.3d" % (XCnt)],dgm.CreateMetaDatumStr("Units2","meters"))
        #dgm.AddMetaDatumWI(outwfmdict["ROI_out%.3d" % (XCnt)],dgm.CreateMetaDatumDbl("IniVal2",y0))
        dgm.AddMetaDatumWI(outwfmdict["ROI_out%.3d" % (XCnt)],dgm.CreateMetaDatumDbl("IniVal2",0.0))
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
    outwfmdict["Xposvecs"].ndim=2
    dgm.AddMetaDatumWI(outwfmdict["Xposvecs"],dgm.CreateMetaDatumStr("Coord1","X Position"))
    dgm.AddMetaDatumWI(outwfmdict["Xposvecs"],dgm.CreateMetaDatumStr("Units1","meters"))
    #dgm.AddMetaDatumWI(outwfmdict["Xposvecs"],dgm.CreateMetaDatumDbl("IniVal1",x0))
    dgm.AddMetaDatumWI(outwfmdict["Xposvecs"],dgm.CreateMetaDatumDbl("IniVal1",0.0))
    dgm.AddMetaDatumWI(outwfmdict["Xposvecs"],dgm.CreateMetaDatumDbl("Step1",dx*dic_scalefactor))
    dgm.AddMetaDatumWI(outwfmdict["Xposvecs"],dgm.CreateMetaDatumStr("Coord2","Image index"))
    dgm.AddMetaDatumWI(outwfmdict["Xposvecs"],dgm.CreateMetaDatumStr("Units2","Unitless"))
    dgm.AddMetaDatumWI(outwfmdict["Xposvecs"],dgm.CreateMetaDatumStr("AmplCoord","X shiftedscaled"))
    dgm.AddMetaDatumWI(outwfmdict["Xposvecs"],dgm.CreateMetaDatumStr("AmplUnits","meters"))


    outwfmdict["Xinivec"]=dg.wfminfo()
    outwfmdict["Xinivec"].Name="Xinivec"
    outwfmdict["Xinivec"].data=Xinivec
    outwfmdict["Xinivec"].dimlen=np.array(Xinivec.shape)
    outwfmdict["Xinivec"].ndim=1
    dgm.AddMetaDatumWI(outwfmdict["Xinivec"],dgm.CreateMetaDatumStr("Coord1","Image index"))
    dgm.AddMetaDatumWI(outwfmdict["Xinivec"],dgm.CreateMetaDatumStr("Units1","Unitless"))
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


    dgfh=dgf.creat(dgs_outfilename)
    dgf.writesnapshot(dgfh,outmetadata,outwfmdict)
    dgf.close(dgfh)

    return (outwfmdict,outmetadata,u_disps,v_disps,ROI_out_arrays,Xposvecs,Xinivec,CrackCenterX,dic_dx,dic_dy)
    
