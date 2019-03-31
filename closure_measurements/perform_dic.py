import sys
import collections
import numpy as np

from matplotlib import pyplot as pl

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

TipCoords1=(.313e-3,3.49e-3) # should have smaller value of y
TipCoords2=(.335e-3,7.20e-3) # Should have larger value of y
XRange=(.15e-3,.8e-3)

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

    dy=ReshapedWfmMetadata["GEV"]["Step1"][0,0] # Not sure if this should be Step1 or Step2...
    dx=ReshapedWfmMetadata["GEV"]["Step2"][0,0] # Not sure if this should be Step1 or Step2...
    assert((ReshapedWfmMetadata["GEV"]["Step1"]==dy).all())
    assert((ReshapedWfmMetadata["GEV"]["Step2"]==dy).all())


    
    x0=ReshapedWfmMetadata["GEV"]["IniVal2"][0,0]
    y0=ReshapedWfmMetadata["GEV"]["IniVal1"][0,0]

    assert(PosnUnits["Y"]=="mm")
    
    Images=ReshapedWfmData["GEV"].data
    # Images is (in one case) 1200x1920x19x2
    # axes are Y,X,...
    # pl.imshow(Images[:,:,2,0].T,vmin=0.3,vmax=0.5)
    
    ny=Images.shape[0]
    nx=Images.shape[1]
    nimages=Images.shape[2]
    nloads=Images.shape[3]
    
    xbase=x0+np.arange(nx,dtype='d')*dx

    YPosns=Posns["Y"]
    StressPosns=Posns["Stress"]
    ActualStressPosns = ActualPosns["Stress"]
    
    return (Images,y0,x0,dy,dx,ny,nx,nimages,nloads,xbase,YPosns,StressPosns,ActualStressPosns)


def dic_plot_click_handler(event):
    print("Selected point (%g,%g)" % (event.xdata/1.e3,event.ydata/1.e3))
    pass

def dic_raw_plots(dgdfilename):
    (Images,y0,x0,dy,dx,ny,nx,nimages,nloads,xbase,YPosns,StressPosns,ActualStressPosns)=load_dgd(dgdfilename)
    
    maxstress_idx=np.argmax(np.abs(StressPosns))
    for Yidx in range(YPosns.shape[0]):
        YPosn = YPosns[Yidx]*1e-3 # Posns is stored in mm

        # *** NOTE: Y axis as defined by motion stages and Y axis from images
        # are flipped in the recorded data. So here we flip the Y axis from the recorded data        
        #use_y0 = y0 - (ny//dic_scalefactor-1)*dy*dic_scalefactor
        use_y0 = y0 - (ny-1)*dy
        Yposvec=YPosn + use_y0 + np.arange(ny,dtype='d')*dy
        #Yposvec=YPosn - y0-np.arange(ny,dtype='d')*dy
        extent=np.array((x0-dx/2.0,x0+nx*dx-dx/2.0,Yposvec[0]-dy/2.0,Yposvec[-1]+dy/2.0,))*1e3
        fig=pl.figure()
        pl.imshow(Images[::-1,:,Yidx,maxstress_idx],origin='lower',extent=extent)
        pl.xlabel('X position')
        pl.ylabel('Y position')
        pl.title("Yidx=%d; YPosn=%f mm" % (Yidx,YPosn*1e3))

        fig.canvas.mpl_connect('button_press_event',dic_plot_click_handler)
        pass
    pass


def execute_one_dic(params):
    (idx2,input1,input2,ROI_buf,dic_scalefactor,dic_radius,n_threads,debug)=params
    
    (v_array,u_array,ROI_out_array) = correlate.correlate(input1,input2,ROI_buf,dic_scalefactor,dic_radius,n_threads=n_threads,debug=debug)

    return (idx2,v_array,u_array,ROI_out_array)

def execute_dic(dgdfilename,dgs_outfilename,dic_scalefactor,dic_radius,TipCoords1,TipCoords2,XRange,n_threads=multiprocessing.cpu_count(),processpool=None,debug=True):
    """Perform DIC on optical microscopy .dgd file. 
     dic_scalefactor and dic_radius parameters to ncorr, given in pixels
     TipCoords1 is an (x,y) tuple indicating the coordinates of the tip with a 
     lower value of y, in meters
     TipCoords2 is an (x,y) tuple indicating the coordinates of the tip with a
     larger value of y, in meters
     XRange is a (x1,x2) tuple indicating the lower and upper bounds of the region of 
     of interest in x, in meters

"""
    

    
    (Images,y0,x0,dy,dx,ny,nx,nimages,nloads,xbase,YPosns,StressPosns,ActualStressPosns)=load_dgd(dgdfilename)

    #dgs_outfilename=os.path.splitext(dgdfilename)[0]+"_dic.dgs"
    

    print("Perform_dic: Got %d %dx%d images at %d loads" % (Images.shape[2],nx,ny,nloads))
    


    CrackCenterY=(TipCoords1[1]+TipCoords2[1])/2.0

    
    ROI=np.zeros((ny,nx),dtype=np.uint8,order='F')
    # Should use better process to locate crack and identify ROI
    #ROI[450:840,:]=1
    
    ROI_xminidx=np.where(xbase > XRange[0])[0][0]
    
    ROI_xmaxidx=np.where(xbase < XRange[1])[0][-1]

    ROI[:,ROI_xminidx:ROI_xmaxidx]=1
    #raise ValueError("Break")
    YRange=(YPosns*1e-3 > TipCoords1[1]) & (YPosns*1e-3-ny*dy < TipCoords2[1])
    YRangeSize=np.count_nonzero(YRange)

    dic_ny = ny//dic_scalefactor
    dic_nx = nx//dic_scalefactor

    dic_dx = dx*dic_scalefactor
    dic_dy = dy*dic_scalefactor
    
    u_disps=np.ones((dic_ny,dic_nx,nloads,nloads,YRangeSize),dtype='f',order='F') # matrix of Python Objects to store u displacements
    u_disps[...]=np.nan
    v_disps=np.ones((dic_ny,dic_nx,nloads,nloads,YRangeSize),dtype='f',order='F') # matrix of Python Objects to store v displacements
    v_disps[...]=np.nan
    ROI_out_arrays=np.ones((dic_ny,dic_nx,nloads,nloads,YRangeSize),dtype='f',order='F')
    ROI_out_arrays[...]=np.nan

    load1=np.zeros((nloads,nloads,YRangeSize),dtype='f',order='F')
    load1[...]=np.nan
    load2=np.zeros((nloads,nloads,YRangeSize),dtype='f',order='F')
    load2[...]=np.nan
    
    YRange_idxs = np.where(YRange)[0]
    Yposvecs=np.ones((dic_ny,YRangeSize),dtype='f',order='F')
    Yposvecs[...]=np.nan

    Yinivec = np.ones(YRangeSize,dtype='f')
    Yinivec[...]=np.nan
    
    for YCnt in range(YRange_idxs.shape[0]):
        #if YCnt != 1:
        #    continue
        Yidx = YRange_idxs[YCnt]
        YPosn = YPosns[Yidx]*1e-3 # Posns is stored in mm
        # *** NOTE: Y axis as defined by motion stages and Y axis from images
        # are flipped in the recorded data. So here we flip the Y axis from the recorded data

        use_y0 = y0 - (ny-1)*dy
        Yinivec[YCnt]=use_y0+YPosn
        Yposvec=Yinivec[YCnt] + np.arange(ny//dic_scalefactor,dtype='d')*dy*dic_scalefactor
        Yposvecs[:,YCnt]=Yposvec
        for idx1 in range(nloads):
            #if idx1 != 0:
            #    continue

            # build up correlate_params so we can use map() (potentially with multiprocessing)
            print("Y=%d/%d; idx1=%d/%d" % (YCnt,YRange_idxs.shape[0],idx1,nloads))

            correlate_params=[]
            
            for idx2 in range(idx1+1,nloads):
                #if idx2 != nloads-1:
                #    continue
                load1[idx1,idx2,YCnt]=ActualStressPosns[Yidx,idx1]
                load2[idx1,idx2,YCnt]=ActualStressPosns[Yidx,idx2]
                load1[idx2,idx1,YCnt]=ActualStressPosns[Yidx,idx2]
                load2[idx2,idx1,YCnt]=ActualStressPosns[Yidx,idx1]
                
                input1=np.asfortranarray(Images[::-1,:,Yidx,idx1].astype(np.float64))
                input2=np.asfortranarray(Images[::-1,:,Yidx,idx2].astype(np.float64))
                ROI_buf=np.asfortranarray(ROI[::-1,:].astype(np.uint8))

                correlate_params.append((idx2,input1,input2,ROI_buf,dic_scalefactor,dic_radius,n_threads,debug))
                pass

            if processpool is None:
                correlate_results=map(execute_one_dic,correlate_params)
                pass
            else:
                correlate_results=processpool.map(execute_one_dic,correlate_params)
                pass
            for (idx2,v_array,u_array,ROI_out_array) in correlate_results:
                
                u_disps[:,:,idx1,idx2,YCnt]=u_array*dx
                v_disps[:,:,idx1,idx2,YCnt]=v_array*dy
                
                u_disps[:,:,idx2,idx1,YCnt]=-u_array*dx
                v_disps[:,:,idx2,idx1,YCnt]=-v_array*dy

                ROI_out_arrays[:,:,idx1,idx2,YCnt] = ROI_out_array
                ROI_out_arrays[:,:,idx2,idx1,YCnt] = ROI_out_array                
                pass
            
            pass
        #break
        pass
    
    # Write output to a .dgs file


    outwfmdict=collections.OrderedDict()

    for YCnt in range(YRange_idxs.shape[0]):
        
        outwfmdict["u_disps%.3d" % (YCnt)]=dg.wfminfo()
        outwfmdict["u_disps%.3d" % (YCnt)].Name="u_disps%.3d" % (YCnt)
        outwfmdict["u_disps%.3d" % (YCnt)].data=u_disps[...,YCnt]
        outwfmdict["u_disps%.3d" % (YCnt)].dimlen=np.array(u_disps.shape[:-1])
        outwfmdict["u_disps%.3d" % (YCnt)].ndim=4
        dgm.AddMetaDatumWI(outwfmdict["u_disps%.3d" % (YCnt)],dgm.CreateMetaDatumStr("Coord1","Y Position"))
        dgm.AddMetaDatumWI(outwfmdict["u_disps%.3d" % (YCnt)],dgm.CreateMetaDatumStr("Units1","meters"))
        dgm.AddMetaDatumWI(outwfmdict["u_disps%.3d" % (YCnt)],dgm.CreateMetaDatumDbl("IniVal1",Yinivec[YCnt]))
        dgm.AddMetaDatumWI(outwfmdict["u_disps%.3d" % (YCnt)],dgm.CreateMetaDatumDbl("Step1",dy*dic_scalefactor))
        dgm.AddMetaDatumWI(outwfmdict["u_disps%.3d" % (YCnt)],dgm.CreateMetaDatumStr("Coord2","X Position"))
        dgm.AddMetaDatumWI(outwfmdict["u_disps%.3d" % (YCnt)],dgm.CreateMetaDatumStr("Units2","meters"))
        dgm.AddMetaDatumWI(outwfmdict["u_disps%.3d" % (YCnt)],dgm.CreateMetaDatumDbl("IniVal2",x0))
        dgm.AddMetaDatumWI(outwfmdict["u_disps%.3d" % (YCnt)],dgm.CreateMetaDatumDbl("Step2",dx*dic_scalefactor))
        dgm.AddMetaDatumWI(outwfmdict["u_disps%.3d" % (YCnt)],dgm.CreateMetaDatumStr("Coord3","Stress Level DIC input 1"))
        dgm.AddMetaDatumWI(outwfmdict["u_disps%.3d" % (YCnt)],dgm.CreateMetaDatumStr("Coord4","Stress Level DIC input 2"))
        dgm.AddMetaDatumWI(outwfmdict["u_disps%.3d" % (YCnt)],dgm.CreateMetaDatumStr("Coord5","Image index"))
        dgm.AddMetaDatumWI(outwfmdict["u_disps%.3d" % (YCnt)],dgm.CreateMetaDatumStr("Units5","Unitless"))
        dgm.AddMetaDatumWI(outwfmdict["u_disps%.3d" % (YCnt)],dgm.CreateMetaDatumStr("AmplUnits","meters"))
        dgm.AddMetaDatumWI(outwfmdict["u_disps%.3d" % (YCnt)],dgm.CreateMetaDatumStr("X Displacement","meters"))
        

        outwfmdict["v_disps%.3d" % (YCnt)]=dg.wfminfo()
        outwfmdict["v_disps%.3d" % (YCnt)].Name="v_disps%.3d" % (YCnt)
        outwfmdict["v_disps%.3d" % (YCnt)].data=v_disps[...,YCnt]
        outwfmdict["v_disps%.3d" % (YCnt)].dimlen=np.array(v_disps.shape[:-1])
        outwfmdict["v_disps%.3d" % (YCnt)].ndim=4
        dgm.AddMetaDatumWI(outwfmdict["v_disps%.3d" % (YCnt)],dgm.CreateMetaDatumStr("Coord1","Y Position"))
        dgm.AddMetaDatumWI(outwfmdict["v_disps%.3d" % (YCnt)],dgm.CreateMetaDatumStr("Units1","meters"))
        dgm.AddMetaDatumWI(outwfmdict["v_disps%.3d" % (YCnt)],dgm.CreateMetaDatumDbl("IniVal1",Yinivec[YCnt]))
        dgm.AddMetaDatumWI(outwfmdict["v_disps%.3d" % (YCnt)],dgm.CreateMetaDatumDbl("Step1",dy*dic_scalefactor))
        dgm.AddMetaDatumWI(outwfmdict["v_disps%.3d" % (YCnt)],dgm.CreateMetaDatumStr("Coord2","X Position"))
        dgm.AddMetaDatumWI(outwfmdict["v_disps%.3d" % (YCnt)],dgm.CreateMetaDatumStr("Units2","meters"))
        dgm.AddMetaDatumWI(outwfmdict["v_disps%.3d" % (YCnt)],dgm.CreateMetaDatumDbl("IniVal2",x0))
        dgm.AddMetaDatumWI(outwfmdict["v_disps%.3d" % (YCnt)],dgm.CreateMetaDatumDbl("Step2",dx*dic_scalefactor))
        dgm.AddMetaDatumWI(outwfmdict["v_disps%.3d" % (YCnt)],dgm.CreateMetaDatumStr("Coord3","Stress Level DIC input 1"))
        dgm.AddMetaDatumWI(outwfmdict["v_disps%.3d" % (YCnt)],dgm.CreateMetaDatumStr("Coord4","Stress Level DIC input 2"))
        dgm.AddMetaDatumWI(outwfmdict["v_disps%.3d" % (YCnt)],dgm.CreateMetaDatumStr("Coord5","Image index"))
        dgm.AddMetaDatumWI(outwfmdict["v_disps%.3d" % (YCnt)],dgm.CreateMetaDatumStr("Units5","Unitless"))
        dgm.AddMetaDatumWI(outwfmdict["v_disps%.3d" % (YCnt)],dgm.CreateMetaDatumStr("AmplUnits","meters"))
        dgm.AddMetaDatumWI(outwfmdict["v_disps%.3d" % (YCnt)],dgm.CreateMetaDatumStr("Y Displacement","meters"))
        

        outwfmdict["ROI_out%.3d" % (YCnt)]=dg.wfminfo()
        outwfmdict["ROI_out%.3d" % (YCnt)].Name="ROI_out%.3d" % (YCnt)
        outwfmdict["ROI_out%.3d" % (YCnt)].data=ROI_out_arrays[...,YCnt]
        outwfmdict["ROI_out%.3d" % (YCnt)].dimlen=np.array(ROI_out_arrays.shape[:-1])
        outwfmdict["ROI_out%.3d" % (YCnt)].ndim=4
        dgm.AddMetaDatumWI(outwfmdict["ROI_out%.3d" % (YCnt)],dgm.CreateMetaDatumStr("Coord1","Y Position"))
        dgm.AddMetaDatumWI(outwfmdict["ROI_out%.3d" % (YCnt)],dgm.CreateMetaDatumStr("Units1","meters"))
        dgm.AddMetaDatumWI(outwfmdict["ROI_out%.3d" % (YCnt)],dgm.CreateMetaDatumDbl("IniVal1",Yinivec[YCnt]))
        dgm.AddMetaDatumWI(outwfmdict["ROI_out%.3d" % (YCnt)],dgm.CreateMetaDatumDbl("Step1",dy*dic_scalefactor))
        dgm.AddMetaDatumWI(outwfmdict["ROI_out%.3d" % (YCnt)],dgm.CreateMetaDatumStr("Coord2","X Position"))
        dgm.AddMetaDatumWI(outwfmdict["ROI_out%.3d" % (YCnt)],dgm.CreateMetaDatumStr("Units2","meters"))
        dgm.AddMetaDatumWI(outwfmdict["ROI_out%.3d" % (YCnt)],dgm.CreateMetaDatumDbl("IniVal2",x0))
        dgm.AddMetaDatumWI(outwfmdict["ROI_out%.3d" % (YCnt)],dgm.CreateMetaDatumDbl("Step2",dx*dic_scalefactor))
        dgm.AddMetaDatumWI(outwfmdict["ROI_out%.3d" % (YCnt)],dgm.CreateMetaDatumStr("Coord3","Stress Level DIC input 1"))
        dgm.AddMetaDatumWI(outwfmdict["ROI_out%.3d" % (YCnt)],dgm.CreateMetaDatumStr("Coord4","Stress Level DIC input 2"))
        dgm.AddMetaDatumWI(outwfmdict["ROI_out%.3d" % (YCnt)],dgm.CreateMetaDatumStr("Coord5","Image index"))
        dgm.AddMetaDatumWI(outwfmdict["ROI_out%.3d" % (YCnt)],dgm.CreateMetaDatumStr("Units5","Unitless"))

        pass
    
    outwfmdict["Yposvecs"]=dg.wfminfo()
    outwfmdict["Yposvecs"].Name="Yposvecs"
    outwfmdict["Yposvecs"].data=Yposvecs
    outwfmdict["Yposvecs"].dimlen=np.array(Yposvecs.shape)
    outwfmdict["Yposvecs"].ndim=2
    dgm.AddMetaDatumWI(outwfmdict["Yposvecs"],dgm.CreateMetaDatumStr("Coord1","Y Position"))
    dgm.AddMetaDatumWI(outwfmdict["Yposvecs"],dgm.CreateMetaDatumStr("Units1","meters"))
    dgm.AddMetaDatumWI(outwfmdict["Yposvecs"],dgm.CreateMetaDatumDbl("IniVal1",y0))
    dgm.AddMetaDatumWI(outwfmdict["Yposvecs"],dgm.CreateMetaDatumDbl("Step1",dy*dic_scalefactor))
    dgm.AddMetaDatumWI(outwfmdict["Yposvecs"],dgm.CreateMetaDatumStr("Coord2","Image index"))
    dgm.AddMetaDatumWI(outwfmdict["Yposvecs"],dgm.CreateMetaDatumStr("Units2","Unitless"))
    dgm.AddMetaDatumWI(outwfmdict["Yposvecs"],dgm.CreateMetaDatumStr("AmplCoord","Y shiftedscaled"))
    dgm.AddMetaDatumWI(outwfmdict["Yposvecs"],dgm.CreateMetaDatumStr("AmplUnits","meters"))


    outwfmdict["Yinivec"]=dg.wfminfo()
    outwfmdict["Yinivec"].Name="Yinivec"
    outwfmdict["Yinivec"].data=Yinivec
    outwfmdict["Yinivec"].dimlen=np.array(Yinivec.shape)
    outwfmdict["Yinivec"].ndim=1
    dgm.AddMetaDatumWI(outwfmdict["Yinivec"],dgm.CreateMetaDatumStr("Coord1","Image index"))
    dgm.AddMetaDatumWI(outwfmdict["Yinivec"],dgm.CreateMetaDatumStr("Units1","Unitless"))
    dgm.AddMetaDatumWI(outwfmdict["Yinivec"],dgm.CreateMetaDatumStr("AmplCoord","Y initial"))
    dgm.AddMetaDatumWI(outwfmdict["Yinivec"],dgm.CreateMetaDatumStr("AmplUnits","meters"))


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
    dgm.AddMetaDatumL(outmetadata,dgm.CreateMetaDatumDbl("CrackCenterY",CrackCenterY))
    dgm.AddMetaDatumL(outmetadata,dgm.CreateMetaDatumDbl("TipCoords1X",TipCoords1[0]))
    dgm.AddMetaDatumL(outmetadata,dgm.CreateMetaDatumDbl("TipCoords1Y",TipCoords1[1]))
    dgm.AddMetaDatumL(outmetadata,dgm.CreateMetaDatumDbl("TipCoords2X",TipCoords2[0]))
    dgm.AddMetaDatumL(outmetadata,dgm.CreateMetaDatumDbl("TipCoords2Y",TipCoords2[1]))

    dgm.AddMetaDatumL(outmetadata,dgm.CreateMetaDatumInt("NumImages",YRange_idxs.shape[0]))

    dgm.AddMetaDatumL(outmetadata,dgm.CreateMetaDatumInt("ROI_xminidx",ROI_xminidx))
    dgm.AddMetaDatumL(outmetadata,dgm.CreateMetaDatumInt("ROI_xmaxidx",ROI_xmaxidx))

    dgm.AddMetaDatumL(outmetadata,dgm.CreateMetaDatumInt("ROI_dic_xminidx",ROI_xminidx//dic_scalefactor))
    dgm.AddMetaDatumL(outmetadata,dgm.CreateMetaDatumInt("ROI_dic_xmaxidx",ROI_xmaxidx//dic_scalefactor))


    dgfh=dgf.creat(dgs_outfilename)
    dgf.writesnapshot(dgfh,outmetadata,outwfmdict)
    dgf.close(dgfh)

    return (outwfmdict,outmetadata,u_disps,v_disps,ROI_out_arrays,Yposvecs,Yinivec,CrackCenterY,dic_dx,dic_dy)
    
