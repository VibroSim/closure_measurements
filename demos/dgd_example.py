import sys
import numpy as np

from matplotlib import pyplot as pl

import scipy
import scipy.optimize

import dg_file as dgf
import dg_dgdread

import pyximport
pyximport.install()

import correlate
import dic_ctod
import initial_fit

#dgdfilename = sys.argv[1]
dgdfilename = "/tmp/C14-UTCB-004F_tortuosity_2014-09-17_collect_optical_data-0027.dgd"
scalefactor=5
r=20


# ***!!!! NOTE: This code is obsolete -- need to swap x and y labels and u and v displacements

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

assert(PosnUnits["Y"]=="mm")

Images=ReshapedWfmData["GEV"].data
# Images is (in one case) 1200x1920x19x2
# axes are Y,X,...
# pl.imshow(Images[:,:,2,0].T,vmin=0.3,vmax=0.5)

ny=Images.shape[0]
nx=Images.shape[1]
nloads=Images.shape[3]

xbase=x0+np.arange(nx,dtype='d')*dx

dic_span=20
dic_window=3

do_raw_plots=False
maxstress_idx=np.argmax(np.abs(Posns["Stress"]))
if do_raw_plots:
    for Yidx in range(Posns["Y"].shape[0]):
        YPosn = Posns["Y"][Yidx]*1e-3 # Posns is stored in mm
        Yposvec=YPosn - np.arange(ny,dtype='d')*dy
        extent=np.array((x0-dx/2.0,x0+nx*dx-dx/2.0,Yposvec[-1]-dy/2.0,Yposvec[0]+dy/2.0,))*1e3
        pl.figure()
        pl.imshow(Images[:,:,Yidx,maxstress_idx],origin='upper',extent=extent)
        pl.xlabel('X position')
        pl.ylabel('Y position')
        pl.title("Yidx=%d; YPosn=%f mm" % (Yidx,YPosn*1e3))
        pass
    pass

TipCoords1=(.313e-3,3.49e-3) # should have smaller value of y
TipCoords2=(.335e-3,7.20e-3) # Should have larger value of y
XRange=(.15e-3,.8e-3)


CrackCenterY=(TipCoords1[1]+TipCoords2[1])/2.0


ROI=np.zeros((ny,nx),dtype=np.uint8,order='F')
# Should use better process to locate crack and identify ROI
#ROI[450:840,:]=1

ROI_xminidx=np.where(xbase > XRange[0])[0][0]

ROI_xmaxidx=np.where(xbase < XRange[1])[0][-1]

ROI[:,ROI_xminidx:ROI_xmaxidx]=1
#raise ValueError("Break")
YRange=(Posns["Y"]*1e-3 > TipCoords1[1]) & (Posns["Y"]*1e-3-ny*dy < TipCoords2[1])
YRangeSize=np.count_nonzero(YRange)


u_disps=np.zeros((YRangeSize,nloads,nloads),dtype='O') # matrix of Python Objects to store u displacements
v_disps=np.zeros((YRangeSize,nloads,nloads),dtype='O') # matrix of Python Objects to store v displacements


YRange_idxs = np.where(YRange)[0]
Yposvecs=np.zeros(YRangeSize,dtype='O')

for YCnt in range(YRange_idxs.shape[0]):
    #if YCnt != 1:
    #    continue
    Yidx = YRange_idxs[YCnt]
    YPosn = Posns["Y"][Yidx]*1e-3 # Posns is stored in mm
    Yposvec=YPosn - np.arange(ny//scalefactor,dtype='d')*dy*scalefactor
    Yposvecs[YCnt]=Yposvec
    for idx1 in range(nloads):
        #if idx1 != 0:
        #    continue
        load1=ActualPosns["Stress"][Yidx,idx1]
        for idx2 in range(idx1+1,nloads):
            print("Y=%d/%d; idx1=%d/%d; idx2=%d/%d" % (YCnt,YRange_idxs.shape[0],idx1,nloads,idx2,nloads))
            #if idx2 != nloads-1:
            #    continue
            load2=ActualPosns["Stress"][Yidx,idx2]
            input1=np.asfortranarray(Images[:,:,Yidx,idx1].astype(np.float64))
            input2=np.asfortranarray(Images[:,:,Yidx,idx2].astype(np.float64))
            (v_array,u_array,ROI_out_array) = correlate.correlate(input1,input2,ROI,scalefactor,r,debug=True)

                        
            u_disps[YCnt,idx1,idx2]=u_array*dx
            v_disps[YCnt,idx1,idx2]=v_array*dy

            u_disps[YCnt,idx2,idx1]=-u_array*dx
            v_disps[YCnt,idx2,idx1]=-v_array*dy


            pass
        pass
    #break
    pass


CTODs=np.zeros((YRangeSize,nloads,nloads),dtype='O')


for YCnt in range(YRange_idxs.shape[0]):
    #if YCnt != 1:
    #    continue
    Yidx = YRange_idxs[YCnt]
    YPosn = Posns["Y"][Yidx]*1e-3 # Posns is stored in mm
    #Yposvec=YPosn - np.arange(ny//scalefactor,dtype='d')*dy*scalefactor
    Yposvec=Yposvecs[YCnt]
    
    for idx1 in range(nloads):
        #if idx1 != 0:
        #    continue
        load1=ActualPosns["Stress"][Yidx,idx1]
        for idx2 in range(idx1+1,nloads):

            #if idx2 != nloads-1:
            #    continue
            load2=ActualPosns["Stress"][Yidx,idx2]
            (CTOD,top_disp,bot_disp) = dic_ctod.dic_ctod(u_disps[YCnt,idx1,idx2],dic_span,dic_window,ROI_out_array,ROI_xminidx//scalefactor,ROI_xmaxidx//scalefactor)
            CTODs[YCnt,idx1,idx2]=CTOD
            CTODs[YCnt,idx2,idx1]=-CTOD


            pass
        pass
    #break
    pass




# First end of crack
#InitialModels=np.zeros((YRangeSize,nloads,nloads),dtype='O')


InitialModels_side1=np.zeros((nloads,nloads),dtype='O')

for idx1 in range(nloads):
    #if idx1 != 0:
    #    continue
    load1=ActualPosns["Stress"][Yidx,idx1]
    for idx2 in range(idx1+1,nloads):
        
        #if idx2 != nloads-1:
        #    continue
        load2=ActualPosns["Stress"][Yidx,idx2]

        YPositions=np.array([],dtype='d')
        CTODValues=np.array([],dtype='d')
        for YCnt in range(YRange_idxs.shape[0]):
            #if YCnt != 1:
            #    continue
            Yidx = YRange_idxs[YCnt]
            YPosn = Posns["Y"][Yidx]*1e-3 # Posns is stored in mm
            #Yposvec=YPosn - np.arange(ny//scalefactor,dtype='d')*dy*scalefactor
            Yposvec=Yposvecs[YCnt]
            CTOD=CTODs[YCnt,idx1,idx2]


            valid_locations = (Yposvec < CrackCenterY) & (~np.isnan(CTOD))
            
            YPositions=np.concatenate((YPositions,Yposvec[valid_locations]))
            CTODValues=np.concatenate((CTODValues,CTOD[valid_locations]))
            pass
        
        y0=(1.0e-13,np.mean(YPositions))
        (c5,yt)=initial_fit.fit_initial_model(y0,YPositions,load1,load2,1,CTODValues)
        InitialModels_side1[idx1,idx2]=initial_fit.initial_model((c5,yt),YPositions,load1,load2,1)
        InitialModels_side1[idx2,idx1]=-initial_fit.initial_model((c5,yt),YPositions,load1,load2,1)
            
        pl.figure()
        pl.plot(YPositions*1e3,CTODValues*1e6,'-',
                YPositions*1e3,InitialModels_side1[idx1,idx2]*1e6,'-')
        pl.xlabel('Y (mm)')
        pl.ylabel('CTOD and initial model (um)')
        pl.title('First end of crack: Load1 = %f MPa; load2 = %f MPa' % (load1/1.e6,load2/1.e6))
        pl.grid()
        pass
    pass




# Second end of crack
#InitialModels=np.zeros((YRangeSize,nloads,nloads),dtype='O')


InitialModels_side2=np.zeros((nloads,nloads),dtype='O')

for idx1 in range(nloads):
    #if idx1 != 0:
    #    continue
    load1=ActualPosns["Stress"][Yidx,idx1]
    for idx2 in range(idx1+1,nloads):
        
        #if idx2 != nloads-1:
        #    continue
        load2=ActualPosns["Stress"][Yidx,idx2]

        YPositions=np.array([],dtype='d')
        CTODValues=np.array([],dtype='d')
        for YCnt in range(YRange_idxs.shape[0]):
            #if YCnt != 1:
            #    continue
            Yidx = YRange_idxs[YCnt]
            YPosn = Posns["Y"][Yidx]*1e-3 # Posns is stored in mm
            #Yposvec=YPosn - np.arange(ny//scalefactor,dtype='d')*dy*scalefactor
            Yposvec=Yposvecs[YCnt]
            CTOD=CTODs[YCnt,idx1,idx2]
            
            
            valid_locations = (Yposvec >= CrackCenterY) & (~np.isnan(CTOD))
            
            YPositions=np.concatenate((YPositions,Yposvec[valid_locations]))
            CTODValues=np.concatenate((CTODValues,CTOD[valid_locations]))
            pass
        
        y0=(1.0e-13,np.mean(YPositions))
        (c5,yt)=initial_fit.fit_initial_model(y0,YPositions,load1,load2,2,CTODValues)
        InitialModels_side2[idx1,idx2]=initial_fit.initial_model((c5,yt),YPositions,load1,load2,2)
        InitialModels_side2[idx2,idx1]=-initial_fit.initial_model((c5,yt),YPositions,load1,load2,2)
            
        pl.figure()
        pl.plot(YPositions*1e3,CTODValues*1e6,'-',
                YPositions*1e3,InitialModels_side2[idx1,idx2]*1e6,'-')
        pl.xlabel('Y (mm)')
        pl.ylabel('CTOD and initial model (um)')
        pl.title('First end of crack: Load1 = %f MPa; load2 = %f MPa' % (load1/1.e6,load2/1.e6))
        pl.grid()
        pass
    pass

