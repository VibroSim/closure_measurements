import sys
import numpy as np

#from matplotlib import pyplot as pl

import scipy
import scipy.integrate

import dg_file as dgf
import dg_metadata as dgm

from . import dic_ctod
from . import initial_fit
from . import full_model
from . import full_model_accel

def load_dgs(dgsfilename):
    (metadatadict,wfmdict)=dgf.loadsnapshot(dgsfilename)

    YRangeSize=wfmdict["Yinivec"].data.shape[0]
    nloads = wfmdict["load1"].data.shape[0]

    dic_ny=wfmdict["u_disps000"].data.shape[0]
    dic_nx=wfmdict["u_disps000"].data.shape[1]

    dic_dx=dgm.GetMetaDatumWIDbl(wfmdict["u_disps000"],"Step2",np.nan)
    dic_dy=dgm.GetMetaDatumWIDbl(wfmdict["u_disps000"],"Step1",np.nan)
    
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

    for YCnt in range(YRangeSize):
        u_disps[:,:,:,:,YCnt] = wfmdict["u_disps%.3d" % (YCnt)].data
        v_disps[:,:,:,:,YCnt] = wfmdict["v_disps%.3d" % (YCnt)].data
        ROI_out_arrays[:,:,:,:,YCnt] = wfmdict["ROI_out%.3d" % (YCnt)].data
        pass

    load1=wfmdict["load1"].data
    load2=wfmdict["load2"].data
    Yposvecs=wfmdict["Yposvecs"].data
    
    Yinivec=wfmdict["Yinivec"].data

    CrackCenterY=dgm.GetMetaDatumLDbl(metadatadict,"CrackCenterY",np.nan)
    TipCoords1=(dgm.GetMetaDatumLDbl(metadatadict,"TipCoords1X",np.nan),
                dgm.GetMetaDatumLDbl(metadatadict,"TipCoords1Y",np.nan))
    TipCoords2=(dgm.GetMetaDatumLDbl(metadatadict,"TipCoords2X",np.nan),
                dgm.GetMetaDatumLDbl(metadatadict,"TipCoords2Y",np.nan))
    ROI_dic_xminidx=dgm.GetMetaDatumLInt(metadatadict,"ROI_dic_xminidx",0)
    ROI_dic_xmaxidx=dgm.GetMetaDatumLInt(metadatadict,"ROI_dic_xmaxidx",0)
    

    return (dic_dy,dic_dx,dic_ny,dic_nx,YRangeSize,nloads,Yinivec,Yposvecs,load1,load2,u_disps,v_disps,ROI_out_arrays,CrackCenterY,TipCoords1,TipCoords2,ROI_dic_xminidx,ROI_dic_xmaxidx)



def Calc_CTODs(dic_ny,nloads,YRangeSize,Yposvecs,u_disps,ROI_out_arrays,ROI_dic_xminidx,ROI_dic_xmaxidx,dic_span,dic_window):
    
    CTODs=np.zeros((dic_ny,nloads,nloads,YRangeSize),dtype='d')
    for YCnt in range(YRangeSize):

        Yposvec=Yposvecs[:,YCnt]
        for idx1 in range(nloads):
            #if idx1 != 0:
            #   ccccccjfh continue
            for idx2 in range(idx1+1,nloads):

                #if idx2 != nloads-1:
                #    continue
                (CTOD,top_disp,bot_disp) = dic_ctod.dic_ctod(u_disps[:,:,idx1,idx2,YCnt],dic_span,dic_window,ROI_out_arrays[:,:,idx1,idx2,YCnt],ROI_dic_xminidx,ROI_dic_xmaxidx)
                CTODs[:,idx1,idx2,YCnt]=CTOD
                CTODs[:,idx2,idx1,YCnt]=-CTOD
                pass
            pass
        pass
    return CTODs
    
    
def CalcInitialModel(nloads,CTODs,load1,load2,Yposvecs,CrackCenterY,side,nominal_length=2e-3,nominal_modulus=100.0e9,nominal_stress=50e6,doplots=False):
    """ side=1 (smaller Y values) or side=2 (larger Y values)"""
    
    InitialModels=np.zeros((nloads,nloads),dtype='O')
    InitialModels[:,:]=None
    InitialCoeffs=np.zeros((2,nloads,nloads),dtype='d') # Zeroth axis is c5, yt
    InitialCoeffs[:,:,:]=np.nan
    
    Error=np.zeros((nloads,nloads),dtype='d')
    Error[:,:] = np.nan

    npoints = np.zeros((nloads,nloads),dtype=np.uint32)

    YPositions = np.zeros((nloads,nloads),dtype='O')
    YPositions[...]=None
    CTODValues = np.zeros((nloads,nloads),dtype='O')
    CTODValues[...]=None

    for idx1 in range(nloads):
        #if idx1 != 0:
        #    continue
        for idx2 in range(idx1+1,nloads):
            
            #if idx2 != nloads-1:
            #    continue
            
            YPositions[idx1,idx2]=np.array([],dtype='f')
            CTODValues[idx1,idx2]=np.array([],dtype='d')
            for YCnt in range(Yposvecs.shape[1]):
                #if YCnt != 1:
                #    continue
                Yposvec=Yposvecs[:,YCnt]
                CTOD=CTODs[:,idx1,idx2,YCnt]
                
                if side==1:
                    valid_locations = (Yposvec < CrackCenterY) & (~np.isnan(CTOD))
                    pass
                else:
                    assert(side==2)
                    valid_locations = (Yposvec > CrackCenterY) & (~np.isnan(CTOD))
                    pass
                
                YPositions[idx1,idx2]=np.concatenate((YPositions[idx1,idx2],Yposvec[valid_locations]))
                CTODValues[idx1,idx2]=np.concatenate((CTODValues[idx1,idx2],CTOD[valid_locations]))
                pass


            if YPositions[idx1,idx2].shape[0] == 0:
                # No Data!
                continue

            y0=(np.sqrt(nominal_length)/nominal_modulus,np.mean(YPositions[idx1,idx2]))
            #y0=(1.0e-13,np.mean(YPositions[idx1,idx2]))
            (c5,yt)=initial_fit.fit_initial_model(y0,YPositions[idx1,idx2],load1[idx1,idx2,YCnt],load2[idx1,idx2,YCnt],side,CTODValues[idx1,idx2],nominal_length,nominal_modulus,nominal_stress)
            print("side=%d; yt=%f" % (side,yt))
            InitialModels[idx1,idx2]=initial_fit.initial_model((c5,yt),YPositions[idx1,idx2],load1[idx1,idx2,YCnt],load2[idx1,idx2,YCnt],side)
            #InitialModels[idx2,idx1]=-initial_fit.initial_model((c5,yt),YPositions,load1[idx1,idx2,YCnt],load2[idx1,idx2,YCnt],side)
            InitialCoeffs[:,idx1,idx2]=(c5,yt)
            npoints[idx1,idx2]=YPositions[idx1,idx2].shape[0]
            Error[idx1,idx2]=np.sum((CTODValues[idx1,idx2]-InitialModels[idx1,idx2])**2.0)/npoints[idx1,idx2]
            #import pdb
            #pdb.set_trace()
            junk = initial_fit.initial_model((c5,yt),YPositions[idx1,idx2],load1[idx1,idx2,YCnt],load2[idx1,idx2,YCnt],side)
            if doplots:
                # Do random sub-percentage
                if np.random.rand() < .05:
                    
                    from matplotlib import pyplot as pl
                    
                    pl.figure()
                    YPositionsSort=np.argsort(YPositions[idx1,idx2])
                    YPositionsSorted=YPositions[idx1,idx2][YPositionsSort]
                    #CTODValuesSorted=CTODValues[idx1,idx2][YPositionsSort]
                    InitialModelValuesSorted=InitialModels[idx1,idx2][YPositionsSort]
                    pl.plot(YPositions[idx1,idx2]*1e3,CTODValues[idx1,idx2]*1e6,'.',
                            YPositionsSorted*1e3,InitialModelValuesSorted*1e6,'-')
                    pl.xlabel('Y (mm)')
                    pl.ylabel('CTOD and initial model (um)')
                    pl.title('First end of crack: Load1 = %f MPa; load2 = %f MPa' % (load1[idx1,idx2,0]/1.e6,load2[idx1,idx2,0]/1.e6))
                    pl.grid()
                    pass
                pass
            
            pass
        pass


    return (InitialModels,
            InitialCoeffs,
            Error,
            npoints,
            YPositions,
            CTODValues)



def EvalEffectiveTip(minload,maxload,full_model_params,sigma):
    # create (t,c,k) for scipy splev
    t=np.array([minload]*4 + [maxload]*4,dtype='d')  # four copies of minload followed by four copies of maxload

    splinecoeff=full_model_params[:4]
    #c5=params[4]

    c=np.concatenate((splinecoeff,[0.0]*4))
    k=3
    tck = (t,c,k)
    
    yt = scipy.interpolate.splev(sigma,tck)
    return yt


def CalcFullModel(load1,load2,InitialCoeffs,Error,npoints,YPositions,CTODValues,side,doplots=True,opencl_ctx=None,opencl_dev=None):
    # Our model is dCOD/dsigma = C5*sqrt(y-yt)u(y-yt) where u(y) is the unit step
    # This integrates to:
    #  COD2-COD1 = integral_sigma1^sigma2 C5*sqrt(y-yt)*u(y-yt) dsigma
    #   where yt is a function of sigma, because the tip location shifts
    #   as the crack opens and closes
    
    # From Suresh, eq. 9.45 @ theta=pi, uy = (K1/(2.0*E))*sqrt((r/(2*pi)))*(1+nu)*(2kappa+1 +1)
    # where kappa = (3-nu)/(1+nu) for plane stress or
    # kappa=(3-4nu) for plane strain... where nu is Poisson's ratio

    # COD=2uy = ((KI/E)/sqrt(2*pi))*(1+nu)*(2kappa+2)  *sqrt(y-yt) 
   # where KI = sigma*sqrt(pi*a)
    # Since this is proportional to sigma, for a fixed length crack
    # Our C5 is expected to equal sqrt(pi*a)*((1/E)/sqrt(2*pi))*(1+nu)*(2kappa+2)
    # or C5=(sqrt(2a)/(E))*(1+nu)*(kappa+1)
    #
    # In our "initial" model we use data from a single pair of loads,
    # assuming the crack tip position yt is fixed. The corresponding
    # delta sigma is the difference of the loads. The corresponding
    # sigma is the average of the two loads. 
    #
    # In our full model we use all data from all load pairs instead.
    # The initial model seeds a spline fit for yt(sigma)

    # Use a first cut C5 estimate to filter out any coefficients that
    # optimized to zero for whatever reason
    min_c5 = np.sqrt(2*50e-6)/1000e9

    (idx1grid,idx2grid)=np.meshgrid(np.arange(InitialCoeffs.shape[1]),np.arange(InitialCoeffs.shape[2]),indexing="ij")
    
    
    # Unwrap the error and yt coefficients
    valid = (~np.isnan(InitialCoeffs[1,:,:].ravel())) & (npoints.ravel() > 20) &  (InitialCoeffs[0,:,:].ravel() >= min_c5)
    
    Error_unwrapped=Error.ravel()[valid]
    c5_unwrapped = InitialCoeffs[0,:,:].ravel()[valid]
    yt_unwrapped = InitialCoeffs[1,:,:].ravel()[valid]


    idx1_unwrapped = idx1grid[:,:].ravel()[valid]
    idx2_unwrapped = idx2grid[:,:].ravel()[valid]
    
    avg_load=(load1+load2).mean(axis=2)/2.0
    avg_load_unwrapped=avg_load.ravel()[valid]
    
    # Now sort them so lowest error comes first
    errsort=np.argsort(Error_unwrapped)

    # Do fit to first 50
    yt_vals = yt_unwrapped[errsort[:50]]
    avg_load_vals=avg_load_unwrapped[errsort[:50]]
    c5_vals = c5_unwrapped[errsort[:50]]
    

    avg_load_vals_sort=np.argsort(avg_load_vals)
    avg_load_vals_sorted=avg_load_vals[avg_load_vals_sort]
    yt_vals_sorted = yt_vals[avg_load_vals_sort]
    
    minload=np.min(load1[~np.isnan(load1)].ravel())
    maxload=np.max(load1[~np.isnan(load1)].ravel())
    
    # Use scipy.interpolate.splrep to do fit.
    # task=-1 means we supply interior knots (there aren't any)
    # k=3 means 3rd order
    (t,c,k) = scipy.interpolate.splrep(avg_load_vals_sorted,yt_vals_sorted,xb=minload,xe=maxload,k=3,task=-1,t=np.array([],dtype='f'))
    assert((t[:4]==minload).all())
    assert((t[4:]==maxload).all())

    assert((c[4:]==0).all()) # This spline only has 4 coefficients. For some reason splrep returns four more that are all zero
    assert(k==3)
    seed_param=(c[0],c[1],c[2],c[3],c5)
    
    # Plot diagnostics
    if doplots:

        sigmarange=np.linspace(minload,maxload,150)
        fittedvals=EvalEffectiveTip(minload,maxload,seed_param,sigmarange)

        from matplotlib import pyplot as pl

        pl.figure()
        pl.plot(yt_unwrapped*1e3,avg_load_unwrapped/1e6,'x',
                yt_vals*1e3,avg_load_vals/1e6,'o',
                fittedvals*1e3,sigmarange/1e6,'-')
        pl.ylabel('Load (MPa)')
        pl.xlabel('Tip position (mm)')
        pl.legend(('All DIC fit data','Best 50','Fit to best 50'))
        pl.title('yt')
        pl.grid()

        
        fig=pl.figure()
        pl.plot(yt_unwrapped*1e3,avg_load_unwrapped/1e6,'x',picker=5)
                #yt_vals*1e3,avg_load_vals/1e6,'o',
                #fittedvals*1e3,sigmarange/1e6,'-')
        pl.ylabel('Load (MPa)')
        pl.xlabel('Tip position (mm)')
        pl.grid()
        pl.title('yt (pickable)')
        def dicfitpick(event):
            thisline=event.artist
            xdata=thisline.get_xdata()
            ydata=thisline.get_ydata()
            indices = event.ind
            print("got indices: %s" % (str(indices)))

            for index in indices:
                pl.figure()

                idx1=idx1_unwrapped[index]
                idx2=idx2_unwrapped[index]
                
                YPositionsSort=np.argsort(YPositions[idx1,idx2])
                YPositionsSorted=YPositions[idx1,idx2][YPositionsSort]
                #CTODValuesSorted=CTODValues[idx1,idx2][YPositionsSort]
                InitialModelValuesSorted=InitialModels[idx1,idx2][YPositionsSort]
                pl.plot(YPositions[idx1,idx2]*1e3,CTODValues[idx1,idx2]*1e6,'.',
                        YPositionsSorted*1e3,InitialModelValuesSorted*1e6,'-')
                pl.xlabel('Y (mm)')
                pl.ylabel('CTOD and initial model (um)')
                pl.title('First end of crack: Load1 = %f MPa; load2 = %f MPa' % (load1[idx1,idx2,0]/1.e6,load2[idx1,idx2,0]/1.e6))
                pl.grid()
                
            
            pass
        fig.canvas.mpl_connect('pick_event',dicfitpick)
        
        pl.figure()
        pl.plot(avg_load_unwrapped/1e6,c5_unwrapped,'x',
                avg_load_vals/1e6,c5_vals,'o')
        pl.xlabel('Load (MPa)')
        pl.title('c5')
        pl.grid()
        pass


    full_model_residual_plot=None

    if doplots:
        from matplotlib import pyplot as pl
        full_model_residual_plot=pl.figure()
        pass
    
    # Perform model fit

    """  full model calculation temporarily commented out
    if opencl_ctx is None:
        full_model_residual=full_model.full_model_residual
        args=(YPositions,CTODValues,np.mean(load1,axis=2),np.mean(load2,axis=2),minload,maxload,side,full_model_residual_plot)
        pass
    else:
        full_model_residual=full_model_accel.full_model_residual_accel
        args=(YPositions,CTODValues,np.mean(load1,axis=2),np.mean(load2,axis=2),minload,maxload,side,full_model_residual_plot,opencl_ctx,opencl_dev)

    full_model_result = scipy.optimize.minimize(full_model_residual,seed_param,args=args,method="nelder-mead",tol=1e-17)
    full_model_params=full_model_result.x
    """
    full_model_result=None
    full_model_params = seed_param

    return (minload,maxload,full_model_params,full_model_result)
                     
