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

    XRangeSize=wfmdict["Xinivec"].data.shape[0]
    nloads = wfmdict["load1"].data.shape[0]

    dic_nx=wfmdict["v_disps000"].data.shape[0]
    dic_ny=wfmdict["v_disps000"].data.shape[1]

    dic_dx=dgm.GetMetaDatumWIDbl(wfmdict["v_disps000"],"Step1",np.nan)
    dic_dy=dgm.GetMetaDatumWIDbl(wfmdict["v_disps000"],"Step2",np.nan)
    
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

    for XCnt in range(XRangeSize):
        u_disps[:,:,:,:,XCnt] = wfmdict["u_disps%.3d" % (XCnt)].data
        v_disps[:,:,:,:,XCnt] = wfmdict["v_disps%.3d" % (XCnt)].data
        ROI_out_arrays[:,:,:,:,XCnt] = wfmdict["ROI_out%.3d" % (XCnt)].data
        pass

    load1=wfmdict["load1"].data
    load2=wfmdict["load2"].data
    Xposvecs=wfmdict["Xposvecs"].data
    
    Xinivec=wfmdict["Xinivec"].data

    CrackCenterX=dgm.GetMetaDatumLDbl(metadatadict,"CrackCenterX",np.nan)
    TipCoords1=(dgm.GetMetaDatumLDbl(metadatadict,"TipCoords1X",np.nan),
                dgm.GetMetaDatumLDbl(metadatadict,"TipCoords1Y",np.nan))
    TipCoords2=(dgm.GetMetaDatumLDbl(metadatadict,"TipCoords2X",np.nan),
                dgm.GetMetaDatumLDbl(metadatadict,"TipCoords2Y",np.nan))
    ROI_dic_yminidx=dgm.GetMetaDatumLInt(metadatadict,"ROI_dic_yminidx",0)
    ROI_dic_ymaxidx=dgm.GetMetaDatumLInt(metadatadict,"ROI_dic_ymaxidx",0)


    relshift_firstimg_lowerleft_corner_x_ref = None
    relshift_firstimg_lowerleft_corner_x_diff = None
    relshift_firstimg_lowerleft_corner_y_ref = None
    relshift_firstimg_lowerleft_corner_y_diff = None

    if "relshift_firstimg_lowerleft_corner_x_ref" in wfmdict:
        relshift_firstimg_lowerleft_corner_x_ref = wfmdict["relshift_firstimg_lowerleft_corner_x_ref"].data
        pass

    if "relshift_firstimg_lowerleft_corner_x_diff" in wfmdict:
        relshift_firstimg_lowerleft_corner_x_diff = wfmdict["relshift_firstimg_lowerleft_corner_x_diff"].data
        pass

    
    
    if "relshift_firstimg_lowerleft_corner_y_ref" in wfmdict:
        relshift_firstimg_lowerleft_corner_y_ref = wfmdict["relshift_firstimg_lowerleft_corner_y_ref"].data
        pass

    if "relshift_firstimg_lowerleft_corner_y_diff" in wfmdict:
        relshift_firstimg_lowerleft_corner_y_diff = wfmdict["relshift_firstimg_lowerleft_corner_y_diff"].data
        pass

    return (dic_dx,dic_dy,
            dic_nx,dic_ny,
            XRangeSize,
            nloads,
            Xinivec,Xposvecs,
            load1,load2,u_disps,v_disps,
            ROI_out_arrays,
            CrackCenterX,TipCoords1,TipCoords2,
            ROI_dic_yminidx,ROI_dic_ymaxidx,
            relshift_firstimg_lowerleft_corner_x_ref,
            relshift_firstimg_lowerleft_corner_x_diff,
            relshift_firstimg_lowerleft_corner_y_ref,
            relshift_firstimg_lowerleft_corner_y_diff)



def Calc_CTODs(dic_nx,nloads,XRangeSize,Xposvecs,v_disps,ROI_out_arrays,ROI_dic_yminidx,ROI_dic_ymaxidx,dic_span,dic_window):
    
    CTODs=np.zeros((dic_nx,nloads,nloads,XRangeSize),dtype='d')
    for XCnt in range(XRangeSize):

        Xposvec=Xposvecs[:,XCnt]
        for idx1 in range(nloads):
            #if idx1 != 0:
            #   continue
            for idx2 in range(idx1+1,nloads):

                #if idx2 != nloads-1:
                #    continue
                (CTOD,top_disp,bot_disp) = dic_ctod.dic_ctod(v_disps[:,:,idx1,idx2,XCnt],dic_span,dic_window,ROI_out_arrays[:,:,idx1,idx2,XCnt],ROI_dic_yminidx,ROI_dic_ymaxidx)
                CTODs[:,idx1,idx2,XCnt]=CTOD
                #CTODs[:,idx2,idx1,XCnt]=-CTOD
                pass
            pass
        pass
    return CTODs
    

def TestRegistration(nloads,Xposvecs,u_disps,v_disps,
                     ROI_out_arrays,
                     relshift_firstimg_lowerleft_corner_x_ref=None,
                     relshift_firstimg_lowerleft_corner_x_diff=None,
                     relshift_firstimg_lowerleft_corner_y_ref=None,
                     relshift_firstimg_lowerleft_corner_y_diff=None):
    for idx1 in range(nloads):
        for idx2 in range(idx1+1,nloads):
            
            for XCnt in range(Xposvecs.shape[1]):
                Xposvec=Xposvecs[:,XCnt]

                u_values = u_disps[:,:,idx1,idx2,XCnt][ROI_out_arrays[:,:,idx1,idx2,XCnt]] # only use values where ROI==1
                meanu = np.mean(u_values.ravel())
                v_values = v_disps[:,:,idx1,idx2,XCnt][ROI_out_arrays[:,:,idx1,idx2,XCnt]] # only use values where ROI==1
                meanv = np.mean(v_values.ravel());

                pred_u = relshift_firstimg_lowerleft_corner_x_diff-relshift_firstimg_lowerleft_corner_x_ref
                pred_v = relshift_firstimg_lowerleft_corner_y_diff-relshift_firstimg_lowerleft_corner_y_ref

                print("idx1 = %d; idx2 = %d; XCnt = %d; meanu = %f um; predu = %f um; meanv= %f um; predv=%f um" % (idx1,idx2,XCnt,meanu*1e6,predu*1e6,meanv*1e6,predv*1e6))
                pass
            pass
        pass
    pass
                
def CalcInitialModel(nloads,CTODs,
                     load1,load2,
                     Xposvecs,CrackCenterX,dic_dy,
                     dic_span,Symmetric_COD,side,
                     relshift_firstimg_lowerleft_corner_x_ref=None,
                     nominal_length=2e-3,nominal_modulus=100.0e9,nominal_stress=50e6,doplots=False):
    """ side=1 (smaller X values) or side=2 (larger X values)"""
    
    InitialModels=np.zeros((nloads,nloads),dtype='O')
    InitialModels[:,:]=None
    InitialCoeffs=np.zeros((2,nloads,nloads),dtype='d') # Zeroth axis is c5, yt
    InitialCoeffs[:,:,:]=np.nan
    
    Error=np.zeros((nloads,nloads),dtype='d')
    Error[:,:] = np.nan

    npoints = np.zeros((nloads,nloads),dtype=np.uint32)

    XPositions = np.zeros((nloads,nloads),dtype='O')
    XPositions[...]=None
    CTODValues = np.zeros((nloads,nloads),dtype='O')
    CTODValues[...]=None

    

    for idx1 in range(nloads):
        #if idx1 != 0:
        #    continue
        for idx2 in range(idx1+1,nloads):
            
            #if idx2 != nloads-1:
            #    continue
            
            XPositions[idx1,idx2]=np.array([],dtype='f')
            CTODValues[idx1,idx2]=np.array([],dtype='d')
            for XCnt in range(Xposvecs.shape[1]):
                #if XCnt != 1:
                #    continue
                Xposvec=Xposvecs[:,XCnt]
                

                if relshift_firstimg_lowerleft_corner_x_ref is not None:
                    # For this load case we have a particular registration correction to apply
                    Xposvec += relshift_firstimg_lowerleft_corner_x_ref[idx1,idx2]
                    pass

                load_diff = load2[idx1,idx2,XCnt]-load1[idx1,idx2,XCnt]
                
                CTOD=CTODs[:,idx1,idx2,XCnt]-(load_diff/nominal_modulus)*(dic_dy*dic_span)
                #CTOD=CTODs[:,idx1,idx2,XCnt]
                if side==1:
                    valid_locations = (Xposvec < CrackCenterX) & (~np.isnan(CTOD))
                    pass
                else:
                    assert(side==2)
                    valid_locations = (Xposvec > CrackCenterX) & (~np.isnan(CTOD))
                    pass
                
                XPositions[idx1,idx2]=np.concatenate((XPositions[idx1,idx2],Xposvec[valid_locations]))
                CTODValues[idx1,idx2]=np.concatenate((CTODValues[idx1,idx2],CTOD[valid_locations]))
                pass


            if XPositions[idx1,idx2].shape[0] == 0:
                # No Data!
                continue

            x0=(np.sqrt(nominal_length)/nominal_modulus,np.mean(XPositions[idx1,idx2]))
            #x0=(1.0e-13,np.mean(XPositions[idx1,idx2]))
            (c5,xt)=initial_fit.fit_initial_model(x0,XPositions[idx1,idx2],load1[idx1,idx2,XCnt],load2[idx1,idx2,XCnt],CrackCenterX,Symmetric_COD,side,CTODValues[idx1,idx2],nominal_length,nominal_modulus,nominal_stress)
            print("side=%d; xt=%f" % (side,xt))
            InitialModels[idx1,idx2]=initial_fit.initial_model((c5,xt),XPositions[idx1,idx2],load1[idx1,idx2,XCnt],load2[idx1,idx2,XCnt],CrackCenterX,Symmetric_COD,side)
            #InitialModels[idx2,idx1]=-initial_fit.initial_model((c5,xt),XPositions,load1[idx1,idx2,XCnt],load2[idx1,idx2,XCnt],side)
            InitialCoeffs[:,idx1,idx2]=(c5,xt)
            npoints[idx1,idx2]=XPositions[idx1,idx2].shape[0]
            Error[idx1,idx2]=np.sum((CTODValues[idx1,idx2]-InitialModels[idx1,idx2])**2.0)/npoints[idx1,idx2]
            #import pdb
            #pdb.set_trace()
            #junk = initial_fit.initial_model((c5,xt),XPositions[idx1,idx2],load1[idx1,idx2,XCnt],load2[idx1,idx2,XCnt],CrackCenterX,Symmetric_COD,side)
            if doplots:
                # Do random sub-percentage
                if np.random.rand() < .05:
                    
                    from matplotlib import pyplot as pl
                    
                    pl.figure()
                    XPositionsSort=np.argsort(XPositions[idx1,idx2])
                    XPositionsSorted=XPositions[idx1,idx2][XPositionsSort]
                    #CTODValuesSorted=CTODValues[idx1,idx2][XPositionsSort]
                    InitialModelValuesSorted=InitialModels[idx1,idx2][XPositionsSort]
                    pl.plot(XPositions[idx1,idx2]*1e3,CTODValues[idx1,idx2]*1e6,'.',
                            XPositionsSorted*1e3,InitialModelValuesSorted*1e6,'-')
                    pl.xlabel('X (mm)')
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
            XPositions,
            CTODValues)


def EvalEffectiveTip(minload,maxload,full_model_params,sigma):
    # create (t,c,k) for scipy splev
    t=np.array([minload]*4 + [maxload]*4,dtype='d')  # four copies of minload followed by four copies of maxload

    splinecoeff=full_model_params[:4]
    #c5=params[4]

    c=np.concatenate((splinecoeff,[0.0]*4))
    k=3
    tck = (t,c,k)
    
    xt = scipy.interpolate.splev(sigma,tck)
    return xt

def InitializeFullModel(load1,load2,TipCoords1,TipCoords2,InitialCoeffs,Error,npoints,XPositions,CTODValues,InitialModels,CrackCenterX,tip_tolerance,Symmetric_COD,side,doplots=True):
    # Perform fit to the results of the Initial models,
    # to seed the full model:
    #
    # NOTE: Actual execution of the full model is no longer done.
    # We just use the result from this initialization

    if Symmetric_COD:
        min_c5 = 0.1/1000e9
        pass
    else:
        min_c5 = np.sqrt(2*50e-6)/1000e9
        pass
    
    (idx1grid,idx2grid)=np.meshgrid(np.arange(InitialCoeffs.shape[1]),np.arange(InitialCoeffs.shape[2]),indexing="ij")

    XPositions_min = np.array( [ np.min(XPosArray) if XPosArray is not None and len(XPosArray) > 0 else np.inf for XPosArray in XPositions.ravel()]).reshape(XPositions.shape)
    XPositions_max = np.array( [ np.max(XPosArray) if XPosArray is not None and len(XPosArray) > 0 else -np.inf for XPosArray in XPositions.ravel()]).reshape(XPositions.shape)


    Signal = np.array( [ np.mean(InitialModelArray) if InitialModelArray is not None else 0.0 for InitialModelArray in InitialModels.ravel() ]).reshape(XPositions.shape)

    


    Noise = np.zeros(XPositions.shape,dtype='d')
    SNR = -np.ones(XPositions.shape,dtype='d') * np.inf
    for id1 in range(XPositions.shape[0]):
        for id2 in range(XPositions.shape[1]):
            if CTODValues[id1,id2] is None or InitialModels[id1,id2] is None:
                continue
            Noise[id1,id2] = np.sqrt(np.mean((InitialModels[id1,id2]-CTODValues[id1,id2])**2.0)) 
            SNR[id1,id2] = Signal[id1,id2]/Noise[id1,id2]
            pass
        pass
    
    
    # Unwrap the error and yt coefficients
    #valid = (~np.isnan(InitialCoeffs[1,:,:].ravel())) & (npoints.ravel() > 20) &  (InitialCoeffs[0,:,:].ravel() >= min_c5)

    c5_or_neginf = InitialCoeffs[0,:,:].ravel()
    c5_or_neginf[np.isnan(c5_or_neginf)] = -np.inf
    valid = (~np.isnan(InitialCoeffs[1,:,:].ravel())) & (npoints.ravel() > 20) &  (c5_or_neginf >= min_c5)
    
    # Use only data points for which xt is inside data range for initial fit and with good SNR

    xt_or_inf = InitialCoeffs[1,:,:].ravel()
    xt_or_inf[np.isnan(xt_or_inf)]=np.inf  # convert NaN to inf to avoid warning

    xt_or_neginf = InitialCoeffs[1,:,:].ravel()
    xt_or_neginf[np.isnan(xt_or_inf)]=-np.inf  # convert NaN to inf to avoid warning
    if side < 1.5: # left side
        for_initial_fit = valid & ( xt_or_neginf >= TipCoords1[0]-tip_tolerance) &  ( xt_or_inf <= CrackCenterX+tip_tolerance) & (SNR.ravel() > 1.0)
        pass
    else: # right side
        for_initial_fit = valid & ( xt_or_neginf >= CrackCenterX-tip_tolerance) &  ( xt_or_inf <= TipCoords2[0]+tip_tolerance) & (SNR.ravel() > 1.0)
        pass
    Error_unwrapped=Error.ravel()[valid]
    c5_unwrapped = InitialCoeffs[0,:,:].ravel()[valid]
    xt_unwrapped = InitialCoeffs[1,:,:].ravel()[valid]


    idx1_unwrapped = idx1grid[:,:].ravel()[valid]
    idx2_unwrapped = idx2grid[:,:].ravel()[valid]
    
    avg_load=(load1+load2).mean(axis=2)/2.0
    avg_load_unwrapped=avg_load.ravel()[valid]

    ## Now sort them so lowest error comes first
    #errsort=np.argsort(Error_unwrapped)    
    
    #raise ValueError("Break")
    ## Do fit to first 50
    #xt_vals = xt_unwrapped[errsort[:50]]
    #avg_load_vals=avg_load_unwrapped[errsort[:50]]
    #c5_vals = c5_unwrapped[errsort[:50]]
    
    xt_vals = InitialCoeffs[1,:,:].ravel()[for_initial_fit]
    avg_load_vals = avg_load.ravel()[for_initial_fit]
    c5_vals = InitialCoeffs[0,:,:].ravel()[for_initial_fit]

    


    avg_load_vals_sort=np.argsort(avg_load_vals)
    avg_load_vals_sorted=avg_load_vals[avg_load_vals_sort]
    xt_vals_sorted = xt_vals[avg_load_vals_sort]
    
    minload=np.min(load1[~np.isnan(load1)].ravel())
    maxload=np.max(load1[~np.isnan(load1)].ravel())
    
    # Use scipy.interpolate.splrep to do fit.
    # task=-1 means we supply interior knots (there aren't any)
    # k=3 means 3rd order
    (t,c,k) = scipy.interpolate.splrep(avg_load_vals_sorted,xt_vals_sorted,xb=minload,xe=maxload,k=3,task=-1,t=np.array([],dtype='f'))
    assert((t[:4]==minload).all())
    assert((t[4:]==maxload).all())

    assert((c[4:]==0).all()) # This spline only has 4 coefficients. For some reason splrep returns four more that are all zero
    assert(k==3)
    seed_param=(c[0],c[1],c[2],c[3],np.median(c5_vals))
    
    # Plot diagnostics
    if doplots:

        sigmarange=np.linspace(minload,maxload,150)
        fittedvals=EvalEffectiveTip(minload,maxload,seed_param,sigmarange)

        from matplotlib import pyplot as pl

        pl.figure()
        pl.plot(xt_unwrapped*1e3,avg_load_unwrapped/1e6,'x',
                xt_vals*1e3,avg_load_vals/1e6,'o',
                fittedvals*1e3,sigmarange/1e6,'-')
        pl.ylabel('Load (MPa)')
        pl.xlabel('Tip position (mm)')
        pl.legend(('All DIC fit data','yt within data range and good SNR','fit to yt within data range and good SNR'))
        pl.title('xt')
        pl.grid()

        
        fig=pl.figure()
        pl.plot(xt_unwrapped*1e3,avg_load_unwrapped/1e6,'x',picker=5)
                #xt_vals*1e3,avg_load_vals/1e6,'o',
                #fittedvals*1e3,sigmarange/1e6,'-')
        pl.ylabel('Load (MPa)')
        pl.xlabel('Tip position (mm)')
        pl.grid()
        pl.title('xt (pickable)')
        def dicfitpick(event):
            thisline=event.artist
            xdata=thisline.get_xdata()
            ydata=thisline.get_ydata()
            indices = event.ind
            print("got indices: %s; side=%d; XPositions[0,1][0]=%f" % (str(indices),side,XPositions[0,1][0]))

            for index in indices:
                pl.figure()

                idx1=idx1_unwrapped[index]
                idx2=idx2_unwrapped[index]
                
                XPositionsSort=np.argsort(XPositions[idx1,idx2])
                XPositionsSorted=XPositions[idx1,idx2][XPositionsSort]
                #CTODValuesSorted=CTODValues[idx1,idx2][XPositionsSort]
                InitialModelValuesSorted=InitialModels[idx1,idx2][XPositionsSort]
                pl.plot(XPositions[idx1,idx2]*1e3,CTODValues[idx1,idx2]*1e6,'.',
                        XPositionsSorted*1e3,InitialModelValuesSorted*1e6,'-')
                pl.xlabel('X (mm)')
                pl.ylabel('CTOD and initial model (um)')
                pl.title('First end of crack: Load1 = %f MPa; load2 = %f MPa' % (load1[idx1,idx2,0]/1.e6,load2[idx1,idx2,0]/1.e6))
                pl.grid()
                
                pass
            pass
        fig.canvas.mpl_connect('pick_event',dicfitpick)
        
        pl.figure()
        pl.plot(avg_load_unwrapped/1e6,c5_unwrapped,'x',
                avg_load_vals/1e6,c5_vals,'o')
        pl.xlabel('Load (MPa)')
        pl.title('c5')
        pl.grid()
        pass


    return (minload,maxload,seed_param)
    
def CalcFullModel(load1,load2,InitialCoeffs,Error,npoints,XPositions,CTODValues,InitialModels,CrackCenterX,Symmetric_COD,side,minload,maxload,seed_param,nominal_length=2e-3,nominal_modulus=100.0e9,nominal_stress=50e6,doplots=True,opencl_ctx=None,opencl_dev=None):
    # Our model (asymmetric case) is dCOD/dsigma = C5*sqrt(x-xt)u(x-xt) where u(x) is the unit step
    # This integrates to:
    #  COD2-COD1 = integral_sigma1^sigma2 C5*sqrt(x-xt)*u(x-xt) dsigma
    #   where xt is a function of sigma, because the tip location shifts
    #   as the crack opens and closes

    # Our model (symmetric case) is dCOD/dsigma = C5*sqrt(x-xt)u(x-xt)*sqrt(2xc-xt-x)*u(2xc-xt-x) where u(x) is the unit step
    # This integrates to:
    #  COD2-COD1 = integral_sigma1^sigma2 C5*sqrt(x-xt)*u(x-xt)*sqrt(2xc-xt-x)*u(2xc-xt-x) dsigma
    #   where xt is a function of sigma, because the tip location shifts
    #   as the crack opens and closes

    
    # From Suresh, eq. 9.45 @ theta=pi, uy = (K1/(2.0*E))*sqrt((r/(2*pi)))*(1+nu)*(2kappa+1 +1)
    # where kappa = (3-nu)/(1+nu) for plane stress or
    # kappa=(3-4nu) for plane strain... where nu is Poisson's ratio


    # In the asymmetric case, 
    # COD=2uy = ((KI/E)/sqrt(2*pi))*(1+nu)*(2kappa+2)  *sqrt(x-xt) 
   # where KI = sigma*sqrt(pi*a)
    # Since this is proportional to sigma, for a fixed length crack
    # Our C5 is expected to equal sqrt(pi*a)*((1/E)/sqrt(2*pi))*(1+nu)*(2kappa+2)
    # or C5=(sqrt(2a)/(E))*(1+nu)*(kappa+1)
    #
    # In our "initial" model we use data from a single pair of loads,
    # assuming the crack tip position xt is fixed. The corresponding
    # delta sigma is the difference of the loads. The corresponding
    # sigma is the average of the two loads. 
    #
    # In our full model we use all data from all load pairs instead.
    # The initial model seeds a spline fit for xt(sigma)

    # Use a first cut C5 estimate to filter out any coefficients that
    # optimized to zero for whatever reason

    c=np.zeros(4,dtype='d')
    (c[0],c[1],c[2],c[3],c5_median) = seed_param
    

    full_model_residual_plot=None

    if doplots:
        from matplotlib import pyplot as pl
        full_model_residual_plot=pl.figure()
        pass
    
    # Perform model fit

    full_model_residual_unaccel_normalized=full_model.full_model_residual_normalized
    args_unaccel=(InitialCoeffs,XPositions,CTODValues,np.mean(load1,axis=2),np.mean(load2,axis=2),minload,maxload,CrackCenterX,Symmetric_COD,side,nominal_length,nominal_modulus,nominal_stress,full_model_residual_plot)
    
    if opencl_ctx is None:
        full_model_residual_normalized = full_model_residual_unaccel_normalized
        args = args_unaccel
        pass
    else:
        full_model_residual_normalized=full_model_accel.full_model_residual_accel_normalized
        args=(InitialCoeffs,XPositions,CTODValues,np.mean(load1,axis=2),np.mean(load2,axis=2),minload,maxload,CrackCenterX,Symmetric_COD,side,nominal_length,nominal_modulus,nominal_stress,full_model_residual_plot,opencl_ctx,opencl_dev)

        test_accel=True
        
        if test_accel:

            # Compare residual output based on seed_param, make sure
            # it matches between unaccelerated and accelerated models
            
            error_unaccel = full_model_residual_unaccel_normalized(seed_param,*args_unaccel)
            error_accel = full_model_residual_normalized(seed_param,*args)

            assert( np.abs(error_accel-error_unaccel)/error_unaccel < .001)
            pass
        
        pass

    if Symmetric_COD:
        c5_normalized = c5_median*nominal_modulus
        pass
    else:
        c5_normalized = c5_median*nominal_modulus/np.sqrt(nominal_length)
        pass
    
    seed_param_normalized = (c[0]/nominal_length,
                             c[1]/nominal_length,
                             c[2]/nominal_length,
                             c[3]/nominal_length,
                             c5_normalized)

    
    #full_model_result = scipy.optimize.minimize(full_model_residual,seed_param,args=args,method="nelder-mead",tol=1e-17)
    #full_model_params=full_model_result.x

    # full_model_result = scipy.optimize.minimize(full_model_residual_normalized,seed_param_normalized,args=args,method="nelder-mead",tol=1e-17)
    full_model_result = scipy.optimize.minimize(full_model_residual_normalized,seed_param_normalized,args=args,method="SLSQP",options={"eps": 1e-6,"ftol": 1e-7},
                                                #bounds=((0.0,None),((np.min(Yposvec)-nominal_length)/nominal_length,(np.max(Yposvec)+nominal_length)/nominal_length))
                                                )
    
    full_model_params_normalized=full_model_result.x

    if Symmetric_COD:
        c5_unnormalized = full_model_params_normalized[4]/nominal_modulus
        pass
    else:
        c5_unnormalized = full_model_params_normalized[4]*np.sqrt(nominal_length)/nominal_modulus
        pass
    
    full_model_params = (full_model_params_normalized[0]*nominal_length,
                         full_model_params_normalized[1]*nominal_length,
                         full_model_params_normalized[2]*nominal_length,
                         full_model_params_normalized[3]*nominal_length,
                         c5_unnormalized)
    #full_model_result=None
    #full_model_params = seed_param

    #sys.modules["__main__"].__dict__.update(globals())
    #sys.modules["__main__"].__dict__.update(locals())
    #raise ValueError("Break!")


        # Plot diagnostics
    if doplots:
        
        sigmarange=np.linspace(minload,maxload,150)
        fittedvals=EvalEffectiveTip(minload,maxload,seed_param,sigmarange)

        from matplotlib import pyplot as pl

        pl.figure()
        pl.plot(xt_unwrapped*1e3,avg_load_unwrapped/1e6,'x',
                xt_vals*1e3,avg_load_vals/1e6,'o',
                fittedvals*1e3,sigmarange/1e6,'-',
                full_model.full_model_xt(full_model_params,sigmarange,minload,maxload)*1e3,sigmarange/1.e6,'-')
        pl.ylabel('Load (MPa)')
        pl.xlabel('Tip position (mm)')
        pl.legend(('All DIC fit data','yt within data range and good SNR','initial fit to yt within data range and good SNR','full model'))
        pl.title('xt')
        pl.grid()


    
    
    return (full_model_params,full_model_result)
                     
