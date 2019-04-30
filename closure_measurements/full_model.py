import sys
import os
import os.path
import numpy as np


import scipy
import scipy.integrate

def full_model_yt(params,sigmas,minload,maxload):
    splinecoeff=params[:4]
    
    # create (t,c,k) for scipy splev
    t=np.array([minload]*4 + [maxload]*4,dtype='d')  # four copies of minload followed by four copies of maxload
    c=np.concatenate((splinecoeff,[0.0]*4))
    k=3
    tck = (t,c,k)

    yt = scipy.interpolate.splev(sigmas,tck)
    return yt


def full_model_kernel(sigma,YPosition,c5,tck,CrackCenterY,Symmetric_COD,side):
    """This is the integrand of: 
          integral_sigma1^sigma2 C5*sqrt(y-yt)*u(y-yt) dsigma (asymmetric COD case)
    or 
          integral_sigma1^sigma2 C5*sqrt(y-yt)*u(y-yt)*sqrt(2yc-yt-y)*u(2yc-yt-y) dsigma (symmetric COD case)

        where yt is a function of sigma given by the spline
        coefficents tck and yc is the coordinate of the crack center
        """
    yt = scipy.interpolate.splev(sigma,tck)
    if side < 1.5: # left side, position > tip posiiton
        sqrtarg = YPosition-yt
        pass
    else: # right side, position < tip position
        sqrtarg = yt-YPosition
        pass
    
    #sqrtarg[sqrtarg < 0.0] = 0.0
    if sqrtarg < 0.0:
        sqrtarg=0.0
        pass

    if Symmetric_COD:
        if side==1:
            sqrtarg2 = 2*CrackCenterY-yt-YPosition
            pass
        else:
            sqrtarg2 = YPosition-2*CrackCenterY+yt
            pass
        if sqrtarg2 < 0.0:
            sqrtarg2=0.0
            pass
        
        modelvals = c5*np.sqrt(sqrtarg)*np.sqrt(sqrtarg2)
        # c5 has units of meters of COD per length per Pascal of load

        pass
    else:        
        # c5 has units of meters of COD per sqrt(length) per Pascal of load
        modelvals = c5*np.sqrt(sqrtarg)
        pass
    return modelvals


def full_model_residual_normalized(params,YPositions,CTODValues,load1,load2,minload,maxload,CrackCenterY,Symmetric_COD,side,nominal_length,nominal_modulus,nominal_stress,full_model_residual_plot):
    splinecoeff_normalized=params[:4]
    c5_normalized=params[4]

    if Symmetric_COD:
        # c5 has units of meters of COD per length per Pascal of load
        c5 = c5_normalized/nominal_modulus
        pass
    else:
        # c5 has units of meters of COD per sqrt(length) per Pascal of load
        c5 = c5_normalized*np.sqrt(nominal_length)/nominal_modulus
        pass
    
    params_unnormalized=(splinecoeff_normalized[0]*nominal_length,
                         splinecoeff_normalized[1]*nominal_length,
                         splinecoeff_normalized[2]*nominal_length,
                         splinecoeff_normalized[3]*nominal_length,
                         c5)

    #  unnormalized result is average over all load pairs of (integral_sigma1^sigma2 C5*sqrt(y-yt)*u(y-yt) dsigma - CTOD)^2... i.e. mean of squared CTODs
    
    nominal_ctod = nominal_length*nominal_stress/nominal_modulus

    
    return full_model_residual(params_unnormalized,YPositions,CTODValues,load1,load2,minload,maxload,CrackCenterY,Symmetric_COD,side,full_model_residual_plot)/(nominal_ctod**2.0)
    
    
def full_model_residual(params,YPositions,CTODValues,load1,load2,minload,maxload,CrackCenterY,Symmetric_COD,side,full_model_residual_plot):

    splinecoeff=params[:4]
    c5=params[4]
    
    # create (t,c,k) for scipy splev
    t=np.array([minload]*4 + [maxload]*4,dtype='d')  # four copies of minload followed by four copies of maxload
    c=np.concatenate((splinecoeff,[0.0]*4))
    k=3
    tck = (t,c,k)

    numpos=0
    
    err=0.0    
    for idx1 in range(load1.shape[0]):
        for idx2 in range(idx1+1,load1.shape[1]):
            # At this load pair, have array of data
            # over y: YPositions[idx1,idx2]
            # and CTODValues[idx1,idx2]
            
            # Calculate the model value over the
            # various Y positions:
            #  integral_sigma1^sigma2 C5*sqrt(y-yt)*u(y-yt) dsigma
            
            for YPosIdx in range(len(YPositions[idx1,idx2])):
                # Evaluate integral at this Y position
                integral = scipy.integrate.quad(full_model_kernel,load1[idx1,idx2],load2[idx1,idx2],(YPositions[idx1,idx2][YPosIdx],c5,tck,CrackCenterY,Symmetric_COD,side))[0]
                err += (integral-CTODValues[idx1,idx2][YPosIdx])**2.0
                pass
            numpos+=len(YPositions[idx1,idx2])
            pass
        pass

    err /= numpos
    
    print("full_model_residual: Calculate at params=%s; err=%g" % (str(params),err))

    if full_model_residual_plot is not None:
        from matplotlib import pyplot as pl
        pl.figure(full_model_residual_plot.number)
        pl.clf()
        loads=np.linspace(minload,maxload,20)
        pl.plot(loads/1e6,scipy.interpolate.splev(loads,tck)*1e3,'-')
        pl.grid()
        pl.title('params=%s\nerr=%g' % (str(params),err))
        pl.xlabel('load (MPa)')
        pl.ylabel('Tip position (mm)')

        #full_model_residual_plot.canvas.draw()
        #full_model_residual_plot.canvas.flush_events()
        pl.savefig('/tmp/loadplot.png',dpi=300)
        pass

        

    
    return err

