import numpy as np
import scipy
import scipy.optimize

# Calculates the crack depth oping based on the elliptical crack opening model proposed and used by Ravichandran (1997) which is able to handling cracks with an aspect ratio beyond both less than and greater than one.  An if statement seperates the two models and the model is capable of conducting optimizations between the two solutions.  Since Ravichandran's model was origional developed for a fully opened fatigue crack, but for this work the closure points for both the surface and depth opening are subsituted for the total crack dimensions in both directions.  This was verified using FEA solutions as a valid technique and produces accurate calcualtions of the crack depth opening.  The model works as a least squared regression where for a given surface closure point and known surface compliance, the unique crack depth opening can be determined.  

def ellipitcal_model(param,ClosurePoint,CrackCenterX,HalfWidth,SampleThickness,YoungsModulus):
    #Elliptical crack opening model proposed by Ravichandran.  Configured to handle cases for a crack of any aspect ratio.  However, there are potential issues with any solution beyond a depth/surface aspect ratio of 3 due to how the optimization works.  Any solution beyond this value should be considered numerical error.  Due to the model development, it is necessary to center the data bout the crack center.  Thie model is currently set up to be able to handle both left and right crack sides.
    DepthOpening=param
    
    aoc = (DepthOpening)/(np.abs(ClosurePoint-CrackCenterX))
    aot = (DepthOpening)/(SampleThickness)
    coa = 1/aoc

    #g = 1+(0.1+0.35*aot**2)
    #M1 = 1.13-0.09*aoc
    #M2 = -0.54+(0.89/(0.2+aoc))
    #M3 = 0.5-(1/(0.65+aoc))+14*(1-aoc)**24
    #Q = 1+1.464*aoc**1.65
    if aoc <= 1:
        lamb = aoc
        g = 1+(0.1+0.35*aot**2)
        M1 = 1.13-0.09*aoc
        M2 = -0.54+(0.89/(0.2+aoc))
        M3 = 0.5-(1/(0.65+aoc))+14*(1-aoc)**24
        Q = 1+1.464*aoc**1.65
    else:
        lamb = np.sqrt(aoc)
        g = 1+(0.1+0.35*coa*aot**2)
        M1 = np.sqrt(coa)*(1+0.04*coa)
        M2 = 0.2*coa**4
        M3 = -0.11*coa**4
        Q = 1+1.464*coa**1.65
        pass
        
    F = (M1+M2*aot**2+M3*aot**4)*np.sqrt(np.pi/Q)
    fw = np.sqrt(np.cos((np.pi*(np.abs(ClosurePoint-CrackCenterX)))/(2*HalfWidth)*np.sqrt(aot))**-1.0)

    xi = 1.6*np.sqrt(8)/np.sqrt(np.pi)*(1-0.3**2)*lamb*F*fw*g

    #c5val = xi/YoungsModulus
    return xi

def elliptical_residual(param,ClosurePoint,CrackCenterX,HalfWidth,SampleThickness,YoungsModulus,c5measured):
    #Calculation of the crack depth openign using a least squared regression.  In order to obtain a solution the value of c5 meaured, the leading coefficient of the model, must be multiplied by the youngs modulus of the material.  In the COD models the C5measured value has units of 1/(Pa) but this results in a term on the order of ~1e-11.  Any optimization attempted within this range of values will return the input value as the least squared residual is on the order of ~1e-22.  To avoid this, Ravichandran proposed a dimensionless value for the surface compliance and the result from the COD model can be converted into this form by multiplying the leading coefficent by the youngs modulus of the material.
    modelval = ellipitcal_model(param,ClosurePoint,CrackCenterX,HalfWidth,SampleThickness,YoungsModulus)
    #from matplotlib import pyplot as pl
    #pl.figure()
    #pl.plot(modelvals,'-',
    #        CTOD,'-')
    residual = (modelval-c5measured*YoungsModulus)**2
    ret = np.mean(residual)
    #print(c5measured*YoungsModulus,modelval,residual,ret)
    #import pdb
    #pdb.set_trace()
    #ret = np.sum(residual[~np.isnan(residual)])
    #ret = np.mean(residual)
    #pl.title('residual=%g' % (ret))
    return ret

def fit_elliptical_model(seed_param,ClosurePoint,CrackCenterX,HalfWidth,SampleThickness,YoungsModulus,c5measured):

    #Calculation of the depth closure point for a 3D elliptlical surface crack
    
    res=scipy.optimize.minimize(elliptical_residual,seed_param,args=(ClosurePoint,CrackCenterX,HalfWidth,SampleThickness,YoungsModulus,c5measured),method="SLSQP",options={"eps": 1e-6,"ftol": 1e-7})
    #print(c5measured)
    #(c5,xt)=res.x

    #import pdb
    #pdb.set_trace()
    return (res.x)
