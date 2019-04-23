import numpy as np
import scipy
import scipy.optimize

# MODEL: CTOD2-CTOD1 = integral_load1^load2 C5*sqrt(x-xt) dsigma



def initial_model(param,Yposvec,load1,load2,side):
    (c5,yt)=param
    
    if side==1:
        sqrtarg = Yposvec-yt
        pass
    else:
        sqrtarg = yt-Yposvec
        pass
    
    sqrtarg[sqrtarg < 0.0] = 0.0
    modelvals = c5*np.sqrt(sqrtarg)*((load2-load1))
    return modelvals


def initial_residual_normalized(param,Yposvec,load1,load2,side,CTOD,nominal_length,nominal_modulus,nominal_stress):
    """Normalized to be unitless and take unitless parameters
    so as to make minimization numerically more stable and accurate"""

    nominal_ctod = nominal_length*nominal_stress/nominal_modulus
    
    #(c5,yt)=param
    c5=param[0]*np.sqrt(nominal_length)/nominal_modulus
    yt=param[1]*nominal_length

    return initial_residual((c5,yt),Yposvec,load1,load2,side,CTOD)/(nominal_ctod**2.0)


def initial_residual(param,Yposvec,load1,load2,side,CTOD):

    modelvals = initial_model(param,Yposvec,load1,load2,side)
    #from matplotlib import pyplot as pl
    #pl.figure()
    #pl.plot(modelvals,'-',
    #        CTOD,'-')
    residual=(modelvals-CTOD)**2
    #import pdb
    #pdb.set_trace()
    #ret = np.sum(residual[~np.isnan(residual)])
    ret = np.mean(residual)
    #pl.title('residual=%g' % (ret))
    return ret

def fit_initial_model(seed_param,Yposvec,load1,load2,side,CTOD,nominal_length,nominal_modulus,nominal_stress):

    #res=scipy.optimize.minimize(initial_residual,seed_param,args=(Yposvec,load1,load2,side,CTOD),method="nelder-mead",tol=1e-17)
    # The value of c5 is like a sqrt(cracklength)/modulus
    # The value of yt is a like a cracklength
    (c5,yt)=seed_param
    seed_param_normalized=(c5*nominal_modulus/np.sqrt(nominal_length),yt/nominal_length)
    
    
    res=scipy.optimize.minimize(initial_residual_normalized,seed_param_normalized,args=(Yposvec,load1,load2,side,CTOD,nominal_length,nominal_modulus,nominal_stress),method="SLSQP",options={"eps": 1e-6,"ftol": 1e-7},bounds=((0.0,None),((np.min(Yposvec)-nominal_length)/nominal_length,(np.max(Yposvec)+nominal_length)/nominal_length)))
    #(c5,yt)=res.x
    c5=res.x[0]*np.sqrt(nominal_length)/nominal_modulus
    yt=res.x[1]*nominal_length
    #import pdb
    #pdb.set_trace()
    return (c5,yt)
