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
    modelvals = c5*np.sqrt(sqrtarg)*(-(load2-load1))
    return modelvals


def initial_residual(param,Yposvec,load1,load2,side,CTOD):

    modelvals = initial_model(param,Yposvec,load1,load2,side)
    #from matplotlib import pyplot as pl
    #pl.figure()
    #pl.plot(modelvals,'-',
    #        CTOD,'-')
    residual=(modelvals-CTOD)**2
    #import pdb
    #pdb.set_trace()
    ret = np.sum(residual[~np.isnan(residual)])
    #pl.title('residual=%g' % (ret))
    return ret

def fit_initial_model(seed_param,Yposvec,load1,load2,side,CTOD):

    res=scipy.optimize.minimize(initial_residual,seed_param,args=(Yposvec,load1,load2,side,CTOD),method="nelder-mead",tol=1e-17)
    (c5,yt)=res.x
    return (c5,yt)
