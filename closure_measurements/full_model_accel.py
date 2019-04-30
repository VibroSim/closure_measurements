import sys
import os
import os.path
import numpy as np


import scipy
import scipy.integrate

# regenerate qagse_fparams.c with:
# f2c -a qagse_fparams.f
# patch -p0 <qagse_fparams.patch
# NOTE: kernelpattern is processed with
# '%' (printf) operator, so
# be careful about percent signs in it...
# '%s's replaced by contents of qagse_fparams.c
kernelpattern=r"""

#ifdef static 
#undef static
#endif

#define static const __constant // f2c generates static when it means const


typedef double doublereal;  
typedef float real; // perform (almost) all computations in single precision
typedef int integer;
typedef int logical;

#ifndef NULL
#define NULL ((char *)0)
#endif

int assert(int a) {
  char *null=NULL;
  if (!a) { 
    if (*null) return 0;// attempt to read from invalid address zero
  }
  return 1;
 }

//typedef real (*E_fp)();
typedef char *E_fp;  // Don't use this anymore... hardwired to funct_(...)

float dabs(float p) { return fabs(p); }
float dmax(float p,float q) { if (p > q) return p;else return q; }
float dmin(float p,float q) { if (p < q) return p;else return q; }

double pow_dd(doublereal *arg1, const __constant doublereal *arg2)
{
  return pow(*arg1,*arg2);
}

/* C source for R1MACH -- remove the * in column 1 */
doublereal r1mach_(const __constant integer *i)
{
	switch(*i){
	  case 1: return FLT_MIN;
	  case 2: return FLT_MAX;
	  case 3: return FLT_EPSILON/FLT_RADIX;
	  case 4: return FLT_EPSILON;
	  case 5: return log10((float)FLT_RADIX);
	  }
        printf("invalid argument: r1mach(%%ld)\n", *i);
        return 0.0f;
//	assert(0); return 0; /* else complaint of missing return value */
}


#define TRUE_ 1
#define FALSE_ 0

#define integrand_sigmaintegral full_model_kernel


//#define LIMIT 7000 // max # of integration intervals

#define qagse_ qagse_sigmaintegral
#define qelg_ qelg_sigmaintegral
#define qk21_ qk21_sigmaintegral
#define qpsrt_ qpsrt_sigmaintegral
#define funct_ integrand_sigmaintegral
#define c__4 c__4_sigmaintegral
#define c__1 c__1_sigmaintegral
#define c__2 c__2_sigmaintegral
#define c_b42 c_b42_sigmaintegral
#define xgk xgk_sigmaintegral
#define wgk wgk_sigmaintegral
#define wg wg_sigmaintegral

%s 

#undef qagse_ 
#undef qelg_ 
#undef qk21_ 
#undef qpsrt_ 
#undef funct_ 
#undef c__4
#undef c__1
#undef c__2
#undef c_b42
#undef xgk
#undef wgk
#undef wg

// See email to Chris Giuffre with MATLAB code for evaluating bsplines, 2019-03-02, 
// and also ~/research/software/bspline_matlab
#if 0 // disable recursive version because OpenCL can't handle recursion (use separate functions for k=4,3, 2, and 1
// NOTE: k here is k+1 as defined for splev
float Bfunc_recursive(float *knots,unsigned i,unsigned k,float xpos)
{
  float Bkm1; // B_{i,k-1} 
  float Bip1km1; // B_{i+1,k-1}
  float Bk;

  if (k==1) {
    if (xpos >= knots[i-1] && xpos <= knots[i]) {  // ***!!! Note <= instead of < because we only allow one segment and want to define all the way out to the end!!!
      // Within segment
      return 1.0;
    } else {
      // outside segment
      return 0.0; 
    }
  }
  
  Bkm1=Bfunc(knots,i,k-1,xpos);

  Bip1km1=Bfunc(knots,i+1,k-1,xpos);

  Bk=0.0;

  if (knots[i+k-1-1] != knots[i-1]) {
    Bk += ((xpos-knots[i-1])/(knots[i+k-1-1] - knots[i-1]))*Bkm1;
  }

  if (knots[i+k-1] != knots[i]) {
    Bk += ((knots[i+k-1]-xpos)/(knots[i+k-1]-knots[i]))*Bip1km1;
  }
  return Bk;

}
#endif // if 0

float Bfunc_k1(float *knots,unsigned i,float xpos)
// k==1 case
{
  float Bkm1; // B_{i,k-1} 
  float Bip1km1; // B_{i+1,k-1}
  float Bk;

  if (xpos >= knots[i-1] && xpos <= knots[i]) {  // ***!!! Note <= instead of < because we only allow one segment and want to define all the way out to the end!!!
    // Within segment
    return 1.0;
  } else {
    // outside segment
    return 0.0; 
  }
  
}


float Bfunc_k2(float *knots,unsigned i,float xpos)
// k==2 case
{
  float Bkm1; // B_{i,k-1} 
  float Bip1km1; // B_{i+1,k-1}
  float Bk;
  unsigned k=2; 

  Bkm1=Bfunc_k1(knots,i,xpos);

  Bip1km1=Bfunc_k1(knots,i+1,xpos);

  Bk=0.0;

  if (knots[i+k-1-1] != knots[i-1]) {
    Bk += ((xpos-knots[i-1])/(knots[i+k-1-1] - knots[i-1]))*Bkm1;
  }

  if (knots[i+k-1] != knots[i]) {
    Bk += ((knots[i+k-1]-xpos)/(knots[i+k-1]-knots[i]))*Bip1km1;
  }
  return Bk;

}

float Bfunc_k3(float *knots,unsigned i,float xpos)
// k==3 case
{
  float Bkm1; // B_{i,k-1} 
  float Bip1km1; // B_{i+1,k-1}
  float Bk;
  unsigned k=3; 

  Bkm1=Bfunc_k2(knots,i,xpos);

  Bip1km1=Bfunc_k2(knots,i+1,xpos);

  Bk=0.0;

  if (knots[i+k-1-1] != knots[i-1]) {
    Bk += ((xpos-knots[i-1])/(knots[i+k-1-1] - knots[i-1]))*Bkm1;
  }

  if (knots[i+k-1] != knots[i]) {
    Bk += ((knots[i+k-1]-xpos)/(knots[i+k-1]-knots[i]))*Bip1km1;
  }
  return Bk;

}

float Bfunc_k4(float *knots,unsigned i,float xpos)
// k==4 case
{
  float Bkm1; // B_{i,k-1} 
  float Bip1km1; // B_{i+1,k-1}
  float Bk;
  unsigned k=4; 

  Bkm1=Bfunc_k3(knots,i,xpos);

  Bip1km1=Bfunc_k3(knots,i+1,xpos);

  Bk=0.0;

  if (knots[i+k-1-1] != knots[i-1]) {
    Bk += ((xpos-knots[i-1])/(knots[i+k-1-1] - knots[i-1]))*Bkm1;
  }

  if (knots[i+k-1] != knots[i]) {
    Bk += ((knots[i+k-1]-xpos)/(knots[i+k-1]-knots[i]))*Bip1km1;
  }
  return Bk;

}



doublereal full_model_kernel(float *sigma, 
                             float *YPosition, 
                             float *minload,float *maxload, 
                             float *c1,float *c2, float *c3, float *c4, 
                             float *c5,float *side)
{
  float t[8];
  float B14,B24,B34,B44;
  float yt; 
  float sqrtarg;
  
  t[0]=t[1]=t[2]=t[3]=*minload;
  t[4]=t[5]=t[6]=t[7]=*maxload;

  B14 = Bfunc_k4(t,1,*sigma);
  B24 = Bfunc_k4(t,2,*sigma);
  B34 = Bfunc_k4(t,3,*sigma);
  B44 = Bfunc_k4(t,4,*sigma);

  yt = (*c1)*B14 + (*c2)*B24 + (*c3)*B34 + (*c4)*B44;

  if (*side < 1.5) { // left side, position > tip position
    sqrtarg = *YPosition - yt;
  } else { // right side, position < tip position
    sqrtarg = yt - *YPosition;
  }

  if (sqrtarg < 0.0) {
    sqrtarg=0.0;
  }

  return (*c5)*sqrt(sqrtarg);
}


__kernel void integrate_full_model_kernel(float load1,
                                          float load2,
                                          __global const float *YPositions,
                                          float minload,float maxload,
                                          float c1,float c2, float c3, float c4,
                                          float c5,
                                          float side,  // positive 1: RHS, negative 1: LHS
                                          int limit,
                                          __global float *flists,
                                          __global int *iords,
                                          __global float *out)
{
  float result=0.0;
  float epsrel=1e-6f; 
  float epsabs=1e-8f; // !!!*** Do we need to set this more sensibly?
  //float alist[LIMIT];
  __global float *alist,*blist,*rlist,*elist;
  __global int *iord;
  //float blist[LIMIT];
  //float rlist[LIMIT];
  //float elist[LIMIT];
  //int iord[LIMIT];
  int last=0;
  //int limit=LIMIT;
  float abserr=0.0;
  int neval=0;
  int ier=0;
  float YPosition;
  size_t Y_idx;

  Y_idx = get_global_id(0);
  YPosition=YPositions[Y_idx];

  alist = flists + Y_idx*4*limit;
  blist = flists + Y_idx*4*limit + limit;
  rlist = flists + Y_idx*4*limit + 2*limit;
  elist = flists + Y_idx*4*limit + 3*limit;
  iord = iords + Y_idx*limit;

  qagse_sigmaintegral(NULL,
                               &YPosition,&minload,&maxload, // this start/end are the bounds of the splines
                               &c1,&c2,&c3,&c4,&c5,
                               &side,
                               &load1,&load2, // this start/end are the bounds of the integration
                               &epsabs,&epsrel,
                               &limit,
                               &result,
                               &abserr,
                               &neval,
                               &ier,
                               alist,blist,rlist,elist,iord,
                               &last);
  if (ier != 0) {
    printf("sigmaintegral: ier=%%d; Y_idx=%%u; YPosition=%%f\n",ier,Y_idx,YPosition);
  }
  out[Y_idx]=result;   
}
"""

# find current module so we can use path to load "qagse_fparams.c"
class dummy(object):
    pass
modpath = sys.modules[dummy.__module__].__file__
moddir = os.path.split(modpath)[0]

kernelprog=None # OpenCL program; this global is used to cache the compiled copy

def get_kernelprog(ctx):
    global kernelprog

    if kernelprog is None:
        import pyopencl as cl
        
        qagse_fparams=open(os.path.join(moddir,"qagse_fparams.c"),"r").read()
        kernelcode=kernelpattern % (qagse_fparams)
        prg=cl.Program(ctx,kernelcode)
        prg.build()
        
        kernelprog=prg
        pass
    return kernelprog

    


def full_model_kernel(sigma,YPosition,c5,tck,side):
    """This is the integrand of: 
          integral_sigma1^sigma2 C5*sqrt(y-yt)*u(y-yt) dsigma
        where yt is a function of sigma given by the spline
        coefficents tck
        """
    yt = scipy.interpolate.splev(sigma,tck)
    if side < 1.5:
        sqrtarg = YPosition-yt
        pass
    else:
        sqrtarg = yt-YPosition
        pass
    
    #sqrtarg[sqrtarg < 0.0] = 0.0
    if sqrtarg < 0.0:
        sqrtarg=0.0
        pass
    
    modelvals = c5*np.sqrt(sqrtarg)
    return modelvals


def plot_full_model_residual(params,YPositions,CTODValues,load1,load2,minload,maxload,side,full_model_residual_plot,err):
    from matplotlib import pyplot as pl

    splinecoeff=params[:4]
    c5=params[4]
    
    # create (t,c,k) for scipy splev
    t=np.array([minload]*4 + [maxload]*4,dtype='d')  # four copies of minload followed by four copies of maxload
    c=np.concatenate((splinecoeff,[0.0]*4))
    k=3
    tck = (t,c,k)


    avg_load=(load1+load2)/2.0
    
    
    fig=pl.figure(full_model_residual_plot.number)
    pl.clf()
    #loads=np.linspace(minload,maxload,20)
    #pl.plot(loads/1e6,scipy.interpolate.splev(loads,tck)*1e3,'-',
    pl.plot(avg_load.ravel()/1e6,scipy.interpolate.splev(avg_load.ravel(),tck)*1e3,'x',picker=5)
    
    pl.grid()
    pl.title('params=%s\nerr=%g (clickable)' % (str(params),err))
    pl.xlabel('load (MPa)')
    pl.ylabel('Tip position (mm)')
    def dicfitpick(event):
        thisline=event.artist
        xdata=thisline.get_xdata()
        ydata=thisline.get_ydata()
        indices = event.ind
        print("got indices: %s; side=%d" % (str(indices),side))
        
        for index in indices:
            
            idx1=index // avg_load.shape[1]
            idx2=index % avg_load.shape[1]

            if YPositions[idx1,idx2] is None:
                # No data here
                continue

            pl.figure()
            
            #load=avg_load[idx1,idx2]

            #sys.modules["__main__"].__dict__.update(globals())
            #sys.modules["__main__"].__dict__.update(locals())
            YPositionsSort=np.argsort(YPositions[idx1,idx2])
            YPositionsSorted=YPositions[idx1,idx2][YPositionsSort]
            #CTODValuesSorted=CTODValues[idx1,idx2][YPositionsSort]
            #InitialModelValuesSorted=InitialModels[idx1,idx2][YPositionsSort]

            integralvals = np.array([ scipy.integrate.quad(full_model_kernel,load1[idx1,idx2],load2[idx1,idx2],(YPosition,c5,tck,side))[0]  for YPosition in YPositionsSorted ],dtype='d')

            pl.plot(YPositions[idx1,idx2]*1e3,CTODValues[idx1,idx2]*1e6,'.',
                    YPositionsSorted*1e3,integralvals*1e6,'-')
            pl.xlabel('Y (mm)')
            pl.ylabel('CTOD and full model (um)')
            pl.title('First end of crack: Load1 = %f MPa; load2 = %f MPa' % (load1[idx1,idx2]/1.e6,load2[idx1,idx2]/1.e6))
            pl.grid()
            
            
            pass
        pass

    # disconnect prior callback handlers
    if "pick_event" in fig.canvas.callbacks.callbacks:
        for cid in fig.canvas.callbacks.callbacks["pick_event"].keys():
            fig.canvas.mpl_disconnect(cid)
            pass
        pass
    
    
    fig.canvas.mpl_connect('pick_event',dicfitpick)

    
    #full_model_residual_plot.canvas.draw()
    #full_model_residual_plot.canvas.flush_events()
    pl.savefig('/tmp/loadplot.png',dpi=300)
    pass




def full_model_residual_accel_normalized(params,YPositions,CTODValues,load1,load2,minload,maxload,side,nominal_length,nominal_modulus,nominal_stress,full_model_residual_plot,opencl_ctx,opencl_dev):
    splinecoeff_normalized=params[:4]
    c5_normalized=params[4]

    params_unnormalized=(splinecoeff_normalized[0]*nominal_length,
                         splinecoeff_normalized[1]*nominal_length,
                         splinecoeff_normalized[2]*nominal_length,
                         splinecoeff_normalized[3]*nominal_length,
                         c5_normalized*np.sqrt(nominal_length)/nominal_modulus)

    #  unnormalized result is average over all load pairs of (integral_sigma1^sigma2 C5*sqrt(y-yt)*u(y-yt) dsigma - CTOD)^2... i.e. mean of squared CTODs
    
    nominal_ctod = nominal_length*nominal_stress/nominal_modulus

    
    return full_model_residual_accel(params_unnormalized,YPositions,CTODValues,load1,load2,minload,maxload,side,full_model_residual_plot,opencl_ctx,opencl_dev)/(nominal_ctod**2.0)
    


def full_model_residual_accel(params,YPositions,CTODValues,load1,load2,minload,maxload,side,full_model_residual_plot,opencl_ctx,opencl_dev):


    import pyopencl as cl
    mf = cl.mem_flags

    err=0.0
    numpos=0
    
    prg=get_kernelprog(opencl_ctx);

    queue=cl.CommandQueue(opencl_ctx,opencl_dev)

    
    c1=params[0]
    c2=params[1]
    c3=params[2]
    c4=params[3]
    c5=params[4]
    

    limit=7000

    for idx1 in range(load1.shape[0]):

        # Lists of buffers that are being used by OpenCL
        YPos_bufs=[]
        CTOD_arrays=[]
        out_arrays=[]
        out_bufs=[]
        out_events=[]
        
        for idx2 in range(idx1+1,load1.shape[1]):
            # At this load pair, have array of data
            # over y: YPositions[idx1,idx2]
            # and CTODValues[idx1,idx2]
            
            # Calculate the model value over the
            # various Y positions:
            #  integral_sigma1^sigma2 C5*sqrt(y-yt)*u(y-yt) dsigma
            assert(len(YPositions[idx1,idx2].shape)==1)
            
            if YPositions[idx1,idx2].shape[0]==0:
                continue  # Empty y list... nothing to do here
                
            #sys.modules["__main__"].__dict__.update(globals())
            #sys.modules["__main__"].__dict__.update(locals())
            #raise ValueError("Break!")
            
            kern = prg.integrate_full_model_kernel
            kern.set_scalar_arg_dtypes([np.float32,np.float32,
                                        None,
                                        np.float32,np.float32,
                                        np.float32,np.float32,np.float32,np.float32,
                                        np.float32,
                                        np.float32,
                                        np.int32,
                                        None,
                                        None,
                                    None])
        
            assert(YPositions[idx1,idx2].flags.contiguous)
            assert(YPositions[idx1,idx2].dtype==np.float32)
        
            YPositions_buf = cl.Buffer(opencl_ctx,mf.READ_ONLY,size=YPositions[idx1,idx2].nbytes)
            YCopyEv=cl.enqueue_copy(queue,YPositions_buf,YPositions[idx1,idx2],is_blocking=False);

            flists_buf = cl.Buffer(opencl_ctx,mf.READ_WRITE,size=4*4*limit*YPositions[idx1,idx2].shape[0]) # 4 LIMIT-length float buffers per work element
            
            iord_buf = cl.Buffer(opencl_ctx,mf.READ_WRITE,size=4*limit*YPositions[idx1,idx2].shape[0]) # 1 LIMIT-length int buffer per work element

        
            out_array=np.empty(YPositions[idx1,idx2].shape,dtype='f')
            out_buf=cl.Buffer(opencl_ctx,mf.WRITE_ONLY,size=out_array.nbytes)
        
            #import pdb
            #pdb.set_trace()
            
            KernEv=kern(queue,YPositions[idx1,idx2].shape,None,
                        load1[idx1,idx2],
                        load2[idx1,idx2],
                        YPositions_buf,
                        minload,maxload,
                        c1,c2,c3,c4,
                        c5,side,
                        limit,
                        flists_buf,
                        iord_buf,
                        out_buf,wait_for=(YCopyEv,))
        
            out_event=cl.enqueue_copy(queue,out_array,out_buf,wait_for=(KernEv,),is_blocking=False)
            # err += (integral-CTODValues[idx1,idx2][YPosIdx])**2.0
            queue.flush()
            
            YPos_bufs.append(YPositions_buf)
            CTOD_arrays.append(CTODValues[idx1,idx2])
            out_arrays.append(out_array)
            out_bufs.append(out_buf)
            out_events.append(out_event)
            pass
            
        # Go through list of buffers, waiting for completion
        for pos in range(len(YPos_bufs)):
            YPos_buf=YPos_bufs[pos]
            CTOD_array=CTOD_arrays[pos]
            out_array=out_arrays[pos]
            out_buf=out_bufs[pos]
            out_event=out_events[pos]
            
            out_event.wait()
        
            # Error gets sum squared residual added
            err += np.sum((out_array-CTOD_array)**2.0)
            numpos += out_array.shape[0]
            
            YPos_buf.release()
            out_buf.release()
            #out_event.release()
        
            pass
        pass

    err /= numpos
    
    print("full_model_residual: Calculate at params=%s; err=%g" % (str(params),err))

    if full_model_residual_plot is not None:
        plot_full_model_residual(params,YPositions,CTODValues,load1,load2,minload,maxload,side,full_model_residual_plot,err)
        pass
        
    
    
    return err
