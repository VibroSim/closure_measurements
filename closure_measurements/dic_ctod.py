import numpy as np
import scipy
import scipy.ndimage


def dic_ctod(v_disp,span,window,DIC_ROI_out,ROI_yminidx,ROI_ymaxidx):
    """  Calculates Crack Tip Opening displacement from DIC data
     Parameters:
       v_disp: DIC displacement array... will trim to maximum 
               area
       span: distance in DIC analysis points between analysis windows
       window: square size of window for averaging, in DIC pixels.. Should be odd
       DIC_ROI_out: ROI output from the DIC algorithm, same shape as v_disp
       ROI_yminidx: starting y index of the ROI divided bx the scalefactor
       ROI_ymaxidx: ending y index of the ROI divided bx the scalefactor
       appl_load: load applied to sample 
       #y_pos:  Positions in meters of the various y positions of u_disp 
     Returns: 
       CTOD:  Arrax of the same length as y_pos
       top_disp:   Top displacement prior to subtraction
       bot_disp:   Bottom displacement prior to subtraction
    """

    # Create storage for return values
    top_disp=np.ones(v_disp.shape[0],dtype='d')*np.NaN
    bot_disp=np.ones(v_disp.shape[0],dtype='d')*np.NaN
    CTOD=np.ones(v_disp.shape[0],dtype='d')*np.NaN

    # window our data down according to ROI
    
    # Shift starting point up while it points at all zeros in the ROI out
    while (DIC_ROI_out[:,ROI_yminidx]==0).all():
        ROI_yminidx+=1
        pass

    # Shift ending point down while its predecessor points at all zeros in the ROI out
    while (DIC_ROI_out[:,ROI_ymaxidx-1]==0).all():
        ROI_ymaxidx-=1
        pass

    ROI_xminidx=0;
    # Shift starting point right while there are any zeros in the ROI out
    while ROI_xminidx < DIC_ROI_out.shape[0] and (DIC_ROI_out[ROI_xminidx,ROI_yminidx:ROI_ymaxidx]==0).any():
        ROI_xminidx+=1
        pass


    ROI_xmaxidx=v_disp.shape[0];
    # Shift ending point left while there are any zeros in the ROI out
    while ROI_xmaxidx >= 0 and (DIC_ROI_out[ROI_xmaxidx-1,ROI_yminidx:ROI_ymaxidx]==0).any():
        ROI_xmaxidx-=1
        pass
    
    # We had a failed correlation if we lost more that 20% of the area
    # from this last operation... So only proceed if we have good
    # data (otherwise we just return our default arrays of NaNs)
    if (ROI_xmaxidx-ROI_xminidx) > 0.8 * v_disp.shape[0]:    
    
        # !!! Switch away from filtering because ncorr sometimes leaves little
        # bits of undefined result in the middle of the data.
        # ... We convert these to NaN, then they contaminate their vicinity
        # but it's OK with the simple averaging filter (or any other short FIR
        # filter)
        #u_disp_fft2=np.fft.fft2(u_disp)
        #filtered_fft2=scipy.ndimage.fourier_gaussian(u_disp_fft2,sigma=window/2.0)
        #u_disp_filtered=np.fft.ifft2(filtered_fft2).real
        
        v_disp_filtered=np.empty(v_disp.shape,dtype='d')
        v_disp_filtered[...]=np.nan
        
        window_edge=(window-1)//2
    
        v_disp_nanmask = v_disp.copy()
        v_disp_nanmask[DIC_ROI_out != 1] = np.nan
        
        v_disp_filtered[window_edge:-window_edge,window_edge:-window_edge]=0.0
    
        # Average over the area of the square window
        for windowshiftx in range(-window_edge,window_edge+1):
            for windowshifty in range(-window_edge,window_edge+1):
                v_disp_filtered[(window_edge+windowshiftx):(v_disp_filtered.shape[0]-window_edge+windowshiftx),(window_edge+windowshifty):(v_disp_filtered.shape[1]-window_edge+windowshifty)] += v_disp_nanmask[window_edge:-window_edge,window_edge:-window_edge]
                pass
            pass
        v_disp_filtered /= window**2.0 # divide by the number of items added together
        
        # Cut off outside of ROI
        v_disp_filtered_cropped=v_disp_filtered[:,(ROI_yminidx+window*2):(ROI_ymaxidx-window*2)] 
        
        
        for xnum in range(window+ROI_xminidx,ROI_xmaxidx-window):
            # Iterate over y-positions (i.e. move parallel to crack)
            # but not closer than a window width to the edge
            
            # Find biggest step
            stepsizes = np.abs(np.diff(v_disp_filtered_cropped[xnum,:]))
            stepsizes[np.isnan(stepsizes)]=-np.inf
            maxstep_idx = np.argmax(stepsizes)
            
            # Lower position is "span"/2 piyels down from biggest step
            lower_idx = int(round(maxstep_idx+0.5-span/2.0))
            # Upper position is "span"/2 piyels up from biggest step
            upper_idx = int(round(maxstep_idx+0.5+span/2.0))
            
            if (upper_idx < v_disp_filtered_cropped.shape[1] and
                lower_idx >= 0):
                top_disp[xnum]=v_disp_filtered_cropped[xnum,upper_idx]
                bot_disp[xnum]=v_disp_filtered_cropped[xnum,lower_idx]
                
                
                CTOD[xnum]=top_disp[xnum]-bot_disp[xnum]
                pass
            
            pass
        pass
        


        
    return (CTOD,top_disp,bot_disp)

