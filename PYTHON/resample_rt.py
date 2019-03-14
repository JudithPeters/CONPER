import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#%%

subj_id = [2, 3, 4, 5, 6 , 7, 8, 10]
for i_task in range(3):
    plt.figure()
    ac_mean = np.zeros(600)
    for i_subj in range(8):
        file = 'CP_fMRI_task%d_s%.2d.txt' % (i_task + 1, subj_id[i_subj])
        data = pd.read_csv(file, sep="\t", header=None)
        data.columns = ["Pp","Pc","Pt","On","Tk","Tn","C1","C2","D","R1","K1","F1"]
        
        '''
        find the moment when scanning started (participants started behavior
        before scan onset and those trials need to be ignore)
        '''
        idx = np.argmax(np.diff(data.Pt)) + 1
        rt = np.array(data.R1[idx::])
        time_idx = (np.array(data.Pt[idx::]) / 1500).astype(int)     # normalize time to reflect TRs
        time_idx = time_idx - time_idx[0]                            # set scan onset to TR = 0
        
        '''
        arrange reaction times into a temporary array according to their time 
        stamps (measured in TRs).
        All TRs for which observations exist are filled in. 
        The others are filled with NaNs.
        Additionally, both ends of the array are padded with the average 
        reaction time
        '''
        tmp = np.zeros(602) * np.nan            
        tmp[time_idx] = rt
        tmp[[0,-1]] = np.mean(rt)
        
        '''
        interpolate NaN values using nearest neighbor interpolation and remove
        the padding.
        The resulting rt array is resampled to match the BOLD time series
        '''
        rt = pd.Series(tmp).interpolate(method='nearest',limit=4,
                      limit_direction='both').values.ravel()
        rt = rt[1:601]




        ## compute autocorrelation
        x = rt - np.mean(rt)
        norm = np.dot(x, x)

        ac = np.correlate(x,x,'full') / norm
        ac = ac[599:]
        ac_mean += (ac - ac_mean) / (i_subj + 1)
        plt.plot(ac)
    plt.plot(ac_mean,'k', linewidth=2)
