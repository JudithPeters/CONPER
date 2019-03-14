import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt

# %matplotlib qt

#%% prepare data

rm_subj = [0, 8]
rm_roi = [82, 83, 84]
TS = np.load('fDATA_7T_3DMC_HP_CSF_WD/TS.npy')
TS = np.delete(TS, rm_subj, axis = 0)
TS = np.delete(TS, rm_roi, axis = 2)
n_subjects, n_tasks, n_regions, n_time = np.shape(TS)


# bandpass filter
nyquist = 0.5 / 1.5
cutoff = np.array([0.01,0.1])
b, a = signal.butter(2, cutoff / nyquist, 'bandpass')
TS_filt = signal.filtfilt(b, a, TS, axis = 3)

# hilbert transform
TS_hilbert = signal.hilbert(TS_filt, axis = 3)


#%% compute DFC and var(DFC)

DFC = np.zeros((n_subjects, n_tasks, n_time, n_time))
DFC_var = np.zeros((n_subjects, n_tasks, n_time))
ones = np.ones(n_regions)
PHI = np.angle(TS_hilbert)

n_links = int(0.5 * n_regions * (n_regions - 1))
mask = np.triu(np.ones((n_regions,n_regions)), 1).astype(bool)
iFC = np.zeros((n_time, n_links))

for i_subj in range(n_subjects):
    for i_task in range(n_tasks):
        for t in range(n_time):
            delta_PHI = np.outer(PHI[i_subj, i_task, :, t], ones) -\
                        np.outer(ones, PHI[i_subj, i_task, :, t])
            iFC[t, :] = np.cos(delta_PHI[mask])
        
        
        iFC_mag = np.sqrt(np.sum(iFC ** 2,  axis = 1))
        iFC_mean = np.mean(iFC, axis = 0)
        DFC[i_subj, i_task, :, :] = np.dot(iFC,iFC.T) /\
                                    np.outer(iFC_mag, iFC_mag)
        DFC_var[i_subj, i_task, :] = np.mean((iFC - iFC_mean) ** 2 , axis = 1)


#%% autocorrelation
n_lags = 100
idx = range(n_time - 1, n_time - 1 + n_lags)
ac_mean = np.zeros((n_subjects, n_tasks, n_lags))
ac_var = np.zeros((n_subjects, n_tasks, n_lags))

for i_subj in range(n_subjects):
    for i_task in range(n_tasks):
        x = DFC_var[i_subj, i_task, :] - np.mean(DFC_var[i_subj, i_task, :])
        norm = np.sum(x ** 2)
        ac_var[i_subj, i_task, :] = (signal.correlate(x, x,'full') / norm)[idx]
        for t in range(n_time):
            dfc = np.roll(DFC[i_subj, i_task, t, :], -t)
            ac_mean[i_subj, i_task, :] += (dfc[0:n_lags] -\
                   ac_mean[i_subj, i_task, :]) / (t + 1)

plt.figure
plt.figure
for i in range(n_tasks):
    plt.figure(1)
    plt.subplot(2,2, i + 1)
    plt.plot([0, n_lags], [0, 0], 'r--', linewidth = 5)
    plt.plot(ac_mean[:, i, :].T)
    plt.plot(np.mean(ac_mean[:,i,:],axis=0),'k', linewidth = 10)
    
    plt.figure(2)
    plt.subplot(2,2, i + 1)
    plt.plot([0, n_lags], [0, 0], 'r--', linewidth = 5)
    plt.plot(np.mean(ac_var[:,i,:],axis=0),'k', linewidth = 10)
    plt.plot(ac_var[:, i, :].T)

