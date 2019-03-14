%% Data & Metrics
%% Data
%% Load & arrange data

tasks =  {'RS';...
    'MR';...
    'OD';...
    'NUM'};

load('fDATA_7T_3DMC_HP_CSF_WD.mat','TS')

n_regions = 82;             % number of ROIs (ignore the last 3 regions)
[n_subjects,...             % retrieve number of subjects,
    n_tasks,...             % number of tasks
    ~,...                   %
    n_time] = size(TS);     % and number of time points from size of TS

TS = permute(TS,[4,3,2,1]); % re-arrange dimensions of data
%% Band-pass filtering
% The BOLD signal of each region of interest (ROI) was bandpass filtered in 
% the range $\left\lbrack 0\ldotp 04\mathrm{Hz},0\ldotp 07\mathrm{Hz}\right\rbrack$ 
% using a second order Butterworth filter.

s_rate = 1.5;               % sampling rate (seconds)
s_freq = 1 / s_rate;        % sampling frequency (Hz)
n_freq = 0.5 * s_freq;      % nyquist frequency
c_high = 0.04;              % high-pass cutoff
c_low = 0.07;               % low-pass cutoff
band = [c_high,c_low] / n_freq;

[b,a] = butter(2,band,'bandpass');

TS_filt = zeros(n_time,...  % pre-allocate memory for filtered data
    n_regions,...
    n_tasks,...
    n_subjects);

for s=1:n_subjects
    for t=1:n_tasks
        TS_filt(:,:,t,s) = filtfilt(b,a,TS(:,1:n_regions,t,s));
    end
end
%% *Metrics*
%% *Grand average functional connectivity*
% This metric reflects the static component of the relationships between  region-specific 
% activation patterns. A grand average FC matrix per task  (and rest) was obtained 
% by calculating individual FC matrices in the  form of pairwise Pearson correlation 
% coefficients of bandpass-filtered  (in the range from .04 Hz to .07 Hz) regional 
% BOLD signals for each subject, and subsequently averaging across individual 
% subject's FCs. 

gFC = zeros(n_regions,...   % pre-allocate global FC matrix
    n_regions,...
    n_tasks);


for t=1:n_tasks             % compute FC per task
    c = 1;
    for s=1:n_subjects      % average over subjects
        ts = TS_filt(:,:,t,s);
        if ~isnan(ts(1,1))  % some data may be missing
            % compute mean using Welford's Online algorithm
            gFC(:,:,t) = gFC(:,:,t) + (corrcoef(ts) - gFC(:,:,t)) / c;
            c = c + 1;
        end
    end
end

fprintf('global functional connectivity')
figure
for i=1:4
    subplot(2,2,i)
    imagesc(gFC(:,:,i))
    title(tasks{i})
    if i==3
        xlabel('ROI index')
        ylabel('ROI index')
    end
end
colormap jet
%% *Dynamic functional connectivity*
% This metric reflects the dynamics of functional couplings as short-lived global 
% network states dissolve and may re-emerge at different moments in time. The 
% DFC for each task (and rest) was given by the cumulative distribution of the 
% pairwise similarity between instantaneous functional connectivity (iFC) matrices 
% across subjects. Specifically, in an individual subject we first computed the 
% analytic signal of each voxel cortical region by applying  the Hilbert transformation 
% to its bandpass-filtered BOLD signal. 

TS_analytic = zeros(n_time,... % pre-allocate memory for analytical signal
    n_regions,...
    n_tasks,...
    n_subjects);

for s=1:n_subjects
    for t=1:n_tasks
        TS_analytic(:,:,t,s) = hilbert(TS_filt(:,:,t,s));
    end
end
%% 
% This  allowed for the calculation of instantaneous phases (i.e. phases  
% observed at each moment in time) of cortical regions. 

PHI = angle(TS_analytic);
%% 
% Next, a phase difference matrix was obtained at each point in time by 
% calculating the instantaneous phase differences between all pairs of cortical 
% regions. 

% exemplary for first point in time of rest in subject 1
phase_diff = meshgrid(PHI(1,:,1,1)) - meshgrid(PHI(1,:,1,1))';

fprintf('phase difference')
figure
imagesc(phase_diff)
xlabel('ROI index')
ylabel('ROI index')
colormap jet
%% 
% These difference matrices were then transformed to similarity matrices 
% by computing the cosine of their entries. A single phase similarity matrix reflects 
% the functional connectivity among cortical regions observed at a single moment 
% in time; i.e. the instantaneous functional connectivity (iFC). 

% exemplary for first point in time of rest in subject 1
iFC = cos(phase_diff);

fprintf('instantaneous functional connectivity')
figure
imagesc(iFC)
xlabel('ROI index')
ylabel('ROI index')
colormap jet
%% 
% To estimate the similarity between functional connectivity observed at 
% different moments in time, we calculated the <https://en.wikipedia.org/wiki/Cosine_similarity 
% cosine similarity> of the upper triangular of iFC matrices between all pairs 
% of time points. This results in a DFC matrix. 

% exemplary for rest in subject 1
mask_fc = logical(triu(ones(n_regions),1)); % mask of upper triungular

DFC = zeros(n_time);
for t1=1:n_time-1
    phase_diff = meshgrid(PHI(t1,:,1,1)) - meshgrid(PHI(t1,:,1,1))';
    iFC1 = cos(phase_diff);
    iFC1 = iFC1(mask_fc);
    M1 = sqrt(sum(iFC1.^2));
    for t2=t1+1:n_time
        phase_diff = meshgrid(PHI(t2,:,1,1)) - meshgrid(PHI(t2,:,1,1))';
        iFC2 = cos(phase_diff);
        iFC2 = iFC2(mask_fc);
        M2 = sqrt(sum(iFC2.^2));
        DFC(t1,t2) = iFC1' * iFC2 / (M1 * M2);
    end
end
DFC = DFC + DFC' + eye(n_time);

fprintf('dynamic functional connectivity')
figure
imagesc(DFC)
xlabel('time')
ylabel('time')
colormap jet

%% 
% The entries falling in the upper triangular of this matrix form the distribution 
% of similarity among pairs of time points in terms of the functional connectivity 
% observed at these moments in time. 

mask_dfc = logical(triu(ones(n_time),1)); % mask of upper triungular
DFC_vec = DFC(mask_dfc);

fprintf('distribution of DFC values')
figure
hist(DFC_vec,100)
xlabel('DFC value')
ylabel('count')
%% 
% This procedure was repeated for all subjects with the final distribution 
% of cosine similarity values being the aggregated distributions observed for 
% individual subjects.

DFC_vec = cell(n_tasks,1);  % pre-allocated DFC distributions (they may have different number of entries per task)
DFC = zeros(n_time);
for i=1:n_tasks             % compute DFC distribution per task
    DFC_vec{i} = [];        % assign an empty vector
    for s=1:n_subjects      % average over subjects
        phi = PHI(:,:,i,s);
        if ~isnan(phi(1,1)) % some data may be missing
            for t1=1:n_time-1
                phase_diff = meshgrid(phi(t1,:)) - meshgrid(phi(t1,:))';
                iFC1 = cos(phase_diff);
                iFC1 = iFC1(mask_fc);
                M1 = sqrt(sum(iFC1.^2));
                for t2=t1+1:n_time
                    phase_diff = meshgrid(phi(t2,:)) - meshgrid(phi(t2,:))';
                    iFC2 = cos(phase_diff);
                    iFC2 = iFC2(mask_fc);
                    M2 = sqrt(sum(iFC2.^2));
                    DFC(t1,t2) = iFC1' * iFC2 / (M1 * M2);
                end
            end
        end
        DFC_vec{i} = [DFC_vec{i};DFC(mask_dfc)];
    end
end

fprintf('dynamic functional connectivity')
figure
for i=1:4
    subplot(2,2,i)
    hist(DFC_vec{i},100)
    title(tasks{i})
    if i==3
        xlabel('DFC value')
        ylabel('count')
    end
end
%% *Metastability*
% The final metric used here reflects the overall variability of network states 
% of the system; i.e. to what extent the system exhibits transient synchronization 
% dynamics. Metastability in each task (and rest) was measured as the standard 
% deviation of the Kuramoto order parameter observed over time. The  Kuramoto 
% order parameter $R\left(t\right)$ reflects the extent of synchronization exhibited 
% among brain regions at a specific moment in time and is given by
% 
% $$R\left(t\right)=\left|\frac{1}{N}\sum_{j=1}^N e^{i\phi_j \left(t\right)} 
% \right|$$
% 
% with $\phi_j \left(t\right)$ being the instantaneous phase of the bandpass-filtered 
% BOLD signal of region $j$ at time $t$ and $N$ the total number of brain regions.

% exemplary for rest in subject 1

R = zeros(n_time,1);
for t=1:n_time
    R(t) = abs(mean(exp(PHI(t,:,1,1) * 1i)));
end
MS = std(R);

fprintf('metastability in rest of subject 1  = %.2f', MS)
%% 
% Metastability was estimated for each subject separately and subsequently 
% averaged.

MS = nan(10,4);             % pre-allocated MS distributions (they may have different number of entries per task)

for i=1:n_tasks             % compute MS per task
    for s=1:n_subjects      % average over subjects
        phi = PHI(:,:,i,s);
        if ~isnan(phi(1,1)) % some data may be missing
            for t=1:n_time
                R(t) = abs(mean(exp(phi(t,:) * 1i)));
            end
            MS(s,i) = std(R);
        end
    end
end

fprintf('dynamic functional connectivity')
figure
boxplot(MS,tasks)

MS_average = nanmean(MS);