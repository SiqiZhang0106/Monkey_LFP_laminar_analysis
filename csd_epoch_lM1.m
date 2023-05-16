clc
clear
%Session = {'0917';'0918';'0919';'0920';'0921';'0923';'0924';'0926';'0927';'1001';'1003';'1004';'1005';};
Session = {'0601';};
inputdir = 'E:\Monkey\w';
outdir = 'E:\Monkey\BetaBurst\individuals\';


addpath('E:\Monkey\burst_detection-main\matlab\');
%load('E:\Monkey\Superlets\WT_chan14.mat');%64+16for the right M1 region
addpath('E:\Monkey\Superlets\');

midline = 100;

fs = 1000;
win_size=.26;
wlen = round(win_size * fs);
half_wlen = round(wlen * .5);
fois = linspace(4,50,100);%frequency of interest
c1 = 4;
srord = [1,40];
mult = 1;
%
 burst_detect_chan=12;

 for s =1:length(Session)
     
     load([inputdir,'\','lfp_trials_2018',num2str(Session{s})])
     %load('E:\Monkey\ob\lfp_trials_20190927.mat');
     n_trials=length(new_trials);
     
     
     %
     % % Compute TF (or can just load it)
     WT={};
     %trial_psds=[];
     trial_right = 0;
    
     for t_idx=1:n_trials
         t_idx
         if new_trials(t_idx).reach_hand ==2
             trial_right = trial_right+1;
             trial = new_trials(t_idx);
             times = trial.lfp(:,end);
             %
             % Center electrode signal
             Xsignal = trial.lfp(:,burst_detect_chan)';
             %     %common average referencing
             %     Xsignalall = trial.lfp(:,1:32)';
             %     %Xsignal = Xsignalall-mean(Xsignalall);
             %     Xsignal = mean(Xsignalall);
             %     % Superlet TF
             [wtresult] = faslt(Xsignal, fs, fois, c1, srord, mult);
             %     %Average over time for trial PSD
             %trial_psds(t_idx,:)=mean(wtresult,2);
             %     %
             %     %WT_allchan{t_idx} = WT_allchan{t_idx}+WT{t_idx};
             %     trial_right = trial_right+1;
             %     trial_psds(trial_right,:)=mean(WT{t_idx},2)';
             WT{t_idx}=wtresult;
             clear wtresult;
             
         end
     end
     save(fullfile(outdir,strcat(Session{s},'_WT_chan12_RH')),'WT');
     %save WT_0927_chan14_RH WT%RH:Right Hand
 end

% 1/f
oof=log10(1./fois);
% Compute PSD by averaging trial PSDs
psd=mean(trial_psds);
% Fit 1/f to spectrum
lm_psd=fitlm(oof,psd,'RobustOpts','on');
% Compute FOOOF threshold
fooof_thresh=lm_psd.Fitted;

% Plot - sanity check
figure();
plot(fois,psd);
hold all;
plot(fois,fooof_thresh);

% Search from 10-33Hz
search_range=knnsearch(fois',[10 33]');
search_foi=fois(search_range(1):search_range(2));

bursts=[];
bursts.trial=[];
bursts.waveform=[];
bursts.peak_freq=[];
bursts.peak_amp_iter=[];
bursts.peak_amp_base=[];
bursts.peak_time=[];
bursts.peak_adjustment=[];
bursts.fwhm_freq=[];
bursts.fwhm_time=[];
bursts.polarity=[];
bursts.waveform_times=[];

lfp=[];
aligned_lfp=[];
all_b_idx=1;
trial_right = 0;
for t_idx=1:n_trials
    if new_trials(t_idx).reach_hand ==2
        trial_right = trial_right+1;
        
        t_idx
        wtresult=WT{t_idx};
        
        trial = new_trials(t_idx);
        times = trial.lfp(:,end);
        
        %common average referencing
        Xsignalall = trial.lfp(:,1:32)';
        %Xsignalall = Xsignalall-mean(Xsignalall);
        
        % Center electrode signal
        Xsignal = Xsignalall(burst_detect_chan,:);
        
        % Get bursts for this trial
        trial_bursts=extract_bursts_trial(Xsignal, wtresult(search_range(1):search_range(2),:),...
            times, search_foi, [10 33], fooof_thresh(search_range(1):search_range(2)), fs,...
            'win_size', .26);
        n_bursts=length(trial_bursts.peak_time);
        % Set trial index
        trial_bursts.trial=ones(1,n_bursts).*t_idx;
        
        % Add to list of all bursts
        bursts.trial(end+1:end+n_bursts)=trial_bursts.trial;
        bursts.waveform(end+1:end+n_bursts,:)=trial_bursts.waveform;
        bursts.peak_freq(end+1:end+n_bursts)=trial_bursts.peak_freq;
        bursts.peak_amp_iter(end+1:end+n_bursts)=trial_bursts.peak_amp_iter;
        bursts.peak_amp_base(end+1:end+n_bursts)=trial_bursts.peak_amp_base;
        bursts.peak_time(end+1:end+n_bursts)=trial_bursts.peak_time;
        bursts.peak_adjustment(end+1:end+n_bursts)=trial_bursts.peak_adjustment;
        bursts.fwhm_freq(end+1:end+n_bursts)=trial_bursts.fwhm_freq;
        bursts.fwhm_time(end+1:end+n_bursts)=trial_bursts.fwhm_time;
        bursts.polarity(end+1:end+n_bursts)=trial_bursts.polarity;
        bursts.waveform_times=trial_bursts.waveform_times;
        
        Xsignalall = Xsignalall-mean(Xsignalall);
        % For each burst in this trial
        for b=1:length(trial_bursts.peak_time)
            
            if trial_bursts.polarity(b)==0
                % Find peak time index
                peak_time=trial_bursts.peak_time(b);
                peak_idx=knnsearch(times,peak_time);
                
                for chan=1:32
                    % Extract LFP signals around peak time
                    chan_lfp=Xsignalall(chan,peak_idx - half_wlen:peak_idx + half_wlen);
                    % Correct for DC offset
                    %chan_lfp=chan_lfp-mean(chan_lfp);
                    % Flip if burst is flipped
                    %                 if trial_bursts.polarity(b)==1
                    %                     chan_lfp=-1*chan_lfp;
                    %                 end
                    
                    lfp(all_b_idx,chan,:)=chan_lfp;
                    aligned_lfp(all_b_idx,chan,:)=chan_lfp;
                    % If not the channel used to detect bursts, align waveform
                    %                 if chan==burst_detect_chan
                    %                     aligned_lfp(all_b_idx,chan,:)=chan_lfp;
                    %                 else
                    %                     % Bandpass filter within frequency span of burst
                    %                     freq_range = [trial_bursts.peak_freq(b)-.5*trial_bursts.fwhm_freq(b),...
                    %                         trial_bursts.peak_freq(b)+.5*trial_bursts.fwhm_freq(b)];
                    %
                    %                     dc=mean(chan_lfp);
                    %                     % Pad with 1s on either side
                    %                     padded_data=[repmat(dc, 1, fs) chan_lfp repmat(dc, 1, fs)];
                    %                     % Bandpass filter
                    %                     filtered = ft_preproc_bandpassfilter(padded_data, fs, freq_range, 6,...
                    %                         'but', 'twopass', 'reduce');
                    %                     % Remove padding
                    %                     filtered=filtered(fs+1:fs+length(chan_lfp));
                    %
                    %                     % Hilbert transform
                    %                     analytic_signal = hilbert(filtered);
                    %                     % Get phase
                    %                     instantaneous_phase = mod(unwrap(angle(analytic_signal)), pi);
                    %
                    %                     % Find phase local minima (near 0)
                    %                     [~,zero_phase_pts]= findpeaks(-1*instantaneous_phase);
                    %
                    %                     % Find local phase minima with negative deflection closest to TF peak
                    %                     center_idx=half_wlen;
                    %                     [~,min_idx]=min(abs(center_idx - zero_phase_pts));
                    %                     closest_pt = zero_phase_pts(min_idx);
                    %                     new_peak_time_idx = peak_idx+(center_idx - closest_pt);
                    %
                    %                     % If not cut off
                    %                     if new_peak_time_idx > half_wlen && new_peak_time_idx + half_wlen < size(trial.lfp,1)
                    %                         % Get aligned LFP
                    %                         aligned_chan_lfp=trial.lfp(new_peak_time_idx - half_wlen:new_peak_time_idx + half_wlen,chan);
                    %                         % Correct for DC offset
                    %                         aligned_chan_lfp=aligned_chan_lfp-mean(aligned_chan_lfp);
                    %                         % Flip if burst is flipped
                    %                         if trial_bursts.polarity(b)==1
                    %                             aligned_chan_lfp=-1*aligned_chan_lfp;
                    %                         end
                    %                         % Can't align - use raw LFP
                    %                     else
                    %                         aligned_chan_lfp=chan_lfp;
                    %                     end
                    %                     aligned_lfp(all_b_idx,chan,:)=aligned_chan_lfp;
                    %                 end
                end
                all_b_idx=all_b_idx+1;
            end
        end
    end
end
%save bursts_chan20 bursts
mean_lfp=squeeze(mean(lfp,1));
mean_aligned_lfp=squeeze(mean(aligned_lfp,1));
%save mean_aligned_lfp mean_aligned_lfp

plot_lfp(mean_lfp,bursts.waveform_times);
plot_lfp(mean_aligned_lfp,bursts.waveform_times);
save mean_aligned_lfp_chan14_test4 mean_aligned_lfp
save aligned_lfp_chan14_test4 aligned_lfp



mean_aligned_lfp = mean_aligned_lfp(:,41:161);

spacing=0.1*10^-3; %%%%%%%%%%% spacing between neiboring electrodes

nd=1;
csdnew = [];
for t = 1:size(mean_aligned_lfp,2)
    phi = mean_aligned_lfp(:,t);
    csdnew(1,t)=mean_aligned_lfp(1,t);
    csdnew(2,t)=mean_aligned_lfp(2,t);
    for z = 3:30
        csdnew(z,t)=(phi(z+2)-2*phi(z)+phi(z-2))/((nd*spacing)^2);
    end
    csdnew(31,t)=mean_aligned_lfp(31,t);
    csdnew(32,t)=mean_aligned_lfp(32,t);
end

csd.csd=csdnew;
csd.time=-0.06:0.001:0.06;
%csd.time=bursts.waveform_times(41:161)';
csd.label=1:1:32;
%csd.dimord='rpt_chan_time';
csd.dimord=[];
cfg = [];
%cfg.baseline = [-0.03 0];
addpath('E:\Monkey\CSD')
plot_csd(cfg, csd, mean_aligned_lfp)

timesnew=bursts.waveform_times(41:161);
[X,Y]=meshgrid(timesnew,[1:32]);
% Only smooth in spatial domain
[Xq,Yq]=meshgrid(timesnew,[1:.1:32]);

csdsmooth=interp2(X,Y,csdnew,Xq,Yq,'spline');
Iblur = imgaussfilt(csdnew, 3);
csdsmooth=Iblur;

figure
imagesc(timesnew,(1:.1:32)',csdsmooth, median(csdsmooth(:))+3*[-iqr(csdsmooth(:)) iqr(csdsmooth(:))]); %%%%%%%%%%%%%%%%%  fixing the color range for comparing different data
%yticklabels(num2cell(chnl_labels));
colormap(redblue); % blue = source; red = sink
xlabel(' time (s)');      title('smoothed CSD (\color{red}sink, \color{blue}source\color{black})');
colorbar