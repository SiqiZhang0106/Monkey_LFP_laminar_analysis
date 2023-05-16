clc
clear
load('D:\Monkey\ob\lfp_trials_20190917.mat');
load('E:\Monkey\Superlets\WT_chan80.mat');%64+16for the right M1 region

addpath('D:\Monkey\burst_detection-main\matlab\');
%addpath('D:\Monkey\Superlets\');
n_trials=length(new_trials);

midline = 100;

fs = 1000;
win_size=.2;
wlen = round(win_size * fs);
half_wlen = round(wlen * .5);
fois = linspace(4,50,100);%frequency of interest
c1 = 4;
srord = [1,40];
mult = 1;

% DO WE WANT 16? ITS IN THE MIDDLE, BUT IS IT IN THE GRANULAR LAYER?
burst_detect_chan=80;

% Compute TF (or can just load it)
%WT={};
trial_psds=[];
for t_idx=1:n_trials
%     trial = new_trials(t_idx);
%     times = trial.lfp(:,end);
%     
%     % Center electrode signal
%     Xsignal = trial.lfp(:,burst_detect_chan)';
%     % Superlet TF
%     [wtresult] = faslt(Xsignal, fs, fois, c1, srord, mult);
%     %Average over time for trial PSD
%     trial_psds(t_idx,:)=mean(wtresult,2);
%     
 trial_psds(t_idx,:)=mean(WT{t_idx},2)';
%     WT{t_idx}=wtresult;
%     clear wtresult;
end
%  cd('D:\Monkey\Superlets');
%  save WT_0918_chan76 WT
%  cd('D:\Monkey\BetaBurst');

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

for t_idx=1:n_trials
    t_idx
    wtresult=WT{t_idx};
    
    trial = new_trials(t_idx);
    times = trial.lfp(:,end);
    
    % Center electrode signal
    Xsignal = trial.lfp(:,burst_detect_chan)';
    
    % Get bursts for this trial
    trial_bursts=extract_bursts_trial(Xsignal, wtresult(search_range(1):search_range(2),:),...
        times, search_foi, [13 30], fooof_thresh(search_range(1):search_range(2)), fs,...
        'win_size', win_size);
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
    
    % For each burst in this trial
    for b=1:length(trial_bursts.peak_time)
        
        if trial_bursts.polarity(b)==0
            % Find peak time index
            peak_time=trial_bursts.peak_time(b);
            peak_idx=knnsearch(times,peak_time);
            
            %Interpolation for the flat electrodes
            allchan_lfp=trial.lfp(peak_idx - half_wlen:peak_idx + half_wlen,65:96);%rM1-65 until 96
            chan_lfp_interpolation = allchan_lfp;
            
            for tp = 1:size(allchan_lfp,1)%timepoint
                y = allchan_lfp(tp,:);
                x = 0:0.1:3.1;
                p=polyfit(x,y,2);
                figure;plot(x,y,'*',x,polyval(p,x))
                
                newlfp = [polyval(p,x(24)),polyval(p,x(25)),polyval(p,x(26)),polyval(p,x(27)),polyval(p,x(28)),polyval(p,x(29)),polyval(p,x(30)),polyval(p,x(31)),polyval(p,x(32))];
                chan_lfp_interpolation(tp,24:32) = newlfp;%25,26,27,28,31
               % chan_lfp_interpolation(tp,31) = polyval(p,x(31));%25,26,27,28,31
                clear p x y
            end
            
            for chan=65:96
                % Extract LFP signals around peak time
                %chan_lfp=trial.lfp(peak_idx - half_wlen:peak_idx + half_wlen,chan)';
                chan_lfp = chan_lfp_interpolation(:,chan-64)';
                
                % Correct for DC offset
                chan_lfp=chan_lfp-mean(chan_lfp);
                % Flip if burst is flipped
%                 if trial_bursts.polarity(b)==1
%                     chan_lfp=-1*chan_lfp;
%                 end
                
                lfp(all_b_idx,chan-64,:)=chan_lfp;
                
                % If not the channel used to detect bursts, align waveform
                if chan==burst_detect_chan
                    aligned_lfp(all_b_idx,chan-64,:)=chan_lfp;
                else
                    % Bandpass filter within frequency span of burst
                    freq_range = [trial_bursts.peak_freq(b)-.5*trial_bursts.fwhm_freq(b),...
                        trial_bursts.peak_freq(b)+.5*trial_bursts.fwhm_freq(b)];
                    
                    dc=mean(chan_lfp);
                    % Pad with 1s on either side
                    padded_data=[repmat(dc, 1, fs) chan_lfp repmat(dc, 1, fs)];
                    % Bandpass filter
                    filtered = ft_preproc_bandpassfilter(padded_data, fs, freq_range, 6,...
                        'but', 'twopass', 'reduce');
                    % Remove padding
                    filtered=filtered(fs+1:fs+length(chan_lfp));
                    
                    % Hilbert transform
                    analytic_signal = hilbert(filtered);
                    % Get phase
                    instantaneous_phase = mod(unwrap(angle(analytic_signal)), pi);
                    
                    % Find phase local minima (near 0)
                    [~,zero_phase_pts]= findpeaks(-1*instantaneous_phase);
                    
                    % Find local phase minima with negative deflection closest to TF peak
                    center_idx=half_wlen;
                    [~,min_idx]=min(abs(center_idx - zero_phase_pts));
                    closest_pt = zero_phase_pts(min_idx);
                    new_peak_time_idx = peak_idx+(center_idx - closest_pt);
                    
                    % If not cut off
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
                        aligned_chan_lfp=chan_lfp;
                    %end
                    aligned_lfp(all_b_idx,chan-64,:)=aligned_chan_lfp;
                end
            end
            all_b_idx=all_b_idx+1;
        end
    end
end
mean_lfp=squeeze(mean(lfp,1));
mean_aligned_lfp=squeeze(mean(aligned_lfp,1));
%save mean_aligned_lfp mean_aligned_lfp

plot_lfp(mean_lfp,bursts.waveform_times);
plot_lfp(mean_aligned_lfp,bursts.waveform_times);

mean_aligned_lfp = mean_aligned_lfp(:,41:161);

spacing=0.1*10^-3; %%%%%%%%%%% spacing between neiboring electrodes

nd=1;
csdnew = [];
for t = 1:size(mean_aligned_lfp,2)%32*201 channels*times
    phi = mean_aligned_lfp(:,t);
    csdnew(1,t)=mean_aligned_lfp(1,t);
    csdnew(2,t)=mean_aligned_lfp(1,t);
    for z = 3:22
        csdnew(z,t)=(phi(z+2)-2*phi(z)+phi(z-2))/((nd*spacing)^2);
    end
    csdnew(23,t)=mean_aligned_lfp(24,t);
    csdnew(24,t)=mean_aligned_lfp(24,t);
end

timesnew=bursts.waveform_times(41:161);
[X,Y]=meshgrid(timesnew,[1:24]);
% Only smooth in spatial domain
[Xq,Yq]=meshgrid(timesnew,[1:.1:24]);

csdsmooth=interp2(X,Y,csdnew,Xq,Yq,'spline');
Iblur = imgaussfilt(csdnew,1);
csdsmooth=Iblur;

figure
imagesc(timesnew,(1:.1:24)',csdsmooth, median(csdsmooth(:))+3*[-iqr(csdsmooth(:)) iqr(csdsmooth(:))]); %%%%%%%%%%%%%%%%%  fixing the color range for comparing different data
%yticklabels(num2cell(chnl_labels));
colormap(redblue); % blue = source; red = sink
xlabel(' time (s)');      title('smoothed CSD (\color{red}sink, \color{blue}source\color{black})');
colorbar

