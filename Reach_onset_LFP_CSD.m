clc
clear
load('E:\Monkey\ob\lfp_trials_20190924.mat');

n_trials = length(new_trials);

Xsignal={};
baseline_epoch=[];
reach_onset_epoch=[];

% Time windows of interest (WOIs)
baseline_woi=[-.5 0];
reach_onset_woi=[-.25 .25];
t_right = 0;

for t_idx = 1:n_trials
    trial = new_trials(t_idx);
    
    if trial.reach_hand==2
        t_right = t_right+1;
        times = trial.lfp(:,end);
        
        Xsignal = trial.lfp(:,1:32)';
        
        target_on_center = knnsearch(times, trial.targ_onset-trial.start_time);
        reach_onset_center = knnsearch(times, trial.reach_onset-trial.start_time);
        
        baseline_win=[target_on_center+baseline_woi(1)*1000:target_on_center+baseline_woi(2)*1000];
        baseline_times=times(baseline_win);
        baseline_epoch(t_right,:,:) = Xsignal(:,baseline_win);
        
        reach_onset_win=[reach_onset_center+reach_onset_woi(1)*1000:reach_onset_center+reach_onset_woi(2)*1000];
        reach_onset_times=times(reach_onset_win)-times(reach_onset_center);
        lfp_reach_onset(t_right,:,:) = Xsignal(:,reach_onset_win);
    end
    %
    %     for i=1:size(lfp_reach_onset,1)
    %         corrected_lfp_reach_onset(i,:) = lfp_reach_onset(i,:) - nanmean(Xsignal(i,baseline_win),2);
    %     end
    %    reach_onset_epoch(t_idx,:,:) = corrected_lfp_reach_onset;
    
    
end

csd_baseline = [];
csd_reach_onset = [];


spacing=0.1*10^-3; %%%%%%%%%%% spacing between neiboring electrodes
nd=2;

%% baseline
mean_lfp=squeeze(mean(baseline_epoch));
csdnew = [];
for t = 1:size(mean_lfp,2)
    phi = mean_lfp(:,t);
%     csdnew(1,t)=mean_lfp(1,t);
%     csdnew(2,t)=mean_lfp(2,t);
    for z = 3:30
        csdnew(z,t)=(phi(z+2)-2*phi(z)+phi(z-2))/((nd*spacing)^2);
    end
%     csdnew(31,t)=mean_lfp(31,t);
%     csdnew(32,t)=mean_lfp(32,t);
end
csdnew(1:2,:) = [];
csd_baseline = csdnew;
clear mean_lfp csdnew



%% reach onset
mean_lfp=squeeze(mean(lfp_reach_onset));


csdnew = [];
for t = 1:size(mean_lfp,2)
    phi = mean_lfp(:,t);
%     csdnew(1,t)=mean_lfp(1,t);
%     csdnew(2,t)=mean_lfp(2,t);
    for z = 3:30
        csdnew(z,t)=(phi(z+2)-2*phi(z)+phi(z-2))/((nd*spacing)^2);
    end
%     csdnew(31,t)=mean_lfp(31,t);
%     csdnew(32,t)=mean_lfp(32,t);
end
csdnew(1:2,:) = [];

for i=1:size(csdnew,1)
    csd_reach_onset(i,:) = (csdnew(i,:) - nanmean(csd_baseline(i,:),2));
end

%% Plotting
mean_lfp = mean_lfp-mean(mean_lfp);
figure;plot_lfp(mean_lfp,reach_onset_times)
figure();

csd=[];
csd.csd= csd_reach_onset;
csd.time=reach_onset_times;
csd.label=1:1:28;
csd.dimord=[];
cfg = [];
plot_csd(cfg, csd, mean_lfp)
colormap(jet)
%colormap(redblue)