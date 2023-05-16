clc
clear
addpath('E:\Monkey\Superlets');
inputdir = 'E:\Monkey\ob';
outdir = 'E:\Monkey\BetaBurst\Left_M1\';
%Session = {'0521';'0524';'0528';'0530';'0601';'0605';'0607'};%7 sessions for monkry w
Session = {'0917';'0918';'0919';'0920';'0921';'0924';'0926';'0927';'1001';'1003';'1004';'1005'};%12 sessions for monkry ob

for s = 1:length(Session)
    
    load([inputdir,'\','lfp_trials_2019',num2str(Session{s})])
    
    n_trials=length(new_trials);
    
    
    fs = 1000;%sampling frequency
    fois = linspace(4,50,100);%frequency of interest
    cl = 4;
    srord = [1,40];
    mult = 1;
    
    for chan = 1:32 %% Left M1
        WT = {};
        for t_idx = 1:n_trials
            if new_trials(t_idx).reach_hand == 2 %%Right hand
                trial = new_trials(t_idx);
                Xsignal = trial.lfp(:,chan)';
                
                [wtresult] = faslt(Xsignal, fs, fois, cl, srord, 1);
                
                WT{t_idx} = wtresult;
            end
        end
        save(fullfile(outdir,strcat(Session{s},'_WT_chan',num2str(chan),'_LH')),'WT');
        
    end
end
                                                                                                                                                                                Xsignal={};
inputdir = 'E:\Monkey\Superlets\';
load('E:\Monkey\ob\lfp_trials_20190917.mat');

for s = 1:32
    
    load([inputdir,'\','WT_chan',num2str(s)]);
    n_trials = length(WT);
    
    Xsignal={};
    baseline_epoch=[];
    target_on_epoch=[];
    go_epoch=[];
    reach_on_epoch=[];
    enter_target_epoch=[];
    return_on_epoch=[];
    return_peak_epoch=[];
    
    % Time windows of interest (WOIs)
    baseline_woi=[-.5 0];
    t_onset_woi=[-.5 .5];
    go_woi=[-.5 .25];
    reach_on_woi=[-.25 .25];
    enter_targ_woi=[-.25 .25];
    return_on_woi=[-.25 .5];
    return_peak_woi = [-.25 .25];
    
    t_num = 0;
    for t_idx = 1:n_trials
        trial = new_trials(t_idx);
        if trial.reach_hand==2
            
            t_num = t_num+1;
            times = trial.lfp(:,end);
            
            Xsignal = WT{t_idx};
            target_on_center = knnsearch(times, trial.targ_onset-trial.start_time);
            
            baseline_win=[target_on_center+baseline_woi(1)*1000:target_on_center+baseline_woi(2)*1000];
            baseline_epoch = Xsignal(:,baseline_win);
            baseline_times=times(baseline_win);
            
            target_on_win=[target_on_center+t_onset_woi(1)*1000:target_on_center+t_onset_woi(2)*1000];
            target_on_times=times(target_on_win)-times(target_on_center);
            target_on_epoch = Xsignal(:,target_on_win);
            
            go_center = knnsearch(times, trial.go_cue-trial.start_time);
            go_win=[go_center+go_woi(1)*1000:go_center+go_woi(2)*1000];
            go_times=times(go_win)-times(go_center);
            go_epoch = Xsignal(:,go_win);
            
            reach_on_center = knnsearch(times, trial.reach_onset-trial.start_time);
            reach_on_win=[reach_on_center+reach_on_woi(1)*1000:reach_on_center+reach_on_woi(2)*1000];
            reach_on_times=times(reach_on_win)-times(reach_on_center);
            reach_on_epoch = Xsignal(:,reach_on_win);
            
            enter_target_center = knnsearch(times, trial.enter_targ-trial.start_time);
            enter_target_win=[enter_target_center+enter_targ_woi(1)*1000:enter_target_center+enter_targ_woi(2)*1000];
            enter_target_times=times(enter_target_win)-times(enter_target_center);
            enter_target_epoch = Xsignal(:,enter_target_win);
            
            return_on_center = knnsearch(times, trial.return_onset-trial.start_time);
            return_on_win=[return_on_center+return_on_woi(1)*1000:return_on_center+return_on_woi(2)*1000];
            return_on_times=times(return_on_win)-times(return_on_center);
            return_on_epoch = Xsignal(:,return_on_win);
            
            return_peak_center = knnsearch(times, trial.return_peak-trial.start_time);
            return_peak_win=[return_peak_center+return_peak_woi(1)*1000:return_peak_center+return_peak_woi(2)*1000];
            return_peak_times=times(return_peak_win)-times(return_peak_center);
            return_peak_epoch = Xsignal(:,return_peak_win);
            
            
            mean_baseline = squeeze(mean(baseline_epoch,2));
            corrected_target_on(:,:,t_num) = 10*log10(target_on_epoch./repmat(mean_baseline,1,size(target_on_epoch,2)));%(target_on_epoch-repmat(mean_baseline,1,size(target_on_epoch,2)))./repmat(mean_baseline,1,size(target_on_epoch,2))*100; %
            corrected_go(:,:,t_num) = 10*log10(go_epoch./repmat(mean_baseline,1,size(go_epoch,2)));%(go_epoch-repmat(mean_baseline,1,size(go_epoch,2)))./repmat(mean_baseline,1,size(go_epoch,2))*100;%
            corrected_reach_on(:,:,t_num) = 10*log10(reach_on_epoch./repmat(mean_baseline,1,size(reach_on_epoch,2)));%(reach_on_epoch-repmat(mean_baseline,1,size(reach_on_epoch,2)))./repmat(mean_baseline,1,size(reach_on_epoch,2))*100;%
            corrected_enter_target(:,:,t_num) = 10*log10(enter_target_epoch./repmat(mean_baseline,1,size(enter_target_epoch,2)));%(enter_target_epoch-repmat(mean_baseline,1,size(enter_target_epoch,2)))./repmat(mean_baseline,1,size(enter_target_epoch,2))*100;%
            corrected_return_on(:,:,t_num) = 10*log10(return_on_epoch./repmat(mean_baseline,1,size(return_on_epoch,2)));%(return_on_epoch-repmat(mean_baseline,1,size(return_on_epoch,2)))./repmat(mean_baseline,1,size(return_on_epoch,2))*100;%
            corrected_return_peak(:,:,t_num) = 10*log10(return_peak_epoch./repmat(mean_baseline,1,size(return_peak_epoch,2)));%(return_peak_epoch-repmat(mean_baseline,1,size(return_peak_epoch,2)))./repmat(mean_baseline,1,size(return_peak_epoch,2))*100;%
            clear baseline_epoch target_on_epoch go_epoch reach_on_epoch enter_target_epoch return_on_epoch return_peak_epoch
            corrected_baseline(:,t_num) = mean_baseline;
        end
    end
    
    mean_target_on=squeeze(mean(corrected_target_on,3));
    mean_go=squeeze(mean(corrected_go,3));
    mean_reach_on=squeeze(mean(corrected_reach_on,3));
    mean_enter_target=squeeze(mean(corrected_enter_target,3));
    mean_return_on=squeeze(mean(corrected_return_on,3));
    mean_return_peak=squeeze(mean(corrected_return_peak,3));
    
    nchans_baseline(s,:) = squeeze(mean(corrected_baseline,2)');
    nchans_target_on(s,:) =mean(mean_target_on(21:57,:),1);
    nchans_go(s,:)=mean(mean_go(21:57,:),1);
    nchans_reach_on(s,:)= mean(mean_reach_on(21:57,:),1);
    nchans_enter_target(s,:)=mean(mean_enter_target(21:57,:),1);
    nchans_return_on(s,:)= mean(mean_return_on(21:57,:),1);
    nchans_return_peak(s,:)= mean(mean_return_peak(21:57,:),1);%%average beta power
    
end
save nchans_LM1_0917_beta nchans_baseline nchans_target_on nchans_go nchans_reach_on nchans_enter_target nchans_return_on nchans_return_peak

all_tf = [nchans_target_on nchans_go nchans_reach_on nchans_enter_target nchans_return_on nchans_return_peak];
c1 = max(all_tf(:));c2 = min(all_tf(:));

load cmap;

figure();
subplot(1,16,[1:4]);
hold all;
for ii = 1:32
    plot(target_on_times,nchans_target_on(ii,:),'Color',cmap(ii,:),'LineWidth',1)
    hold on
end
hold off
xlim(target_on_times([1 end]));
ylim([c2 c1]);
xlabel('Time (s)');
ylabel('Beta power');
title('target onset');

subplot(1,16,[5 6 7]);
hold all;
for ii = 1:32
    plot(go_times,nchans_go(ii,:),'Color',cmap(ii,:),'LineWidth',1)
    hold on
end
hold off
xlim(go_times([1 end]));
ylim([c2 c1]);
xlabel('Time (s)');
ylabel('Beta power');
title('go');

subplot(1,16,[8 9]);
hold all;
for ii = 1:32
    plot(reach_on_times,nchans_reach_on(ii,:),'Color',cmap(ii,:),'LineWidth',1)
    hold on
end
hold off
xlim(reach_on_times([1 end]));
ylim([c2 c1]);
xlabel('Time (s)');
ylabel('Beta power');
title('reach onset');

subplot(1,16,[10 11]);
hold all;
for ii = 1:32
    plot(enter_target_times,nchans_enter_target(ii,:),'Color',cmap(ii,:),'LineWidth',1)
    hold on
end
hold off
xlim(enter_target_times([1 end]));
ylim([c2 c1]);
xlabel('Time (s)');
ylabel('Beta power');
title('enter target');

subplot(1,16,[12 13 14]);
hold all;
for ii = 1:32
    plot(return_on_times,nchans_return_on(ii,:),'Color',cmap(ii,:),'LineWidth',1)
    hold on
end
hold off
xlim(return_on_times([1 end]));
ylim([c2 c1]);
xlabel('Time (s)');
ylabel('Beta power');
title('return onset');

subplot(1,16,[15 16]);
hold all;
for ii = 1:32
    plot(return_peak_times,nchans_return_peak(ii,:),'Color',cmap(ii,:),'LineWidth',1)
    hold on
end
hold off
xlim(return_peak_times([1 end]));
ylim([c2 c1]);
xlabel('Time (s)');
ylabel('Beta power');
title('return peak');

pos=get(gca,'Position');
%colorbar();
set(gca,'Position', pos);
