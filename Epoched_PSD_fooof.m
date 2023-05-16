clc
clear
Session = {'0917';'0918';'0919';'0920';'0921';'0924';'0926';'0927';'1001';'1003';'1004';'1005'};%12 sessions for monkey ob
%Session = {'0521';'0524';'0528';'0530';'0601';'0605';'0607'};%7 sessions
%for monkey w
inputdir = 'D:\Monkey\ob';
outdir1 = 'D:\Monkey\O_Left_M1\Baseline\';
outdir2 = 'D:\Monkey\O_Left_M1\Target_on\';
outdir3 = 'D:\Monkey\O_Left_M1\Go_cue\';
outdir4 = 'D:\Monkey\O_Left_M1\Reach_onset\';
outdir5 = 'D:\Monkey\O_Left_M1\Enter_target\';
outdir6 = 'D:\Monkey\O_Left_M1\Return_on\';
outdir7 = 'D:\Monkey\O_Left_M1\Return_peak\';

for s = 1:length(Session)
    
    load([inputdir,'\','lfp_trials_2019',num2str(Session{s})]);
    
    n_trials=length(new_trials);
    
    fs = 1000;%sampling frequency
    
    baseline_woi=[-.5 0];
    t_onset_woi = [-.5 .5];
    go_woi = [-.5 .25];
    reach_on_woi=[-.25 .25];
    enter_targ_woi=[-.25 .25];
    return_on_woi=[-.25 .5];
    return_peak_woi = [-.25 .25];
    
    WT_baseline = [];
    WT_targon = [];
    WT_gocue = [];
    WT_reachon = [];
    WT_entertarg = [];
    WT_returnon = [];
    WT_returnpeak = [];
    
    for t_idx = 1:n_trials

        trial = new_trials(t_idx);
        
        if ~isempty(trial.lfp)
            
            times = trial.lfp(:,end);
            
            Xsignal = trial.lfp(:,1:32)';%%left M1
            
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
            
            for chan = 1:32
                
                ch_dat_baseline = baseline_epoch(chan,:);%channels number
                [psdd1, freqs] = pwelch(ch_dat_baseline, 250, 125, 1024, fs);
                
                ch_dat_targon = target_on_epoch(chan,:);
                [psdd2, freqs] = pwelch(ch_dat_targon, 250, 125, 1024, fs);

                ch_dat_gocue = go_epoch(chan,:);
                [psdd3, freqs] = pwelch(ch_dat_gocue, 250, 125, 1024, fs);
                            
                ch_dat_reachon = reach_on_epoch(chan,:);%channels number
                [psdd4, freqs] = pwelch(ch_dat_reachon, 250, 125, 1024, fs);
                
                ch_dat_entertarg = enter_target_epoch(chan,:);%channels number
                [psdd5, freqs] = pwelch(ch_dat_entertarg, 250, 125, 1024, fs);
                
                ch_dat_returnon = return_on_epoch(chan,:);%channels number
                [psdd6, freqs] = pwelch(ch_dat_returnon, 250, 125, 1024, fs);
                
                ch_dat_returnpeak = return_peak_epoch(chan,:);%channels number
                [psdd7, freqs] = pwelch(ch_dat_returnpeak, 250, 125, 1024, fs);

                WT_baseline(t_idx,chan,:) = psdd1';
                WT_targon(t_idx,chan,:) = psdd2';
                WT_gocue(t_idx,chan,:) = psdd3';
                WT_reachon(t_idx,chan,:) = psdd4';
                WT_entertarg(t_idx,chan,:) = psdd5';
                WT_returnon(t_idx,chan,:) = psdd6';
                WT_returnpeak(t_idx,chan,:) = psdd7';
            end
        end
    end
    save(fullfile(outdir1,strcat(Session{s},'_PSD_lM1')),'WT_baseline');
    save(fullfile(outdir2,strcat(Session{s},'_PSD_lM1')),'WT_targon');
    save(fullfile(outdir3,strcat(Session{s},'_PSD_lM1')),'WT_gocue');
    save(fullfile(outdir4,strcat(Session{s},'_PSD_lM1')),'WT_reachon');
    save(fullfile(outdir5,strcat(Session{s},'_PSD_lM1')),'WT_entertarg');
    save(fullfile(outdir6,strcat(Session{s},'_PSD_lM1')),'WT_returnon');
    save(fullfile(outdir7,strcat(Session{s},'_PSD_lM1')),'WT_reachpeak');
    
    %     WT = WT_reachon(:,:,5:52);
    %     meanWT = squeeze(mean(WT));
    %     figure
    %     imagesc(meanWT);figure(gcf);
    %
    %     colormap(jet)
    %     colorbar;
    %
    %     %power_allsession(s,:,:)=meanWT;
    %
    %     WT=mean(meanWT(:,10:27),2);
    %     R = 1:1:32;
    %     figure; plot(WT,R)
    %     set(gca,'ydir','reverse');
end

for s = 1:length(Session)
    
    load([inputdir2,'\',num2str(Session{s}),'_PSD_lM1']);
    WT = WT_reachon; %% or any session
    
    for chan = 1:32
        
        % FOOOF settings
        settings = fooof_check_settings([]);  % Use defaults
        %settings.background_mode = 'knee';
        %settings.max_n_peaks = 3;
        
        f_range = [4,100];
        psd=squeeze(mean(WT(:,chan,:)));
        
        % Run FOOOF
        fooof_results = fooof(f, psd, f_range, settings, true);
        
        
        %         fooof_plot(fooof_results)
        %         figure;plot(fooof_results.freqs,  fooof_results.power_spectrum-fooof_results.ap_fit);
        %         r_squared(s,i) = fooof_results.r_squared;
        %         error(s,i) = fooof_results.error;
        bg_params(chan,:) = fooof_results.aperiodic_params;
        %         full(i,:) = fooof_results.power_spectrum;
        %         background(i,:) = fooof_results.ap_fit;
        relati(chan,:) = fooof_results.power_spectrum-fooof_results.ap_fit;
        freqs = fooof_results.freqs;
        
    end
    save(fullfile(outputdir,strcat(Session{s},'_baseline_fooof')),'bg_params');
    save(fullfile(outputdir,strcat(Session{s},'_baseline_relati')),'relati');
    
end