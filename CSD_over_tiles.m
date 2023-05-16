clc
clear
session = {'80to100';'80to90';'70to80';'60to70';'50to60';'40to50';'30to40';'20to30';'10to20';'0to10'};
%session = {'80to100';'60to80';'40to60';'20to40';'0to20'};
inputdir = 'E:\Monkey\New_CSD_withPCA';
addpath('E:\Monkey\CSD');

for s = 1:length(session)
    load([inputdir,'\','chan_lfp_PC6_',num2str(session{s}),'tile'])
%    load([inputdir,'\','session2_PC6_',num2str(session{s}),'tile'])
%     chan_lfp_enter_target(166:169,:,:)=[];
%     chan_lfp_enter_target(107:122,:,:)=[];
%     
%     chan_lfp_target_on(1527:1561,:,:)=[];
%     chan_lfp_target_on(1022:1230,:,:)=[];
%     chan_lfp_target_on(707:777,:,:)=[];
    
    chan_lfp_enter_target(629:650,:,:)=[];
    chan_lfp_enter_target(434:496,:,:)=[];
    chan_lfp_enter_target(341:354,:,:)=[];
    
    chan_lfp_target_on(2596:2681,:,:)=[];
    chan_lfp_target_on(1752:2070,:,:)=[];
    chan_lfp_target_on(1208:1361,:,:)=[];

      %%%average LFP among trials to get CSD
%     mean1=squeeze(mean(chan_lfp_target_on));
%     mean2=squeeze(mean(chan_lfp_enter_target));
%     mean_aligned_lfp = (mean1+mean2)/2;
%     spacing=0.1*10^-3; %%%%%%%%%%% spacing between neiboring electrodes
%     nd=2;
%    % mean_aligned_lfp = mean_aligned_lfp(:,71:191);
%     
%     csdnew = [];
%     for t = 1:size(mean_aligned_lfp,2)
%         phi = mean_aligned_lfp(:,t);
% %         csdnew(1,t)=mean_aligned_lfp(1,t);
% %         csdnew(2,t)=mean_aligned_lfp(2,t);
%         for z = 3:30
%             csdnew(z,t)=(phi(z+2)-2*phi(z)+phi(z-2))/((nd*spacing)^2);
%         end
% %         csdnew(31,t)=mean_aligned_lfp(31,t);
% %         csdnew(32,t)=mean_aligned_lfp(32,t);
%     end
%     csdnew(1:2,:) = [];
% 
%     for i=1:size(csdnew,1)
%         csdnew(i,:) = (csdnew(i,:) - nanmean(csdnew(i,:),2))./csd_var(i);
%     end
%     
%     
%     csd=[];
%     csd.csd= csdnew;
%     csd.time=[-0.06:0.001:0.06];
%     csd.label=1:1:28;
%     csd.dimord=[];
%     cfg = [];
%     plot_csd(cfg, csd, mean_aligned_lfp(3:30,71:191))
%     colormap(jet)   
%     
    
    %%Average trials after CSD calculation
    bursts_lfp = cat(1,chan_lfp_target_on,chan_lfp_enter_target);
    spacing=0.1*10^-3; %%%%%%%%%%% spacing between neiboring electrodes
    nd=2;
    
    for t_idx = 1:size(bursts_lfp,1)
        mean_lfp=squeeze(bursts_lfp(t_idx,:,:));%squeeze(mean(baseline_epoch));
        csdnew = [];
        for t = 1:size(mean_lfp,2)
            phi = mean_lfp(:,t);
            %     csdnew(1,t)=mean_lfp(1,t);
            %     csdnew(2,t)=mean_lfp(2,t);
            for z = 3:30
                csdnew(z-2,t)=(phi(z+2)-2*phi(z)+phi(z-2))/((nd*spacing)^2);
            end
            %     csdnew(31,t)=mean_lfp(31,t);
            %     csdnew(32,t)=mean_lfp(32,t);
        end
       
        %csdnew(1:2,:) = [];
        csd_trials(t_idx,:,:) = csdnew;clear csdnew
        
    end
    %save csd_trials_PC6_80to100 csd_trials
    
%     csd_mean = mean(squeeze(mean(csd_trials,3)));
%     csd_var = std(squeeze(mean(csd_trials,3)));
    
   
    for t_idx = 1:size(bursts_lfp,1)
        csdnew = squeeze(csd_trials(t_idx,:,:));
        for i=1:size(csdnew,1)
            csd_norm(i,:) = (csdnew(i,:)-mean(csdnew(i,:)))./max(abs(csdnew(i,:)));%(csdnew(i,:) - csd_mean(i))/csd_var(i);
        end
        csd_norm_trials(t_idx,:,:) = csd_norm;
    end
    save csd_norm_trials_PC6_80to100 csd_norm_trials
    
    csd_norm_trials = reshape(csd_norm_trials,4305,28*261);
    [H,P,CI,STATS] = ttest(csd_norm_trials);
    p = reshape(P,28,261);
    t = reshape(STATS.tstat,28,261);
    figure
    imagesc(squeeze(p(:,71:191)));figure(gcf);
    temp=[0 0.2];
    colormap(jet)
    colorbar;
    caxis(temp)
    colormap(redblue)
    
        figure
    imagesc(squeeze(p(:,71:191)));figure(gcf);
    temp=[0 0.2];
    %colormap(jet)
    colorbar;
    caxis(temp)
    colormap(redblue)
    set(gca,'yticklabel',{'7','12','17','22','27'});
    
    csd_norm_trials = reshape(csd_norm_trials,4305,28,261);
    csd_norm_average = squeeze(mean(csd_norm_trials));
    %csd_zscore_average = squeeze(mean(csd_zscore_trials));
    mean_lfp=squeeze(mean(bursts_lfp));
    %mean_lfp = mean_lfp-mean(mean_lfp);
    %mean_lfp = mean_lfp(3:30,:);

    csd=[];
    csd.csd= csd_norm_average(:,71:191);
    csd.time=[-0.06:0.001:0.06];
    csd.label=1:1:28;
    csd.dimord=[];
    cfg = [];
    plot_csd(cfg, csd, mean_lfp(:,71:191))
    colormap(redblue)
    
p_oneTailed = normcdf(csd_zscore_average);   
%     csd1 = csdnew(13:14,35:45);
%     csd2 = csdnew(19:20,35:45);
%     csd3 = csdnew(30,35:45);
%     csd4 = csdnew(5:6,55:65);
%     csd5 = csdnew(13:14,55:65);
%     csd6 = csdnew(19:20,55:65);
%     csd7 = csdnew(30,55:65);
%     csd8 = csdnew(13:14,75:85);
%     csd9 = csdnew(19:20,75:85);
%     csdmean1(s,1) = mean(mean(csd1,2));
%     csdmean2(s,1) = mean(mean(csd2,2));
%     csdmean3(s,1) = mean(mean(csd3,2));
%     csdmean4(s,1) = mean(mean(csd4,2));
%     csdmean5(s,1) = mean(mean(csd5,2));
%     csdmean6(s,1) = mean(mean(csd6,2));
%     csdmean7(s,1) = mean(mean(csd7,2));
%     csdmean8(s,1) = mean(mean(csd8,2));
%     csdmean9(s,1) = mean(mean(csd9,2));
% 
%     clear csdnew csd1 csd2 csd3 csd4 csd5 csd6 csd7 csd8 csd9
end


% csdmean_winall=[csdmean1,csdmean2,csdmean3,csdmean4,csdmean5,csdmean6,csdmean7,csdmean8,csdmean9];
% save csdmean_winall_pc11 csdmean_winall
