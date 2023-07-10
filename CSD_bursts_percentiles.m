clc
clear
inputdir = 'E:\Monkey\New_CSD_withPCA';
outputdir = 'E:\Monkey\New_CSD_withPCA';
fs = 1000;
win_size=.26;
wlen = round(win_size * fs);
half_wlen = round(wlen * .5);
load PCA_scores_bursts.mat
load PCA_textdata.mat
load monkey_waveforms.mat
load polarity
for pc = 6
    for percentile = 1:10
        chan_lfp_enter_target = [];
        chan_lfp_target_on = [];
        waveform_enter_target = [];
        waveform_target_on = [];
        
        Spc = M(:,pc);
        [I,J]=find(Spc>prctile(Spc,10*(percentile-1))&Spc<=prctile(Spc,10*percentile));
        bursts_enter_target = 0;
        bursts_target_on = 0;
        
        for row = 1:length(I)
            if polarity(I(row))==0
                if length(textdata{I(row)+1,8})==3 %check if the session number begin with 0
                    load([inputdir,'\','Monkey_',textdata{I(row)+1,7},'\','epochs_',textdata{I(row)+1,9},'_session0',textdata{I(row)+1,8}])
                else if length(textdata{I(row)+1,8})==4
                        load([inputdir,'\','Monkey_',textdata{I(row)+1,7},'\','epochs_',textdata{I(row)+1,9},'_session',textdata{I(row)+1,8}])
                    end
                end
                
                peak_idx = round((str2num(textdata{I(row)+1,5})-(-0.5))*1000+1);%extract burst peak time
                trial = str2num(textdata{I(row)+1,6});%extract trial number
                
                file = textdata{I(row)+1,9};
                if ~isempty(strfind(file,'enter_target'))
                    bursts_enter_target = bursts_enter_target + 1;
                    Xsignal = squeeze(enter_target_epoch(trial,:,:));
                    clear enter_target_epoch
                    chan_lfp_enter_target(bursts_enter_target,:,:) = Xsignal(:,peak_idx - half_wlen:peak_idx + half_wlen);
                    waveform_enter_target(bursts_enter_target,:) = waveforms(I(row),:);
                else if ~isempty(strfind(file,'target_on'))
                        bursts_target_on = bursts_target_on + 1;
                        Xsignal = squeeze(target_on_epoch(trial,:,:));
                        clear target_on_epoch
                        chan_lfp_target_on(bursts_target_on,:,:) = Xsignal(:,peak_idx - half_wlen:peak_idx + half_wlen);
                        waveform_target_on(bursts_target_on,:) = waveforms(I(row),:);
                    end
                end
            end
            
        end
        save(fullfile(outputdir,strcat('chan_lfp_PC',num2str(pc),'_',num2str(10*(percentile-1)),'to',num2str(10*(percentile)),'percentile')),'chan_lfp_enter_target','chan_lfp_target_on');
        save(fullfile(outputdir,strcat('waveform_PC',num2str(pc),'_',num2str(10*(percentile-1)),'to',num2str(10*(percentile)),'percentile')),'waveform_enter_target','waveform_target_on')
%         save chan_lfp_PC9_90to100tile chan_lfp_enter_target chan_lfp_target_on
%         save waveform_PC9_90to100tile waveform_enter_target waveform_target_on
    end
end

clc
clear
session = {'90to100';'80to90';'70to80';'60to70';'50to60';'40to50';'30to40';'20to30';'10to20';'0to10'};
inputdir = 'E:\Monkey\New_CSD_withPCA';
addpath('E:\Monkey\CSD');

for pc = 6
    for s = 1:length(session)
        load([inputdir,'\','chan_lfp_PC',num2str(pc),'_',num2str(session{s}),'percentile'])
        
        %%Average trials before CSD calculation
        bursts_lfp = cat(1,chan_lfp_target_on,chan_lfp_enter_target);
        spacing=0.1*10^-3; %%%%%%%%%%% spacing between neiboring electrodes
        nd=2;        
        mean_lfp=squeeze(mean(bursts_lfp));
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
        
        csd_bursts(t_idx,:,:) = csdnew;
    end
end