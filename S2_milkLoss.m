%% S2_milkLoss
% This script aims at calculating milk losses and aswering research
% question about the milk losses detected. It starts from the different
% case datasets. Losses are calculated either via milk/quarter losses or
% via daily losses, dependent on what data are available.
% We have 3 types of case detection / datasets:
%       MPR data (milk recording, DHI)
%       Treatment registers
%       QL detection
%
% Research questions answered:
%   1: TR   treatment registers
%   2: MR   milk recording
%   3: QD   quarter dynamics
%
%
%
% 
%% LOAD DATA
clear variables 
close all
clc

cd C:\Users\u0084712\Documents\MastiMan\Research\mastitisMilkLoss


%% MODEL CHARACTERISTICS / difference in expected production @ QL

% differences quarter productions / 

%% TREATMENT REGISTERS
% number of cases matched with robot cow ID per farm
% number of cases for which lactation info is avaible at time of detection
%   and with enough info to estimate losses at DAY level
% number of cases for which lactation info is avaible at time of detection
%   and with enough info to estimate losses at MILK level
% milk production residuals 5 days before treatment
%       at milk level
%       at day level
%           per farm
%           per parity
%     ???   per lactation stage   ???
% proportion of cases with no loss before detection
% milk production residuals ??? 30 ??? days after treatment
%       at milk level
%       at day level
%           per farm
%           per parity
%     ???   per lactation stage   ???
% perturbation-based
%   number of cases associated with detected perturbation
%   start of the perturbation before treatment
%   end of the perturbation after treatment
%   losses per quarter
%   losses of quarter with highest relative loss compared to other losses
%       in the same time period
%   average difference between quarter losses in the same time periods
%       around a treatment
% summaries of losses and results across farms

% load day data via S1 case_selection > moddataDAY
load('D9_TR_cases.mat')
% >> 1 << Number of cases matched with 'complete' lactation data
farms = fieldnames(CLMAS);  % farms with CM data
ML.TR = table(zeros(0,1),'VariableNames',{'B_ID'}); % prepare case table
ML.TR.FarmName = repmat({''},0);
ML.TR = movevars(ML.TR, 'FarmName','Before','B_ID');
mod.Wood = @(p,t) p(1).*t.^p(2).*exp(-p(3).*t);     % wood model - prepare
T = 1;  % counter to fill table
for i = 1:length(farms)
    
    % first select cases not too close to eachother (at least 10 days
    % apart)
    data = CLMAS.(farms{i})(~isnan(CLMAS.(farms{i}).B_ID),:);
    data.diff(1,1) = 10;
    data.diff(2:end,1) = diff(datenum(data.Date));
    [~,indx] = unique(data.B_ID); % find first occurence
    data.first(indx) = 1;
    data(data.first == 0 & data.diff < 10,:) = []; % delete
    % re-add to CLMAS
    data = removevars(data,{'diff','first'});
    CLMAS.(farms{i}) = data; 
    clear data
    
    ind = find(~isnan(CLMAS.(farms{i}).B_ID));
        
    L = length(ind); % 
    ML.TR.FarmName(T:T+L-1,1) = repmat(farms(i),L,1);
    ML.TR.B_ID(T:T+L-1) = CLMAS.(farms{i}).B_ID(ind);
    ML.TR.Date(T:T+L-1) = CLMAS.(farms{i}).Date(ind);
    for ii = 1:length(ind)
        % find corresponding lactation and DIM
        try
            % find the corresponding data in a lactation within 305 dim
            idx = find(moddataDAY.(farms{i}).B_ID == ML.TR.B_ID(T+ii-1) & ...
                       floor(datenum(moddataDAY.(farms{i}).Date)) == floor(datenum(ML.TR.Date(T+ii-1))));     
        catch
            % fill in empty
            idx = [];
        end
        % fill in
        if ~isempty(idx)
            
            ML.TR.Lac(T+ii-1) = moddataDAY.(farms{i}).Lac(idx(1));
            ML.TR.DIM(T+ii-1) = moddataDAY.(farms{i}).DIM(idx(1));
            ML.TR.Q(T+ii-1) = CLMAS.(farms{i}).ms(ind(ii));
            ML.TR.Month(T+ii-1) = month(moddataDAY.(farms{i}).Date(idx(1)));
           
            % data to calculate losses
            data = moddataDAY.(farms{i})(moddataDAY.(farms{i}).B_ID ==  ML.TR.B_ID(T+ii-1) & ...
                                         moddataDAY.(farms{i}).Lac ==  ML.TR.Lac(T+ii-1),[2 8 9 10 13 14]);
                                     
            % summary expected lactation 305
            ML.TR.LacLength(T+ii-1) = max(data.DIM);   % lactation length
            ML.TR.PeakMYexp(T+ii-1) = max(data.ModIT); % expected peak yield
            ML.TR.TotMYexp(T+ii-1) = sum(data.ModIT);  % expected total yield
            % estimate 305 day MY
            p = modelDAY.(farms{i}).pFin{modelDAY.(farms{i}).pFin.B_ID ==  ML.TR.B_ID(T+ii-1) & ...
                                         modelDAY.(farms{i}).pFin.Lac ==  ML.TR.Lac(T+ii-1),3:5};
            ML.TR.TotMYexp305(T+ii-1) = sum(mod.Wood(p,0:305));
            
            % prepare milk loss calculations
            data.RES = data.MYs-data.ModIT;
            data.RELRES = 100*(data.RES./data.ModIT);
            data.RELRES(data.RELRES>100) = -0.1;
            
            % calculate losses 5 days before the detection
            if ML.TR.DIM(T+ii-1) > 10
                sub = data(data.DIM >= ML.TR.DIM(T+ii-1)-5 & data.DIM < ML.TR.DIM(T+ii-1),:);
                ML.TR.ML5DaysBefore(T+ii-1) = sum(sub.RES);     % total milk loss kg
                ML.TR.ML5DaysBeforeREL(T+ii-1) = mean(sub.RES); % average relative loss
            else
                ML.TR.ML5DaysBefore(T+ii-1) = NaN;   %
                ML.TR.ML5DaysBeforeREL(T+ii-1) = NaN; % average relative loss            
            end
            
            % calculate losses 21 days after detection
            sub = data(data.DIM >= ML.TR.DIM(T+ii-1) & data.DIM < ML.TR.DIM(T+ii-1)+21,:);
            ML.TR.ML21NoDays(T+ii-1) = height(sub);         % no of days

            if sub.DIM(end) < 305 && height(sub) < 21
                DIM = table((ML.TR.DIM(T+ii-1):ML.TR.DIM(T+ii-1)+20)','VariableNames',{'DIM'}); % DIM
                p = modelDAY.(farms{i}).pFin{modelDAY.(farms{i}).pFin.B_ID == sub.B_ID(1) & ...
                                             modelDAY.(farms{i}).pFin.Lac == sub.Lac(1),3:5};
                ITW = mod.Wood(p,DIM{:,1}); % ITW
                sub = outerjoin(sub,DIM,'MergeKeys',1);
                sub.MYs(isnan(sub.MYs)) = 0;
                sub.ModIT(:,1) = ITW;
                sub.RES = sub.MYs-sub.ModIT;
                sub.RELRES = 100*(sub.RES./sub.ModIT);
                sub.RELRES(sub.RELRES>100) = -0.1;
            end
            % calculate losses
            ML.TR.ML21DaysAfter(T+ii-1) = sum(sub.RES);     % total milk loss kg
            ML.TR.ML21DaysAfterREL(T+ii-1) = mean(sub.RES); % average relative loss
            
            
            % calculate % loss from day -5 to day 30
            for iii = -5:30
                loss = data.RELRES(data.DIM == ML.TR.DIM(T+ii-1)+iii);
                try
                    ML.TR.loss(T+ii-1,iii+6) = loss;
                catch
                    ML.TR.loss(T+ii-1,iii+6) = NaN;
                end
            end
            
            
        else
            ML.TR.Lac(T+ii-1) = NaN;
            ML.TR.DIM(T+ii-1) = NaN;
        end
    end
    T = T+L;
end
clear sub T L i ii idx ind data DIM ITW ans p test h

disp(['No. of unmatched cases = ' num2str(length(find(isnan(ML.TR.Lac)))) ...
       ', = ' num2str(length(find(isnan(ML.TR.Lac)))./height(ML.TR)*100) '%'])
disp(['No. of matched cases = ' num2str(length(find(~isnan(ML.TR.Lac)))) ...
       ', = ' num2str(length(find(~isnan(ML.TR.Lac)))./height(ML.TR)*100) '%'])

% delete unmatched cases
ML.TR = ML.TR(~isnan(ML.TR.DIM),:);

% summarize losses in plots
figure; subplot(1,3,1); histogram(ML.TR.Lac,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.7 0 0],...
                        'FaceColor',[0.7 0 0])
        title('Parity'); xlabel('Parity number')
        subplot(1,3,2); histogram(ML.TR.DIM,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0 0.7 0],...
                        'FaceColor',[0 0.7 0])
        title('Lactation stage'); xlabel('DIM (days)')
        subplot(1,3,3); histogram(ML.TR.Month,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.4 0 0.8],...
                        'FaceColor',[0.4 0 0.8])
        title('Month'); xlabel('Month')
disp(['No. of cases in M6/7/8 cases = ' num2str(length(find(ismember(ML.TR.Month,[6 7 9])))) ...
       ', = ' num2str(length(find(ismember(ML.TR.Month,[6 7 9])))./height(ML.TR)*100) '%'])
disp(['% of cases in first cases = ' num2str(length(find(ismember(ML.TR.Month,[6 7 9])))) ...
       ', = ' num2str(length(find(ismember(ML.TR.Month,[6 7 9])))./height(ML.TR)*100) '%'])
        
        
% visualize milk losses
ind11 = find(ML.TR.Lac == 1 & ML.TR.DIM <= 63);
ind12 = find(ML.TR.Lac == 1 & ML.TR.DIM > 63 & ML.TR.DIM <= 138);
ind13 = find(ML.TR.Lac == 1 & ML.TR.DIM > 138 & ML.TR.DIM <= 216);
ind14 = find(ML.TR.Lac == 1 & ML.TR.DIM > 216);
ind21 = find(ML.TR.Lac > 1 & ML.TR.DIM <= 63);
ind22 = find(ML.TR.Lac > 1 & ML.TR.DIM > 63 & ML.TR.DIM <= 138);
ind23 = find(ML.TR.Lac > 1 & ML.TR.DIM > 138 & ML.TR.DIM <= 216);
ind24 = find(ML.TR.Lac > 1 & ML.TR.DIM > 216);

figure('Units','Normalized','OuterPosition',[0.1 0.1 0.85 0.8]); 
subplot(2,4,1); histogram(ML.TR.ML5DaysBefore(ind11),'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0.8 1],...
                        'FaceColor',[0.2 0.8 1])
                hold on; plot([0 0], [0 0.4],'r-.','LineWidth',1.4);xlim([-150 60]);
        title('Lac = 1, LS = 1'); xlabel('ML 5 days before detection (kg)')
        subplot(2,4,2); xlim([-150 60]);histogram(ML.TR.ML5DaysBefore(ind12),'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0.6 1],...
                        'FaceColor',[0.2 0.6 1])
                hold on; plot([0 0], [0 0.4],'r-.','LineWidth',1.4)
                    title('Lac = 1, LS = 2'); xlabel('ML 5 days before detection (kg)');xlim([-150 60]);
        subplot(2,4,3);xlim([-150 60]); histogram(ML.TR.ML5DaysBefore(ind13),'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0 1],...
                        'FaceColor',[0.2 0 1])                
                    hold on; plot([0 0], [0 0.4],'r-.','LineWidth',1.4)
        title('Lac = 1, LS = 3'); xlabel('ML 5 days before detection (kg)');xlim([-150 60]);
        subplot(2,4,4); xlim([-150 60]);histogram(ML.TR.ML5DaysBefore(ind14),'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0 0 0.6],...
                        'FaceColor',[0 0 0.6])
                    hold on; plot([0 0], [0 0.4],'r-.','LineWidth',1.4)
        title('Lac = 1, LS = 4'); xlabel('ML 5 days before detection (kg)');xlim([-150 60]);
        subplot(2,4,5); xlim([-150 60]);histogram(ML.TR.ML5DaysBefore(ind21),'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0.8 1],...
                        'FaceColor',[0.2 0.8 1])
                hold on; plot([0 0], [0 0.4],'r-.','LineWidth',1.4)
        title('Lac > 1, LS = 1'); xlabel('ML 5 days before detection (kg)');xlim([-150 60]);
        subplot(2,4,6); xlim([-150 60]);histogram(ML.TR.ML5DaysBefore(ind22),'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0.6 1],...
                        'FaceColor',[0.2 0.6 1])
                hold on; plot([0 0], [0 0.4],'r-.','LineWidth',1.4)
        title('Lac > 1, LS = 2'); xlabel('ML 5 days before detection (kg)');xlim([-150 60]);
        subplot(2,4,7); xlim([-150 60]);histogram(ML.TR.ML5DaysBefore(ind23),'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0 1],...
                        'FaceColor',[0.2 0 1])                
                    hold on; plot([0 0], [0 0.4],'r-.','LineWidth',1.4)
        title('Lac > 1, LS = 3'); xlabel('ML 5 days before detection (kg)');xlim([-150 60]);
        subplot(2,4,8); xlim([-150 60]);histogram(ML.TR.ML5DaysBefore(ind24),'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0 0 0.6],...
                        'FaceColor',[0 0 0.6])
                    hold on; plot([0 0], [0 0.4],'r-.','LineWidth',1.4)
        title('Lac > 1, LS = 4'); xlabel('ML 5 days before detection (kg)');xlim([-150 60]);
                    
figure('Units','Normalized','OuterPosition',[0.1 0.1 0.85 0.8]);                   
        subplot(2,4,1); histogram(ML.TR.ML5DaysBeforeREL(ind11),'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0.8 1],...
                        'FaceColor',[0.2 0.8 1])
                hold on; plot([0 0], [0 0.5],'r-.','LineWidth',1.4); ylim([0 0.5])
        title('Lac = 1, LS = 1'); xlabel('ML 5 days before detection (%)');xlim([-30 10]);
        subplot(2,4,2);histogram(ML.TR.ML5DaysBeforeREL(ind12),'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0.6 1],...
                        'FaceColor',[0.2 0.6 1])
                hold on; plot([0 0], [0 0.5],'r-.','LineWidth',1.4);ylim([0 0.5])
        title('Lac = 1, LS = 2'); xlabel('ML 5 days before detection (%)');xlim([-30 10]);
        subplot(2,4,3);histogram(ML.TR.ML5DaysBeforeREL(ind13),'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0 1],...
                        'FaceColor',[0.2 0 1])                
                    hold on; plot([0 0], [0 0.5],'r-.','LineWidth',1.4);ylim([0 0.5])
        title('Lac = 1, LS = 3'); xlabel('ML 5 days before detection (%)');xlim([-30 10]);
        subplot(2,4,4);histogram(ML.TR.ML5DaysBeforeREL(ind14),'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0 0 0.6],...
                        'FaceColor',[0 0 0.6])
                    hold on; plot([0 0], [0 0.5],'r-.','LineWidth',1.4);ylim([0 0.5])
        title('Lac = 1, LS = 4'); xlabel('ML 5 days before detection (%)');xlim([-30 10]);
        subplot(2,4,5); histogram(ML.TR.ML5DaysBeforeREL(ind21),'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0.8 1],...
                        'FaceColor',[0.2 0.8 1])
                hold on; plot([0 0], [0 0.4],'r-.','LineWidth',1.4);ylim([0 0.30])
        title('Lac > 1, LS = 1'); xlabel('ML 5 days before detection (%)');xlim([-30 10]);
        subplot(2,4,6);histogram(ML.TR.ML5DaysBeforeREL(ind22),'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0.6 1],...
                        'FaceColor',[0.2 0.6 1])
                hold on; plot([0 0], [0 0.4],'r-.','LineWidth',1.4);ylim([0 0.30])
        title('Lac > 1, LS = 2'); xlabel('ML 5 days before detection (%)');xlim([-30 10]);
        subplot(2,4,7); histogram(ML.TR.ML5DaysBeforeREL(ind23),'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0 1],...
                        'FaceColor',[0.2 0 1])                
                    hold on; plot([0 0], [0 0.4],'r-.','LineWidth',1.4);ylim([0 0.30])
        title('Lac > 1, LS = 3'); xlabel('ML 5 days before detection (%)');xlim([-30 10]);
        subplot(2,4,8);histogram(ML.TR.ML5DaysBeforeREL(ind24),'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0 0 0.6],...
                        'FaceColor',[0 0 0.6])
                    hold on; plot([0 0], [0 0.4],'r-.','LineWidth',1.4);ylim([0 0.30])
        title('Lac > 1, LS = 4'); xlabel('ML 5 days before detection (%)');xlim([-30 10]);
                    
figure('Units','Normalized','OuterPosition',[0.1 0.1 0.85 0.8]);
binedges = -800:100:100;
        subplot(2,4,1); histogram(ML.TR.ML21DaysAfter(ind11),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0.8 1],...
                        'FaceColor',[0.2 0.8 1])
        subplot(2,4,2); histogram(ML.TR.ML21DaysAfter(ind12),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0.6 1],...
                        'FaceColor',[0.2 0.6 1])
        subplot(2,4,3); histogram(ML.TR.ML21DaysAfter(ind13),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0 1],...
                        'FaceColor',[0.2 0 1])                
        subplot(2,4,4); histogram(ML.TR.ML21DaysAfter(ind14),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0 0 0.6],...
                        'FaceColor',[0 0 0.6])
        subplot(2,4,5); histogram(ML.TR.ML21DaysAfter(ind21),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0.8 1],...
                        'FaceColor',[0.2 0.8 1])
        subplot(2,4,6); histogram(ML.TR.ML21DaysAfter(ind22),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0.6 1],...
                        'FaceColor',[0.2 0.6 1])
        subplot(2,4,7); histogram(ML.TR.ML21DaysAfter(ind23),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0 1],...
                        'FaceColor',[0.2 0 1])                
        subplot(2,4,8); histogram(ML.TR.ML21DaysAfter(ind24),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0 0 0.6],...
                        'FaceColor',[0 0 0.6])
    
pl.titles = {'Lac = 1, LS = 1','Lac = 1, LS = 2','Lac = 1, LS = 3','Lac = 1, LS = 4',...
             'Lac > 1, LS = 1','Lac > 1, LS = 2','Lac > 1, LS = 3','Lac > 1, LS = 4'};
pl.xlabel = {'sum ML 21 days after detection (kg)'};      
% pl.xlim = [-800 100];
pl.ylim = [0 0.8];
for i = 1:8
    subplot(2,4,i);
    title(pl.titles{i})
%     xlim(pl.xlim)
    ylim(pl.ylim)
    xlabel(pl.xlabel)
    hold on; plot([0 0], [0 1],'r-.','LineWidth',1.4);%xlim([-150 60]);
end
        
        
figure('Units','Normalized','OuterPosition',[0.1 0.1 0.85 0.8]);
binedges = -50:5:5;
        subplot(2,4,1); histogram(ML.TR.ML21DaysAfterREL(ind11),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0.8 1],...
                        'FaceColor',[0.2 0.8 1])
        subplot(2,4,2); histogram(ML.TR.ML21DaysAfterREL(ind12),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0.6 1],...
                        'FaceColor',[0.2 0.6 1])
        subplot(2,4,3); histogram(ML.TR.ML21DaysAfterREL(ind13),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0 1],...
                        'FaceColor',[0.2 0 1])                
        subplot(2,4,4); histogram(ML.TR.ML21DaysAfterREL(ind14),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0 0 0.6],...
                        'FaceColor',[0 0 0.6])
        subplot(2,4,5); histogram(ML.TR.ML21DaysAfterREL(ind21),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0.8 1],...
                        'FaceColor',[0.2 0.8 1])
        subplot(2,4,6); histogram(ML.TR.ML21DaysAfterREL(ind22),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0.6 1],...
                        'FaceColor',[0.2 0.6 1])
        subplot(2,4,7); histogram(ML.TR.ML21DaysAfterREL(ind23),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0 1],...
                        'FaceColor',[0.2 0 1])                
        subplot(2,4,8); histogram(ML.TR.ML21DaysAfterREL(ind24),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0 0 0.6],...
                        'FaceColor',[0 0 0.6])
    
pl.titles = {'Lac = 1, LS = 1','Lac = 1, LS = 2','Lac = 1, LS = 3','Lac = 1, LS = 4',...
             'Lac > 1, LS = 1','Lac > 1, LS = 2','Lac > 1, LS = 3','Lac > 1, LS = 4'};
pl.xlabel = {'average ML 21 days after detection %'};      
pl.xlim = [-55 7];
pl.ylim = [0 0.8];
for i = 1:8
    subplot(2,4,i);
    title(pl.titles{i})
%     xlim(pl.xlim)
    ylim(pl.ylim)
    xlabel(pl.xlabel)
    hold on; plot([0 0], [0 1],'r-.','LineWidth',1.4);%xlim([-150 60]);
end


% summarize per day relative losses and plot
medians(:,1) = nanmedian(ML.TR.loss(ind11,:))';
medians(:,2) = nanmedian(ML.TR.loss(ind12,:))';
medians(:,3) = nanmedian(ML.TR.loss(ind13,:))';
medians(:,4) = nanmedian(ML.TR.loss(ind14,:))';
medians(:,5) = nanmedian(ML.TR.loss(ind21,:))';
medians(:,6) = nanmedian(ML.TR.loss(ind22,:))';
medians(:,7) = nanmedian(ML.TR.loss(ind23,:))';
medians(:,8) = nanmedian(ML.TR.loss(ind24,:))';
qrtls25(:,1) = quantile(ML.TR.loss(ind11,:),0.25)';
qrtls25(:,2) = quantile(ML.TR.loss(ind12,:),0.25)';
qrtls25(:,3) = quantile(ML.TR.loss(ind13,:),0.25)';
qrtls25(:,4) = quantile(ML.TR.loss(ind14,:),0.25)';
qrtls25(:,5) = quantile(ML.TR.loss(ind21,:),0.25)';
qrtls25(:,6) = quantile(ML.TR.loss(ind22,:),0.25)';
qrtls25(:,7) = quantile(ML.TR.loss(ind23,:),0.25)';
qrtls25(:,8) = quantile(ML.TR.loss(ind24,:),0.25)';
qrtls75(:,1) = quantile(ML.TR.loss(ind11,:),0.75)';
qrtls75(:,2) = quantile(ML.TR.loss(ind12,:),0.75)';
qrtls75(:,3) = quantile(ML.TR.loss(ind13,:),0.75)';
qrtls75(:,4) = quantile(ML.TR.loss(ind14,:),0.75)';
qrtls75(:,5) = quantile(ML.TR.loss(ind21,:),0.75)';
qrtls75(:,6) = quantile(ML.TR.loss(ind22,:),0.75)';
qrtls75(:,7) = quantile(ML.TR.loss(ind23,:),0.75)';
qrtls75(:,8) = quantile(ML.TR.loss(ind24,:),0.75)';

figure('Units','Normalized','OuterPosition',[0.05 0.05 0.8 0.9]); box on;
pl.cols = [0.2 0.8 1; 0.2 0.6 1;0.2 0 1;0 0 0.6;0.2 0.8 1; 0.2 0.6 1;0.2 0 1;0 0 0.6];
for i = 1:4
    subplot(4,2,2*i-1); box on; hold on
    h = area(-5:30,qrtls75(:,2*i-1),-60);
    h.FaceColor = pl.cols(i,:);
    h.EdgeColor = pl.cols(i,:);
    h.FaceAlpha = 0.5;
    h = area(-5:30,qrtls25(:,2*i-1),-60);
    h.FaceColor = [1 1 1];
    h.EdgeColor = pl.cols(i,:);
    plot(-5:30,medians(:,i),'LineWidth',3,'Color',pl.cols(i,:));
    plot(-5:30, zeros(length(-5:30),1),'k-.','LineWidth',1)
    subplot(4,2,2*i); box on; hold on
    h = area(-5:30,qrtls75(:,2*i),-60);
    h.FaceColor = pl.cols(i,:);
    h.EdgeColor = pl.cols(i,:);
    h.FaceAlpha = 0.5;
    h = area(-5:30,qrtls25(:,2*i),-60);
    h.FaceColor = [1 1 1];
    h.EdgeColor = pl.cols(i,:);
    plot(-5:30,medians(:,2*i),'LineWidth',3,'Color',pl.cols(i,:));
    plot(-5:30, zeros(length(-5:30),1),'k-.','LineWidth',1)
end

pl.titles = {'Lac = 1, LS = 1','Lac > 1, LS = 1','Lac = 1, LS = 2','Lac > 1, LS = 2',...
             'Lac = 1, LS = 3','Lac > 1, LS = 3','Lac = 1, LS = 4','Lac > 1, LS = 4'};
pl.xlabel = {'Days from detection'};      
pl.xlim = [-4.9 29.9];
pl.ylim = [-60 20];
for i = 1:8
    subplot(4,2,i);
    title(pl.titles{i})
%     xlim(pl.xlim)
    ylim(pl.ylim)
    xlabel(pl.xlabel)
    plot([0 0],[-60 20],'r-','LineWidth',1)
end

% plot groups LS/P
figure('Units','Normalized','OuterPosition',[0.1 0.1 0.5 0.4]);
binedges = [0 64 138 216 305]; 
subplot(1,2,1);box on; hold on; title('Proportion of cases per LS')
histogram(ML.TR.DIM(ML.TR.Lac > 1),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0 0 0.6],...
                        'FaceColor',[0 0 0.6])
histogram(ML.TR.DIM(ML.TR.Lac == 1),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0.8 1],...
                        'FaceColor',[0.2 0.8 1])
                    legend({'Lac > 1','Lac = 1'}); xlabel('DIM [days]')
subplot(1,2,2);box on; hold on; title('Raw number of cases per LS')
histogram(ML.TR.DIM(ML.TR.Lac > 1),binedges,...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0 0 0.6],...
                        'FaceColor',[0 0 0.6])
histogram(ML.TR.DIM(ML.TR.Lac == 1),binedges,...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0.8 1],...
                        'FaceColor',[0.2 0.8 1])
                    legend({'Lac > 1','Lac = 1'}); xlabel('DIM [days]')

disp(['No. cases < DIM10 = ' num2str(sum(ML.TR.DIM < 10)) ', ' num2str(sum(ML.TR.DIM < 10)./length(ML.TR.DIM < 10)*100) '%'])



% no farm analysis add, too different numbers of cases...


disp(['Average milk loss = ' num2str(nanmean(sum(ML.TR{:,[12 15]},2))) ' +/- ' num2str(nanstd(sum(ML.TR{:,[12 15]},2)))])
disp(['Average milk loss rel = ' num2str(nanmean(nanmean(ML.TR{:,[13 16]},2))) ' +/- ' num2str(nanstd(nanmean(ML.TR{:,[13 16]},2)))])

figure; subplot(1,2,1)
binedges = -800:100:0;
histogram(sum(ML.TR{:,[12 15]},2),binedges,'LineWidth',2,'EdgeColor',[0.2 0.8 1],...
                        'FaceColor',[0.2 0.8 1])
title('Absolute milk loss in -5 to +21 days from treatment')
xlabel('kg milk loss'); ylabel('No. of cases')
subplot(1,2,2)
binedges = -30:2:0;
histogram(mean(ML.TR{:,16},2),binedges,'LineWidth',2,'EdgeColor',[0.2 0.8 1],...
                        'FaceColor',[0.2 0.8 1])
title('Relative daily milk loss in 21 days after treatment')
xlabel('% daily milk loss'); ylabel('No. of cases')                    


        
%% TREATMENT REGISTERS - QL DATA
                    
% >> 1 << Number of cases matched with 'complete' lactation data
farms = fieldnames(CLMAS);  % farms with CM data
ML.TRQ = table(zeros(0,1),'VariableNames',{'B_ID'}); % prepare case table
ML.TRQ.FarmName = repmat({''},0);
ML.TRQ = movevars(ML.TRQ, 'FarmName','Before','B_ID');
mod.Wood = @(p,t) p(1).*t.^p(2).*exp(-p(3).*t);     % wood model - prepare
T = 1;  % counter to fill table
for i = 1:length(farms)
    
    % first select cases not too close to eachother (at least 10 days
    % apart)
    data = CLMAS.(farms{i})(~isnan(CLMAS.(farms{i}).B_ID),:);
    data.diff(1,1) = 10;
    data.diff(2:end,1) = diff(datenum(data.Date));
    [~,indx] = unique(data.B_ID); % find first occurence
    data.first(indx) = 1;
    data(data.first == 0 & data.diff < 10,:) = []; % delete
    % re-add to CLMAS
    data = removevars(data,{'diff','first'});
    CLMAS.(farms{i}) = data; 
    clear data

    
    ind = find(~isnan(CLMAS.(farms{i}).B_ID));
        
    L = length(ind); % 
    ML.TRQ.FarmName(T:T+L-1,1) = repmat(farms(i),L,1);
    ML.TRQ.B_ID(T:T+L-1) = CLMAS.(farms{i}).B_ID(ind);
    ML.TRQ.Date(T:T+L-1) = CLMAS.(farms{i}).Date(ind);
    
    for ii = 1:length(ind)
        % find corresponding lactation and DIM
        try
            % find the corresponding data in a lactation within 305 dim
            idx = find(moddata.(farms{i}).B_ID == ML.TRQ.B_ID(T+ii-1) & ...
                       floor(datenum(moddata.(farms{i}).Date)) == floor(datenum(ML.TRQ.Date(T+ii-1))));     
        catch
            % fill in empty
            idx = [];
        end
        % fill in
        if ~isempty(idx)

            ML.TRQ.Lac(T+ii-1) = moddata.(farms{i}).Lac(idx(1));
            ML.TRQ.DIM(T+ii-1) = moddata.(farms{i}).DIM(idx(1));
            ML.TRQ.Q(T+ii-1) = CLMAS.(farms{i}).ms(ind(ii));
            ML.TRQ.Month(T+ii-1) = month(moddata.(farms{i}).Date(idx(1)));

            % data to calculate losses
            data = moddata.(farms{i})(moddata.(farms{i}).B_ID ==  ML.TRQ.B_ID(T+ii-1) & ...
                                         moddata.(farms{i}).Lac ==  ML.TRQ.Lac(T+ii-1),[2 6 7 8 10 11 13:25 ]);
                                     
            % summary expected lactation 305
            ML.TRQ.LacLength(T+ii-1) = max(data.DIM);     % lactation length
            ML.TRQ.PeakMYexp1(T+ii-1) = max(data.ModIT1); % expected peak yield
            ML.TRQ.PeakMYexp2(T+ii-1) = max(data.ModIT2); % expected peak yield
            ML.TRQ.PeakMYexp3(T+ii-1) = max(data.ModIT3); % expected peak yield
            ML.TRQ.PeakMYexp4(T+ii-1) = max(data.ModIT4); % expected peak yield
            ML.TRQ.TotMYexp1(T+ii-1) = nansum(data.ModIT1);  % expected total yield
            ML.TRQ.TotMYexp2(T+ii-1) = nansum(data.ModIT2);  % expected total yield
            ML.TRQ.TotMYexp3(T+ii-1) = nansum(data.ModIT3);  % expected total yield
            ML.TRQ.TotMYexp4(T+ii-1) = nansum(data.ModIT4);  % expected total yield      
            
            % prepare milk loss calculations
            data.RES1 = data.MYLF-data.ModIT1;  % residual milk yield
            data.RES2 = data.MYRF-data.ModIT2;  % residual milk yield
            data.RES3 = data.MYLH-data.ModIT3;  % residual milk yield
            data.RES4 = data.MYRH-data.ModIT4;  % residual milk yield
            
            data.RELRES1 = 100*(data.RES1./data.ModIT1); data.RELRES1(data.RELRES1>100) = -0.1;
            data.RELRES2 = 100*(data.RES2./data.ModIT2); data.RELRES2(data.RELRES2>100) = -0.1;
            data.RELRES3 = 100*(data.RES3./data.ModIT3); data.RELRES3(data.RELRES3>100) = -0.1;
            data.RELRES4 = 100*(data.RES4./data.ModIT4); data.RELRES4(data.RELRES4>100) = -0.1;
            
            
            % expected production in period around detection

            
            % Q1
            ML.TRQ.ExpMY1(T+ii-1)  = nansum(data.ModIT1(floor(data.DIM) >= floor(ML.TRQ.DIM(T+ii-1)-5) & ...
                                               floor(data.DIM) <= floor(ML.TRQ.DIM(T+ii-1)+30)));
            ML.TRQ.ProdMY1(T+ii-1) = nansum(data.MYLF(floor(data.DIM) >= floor(ML.TRQ.DIM(T+ii-1)-5) & ...
                                               floor(data.DIM) <= floor(ML.TRQ.DIM(T+ii-1)+30)));
            ML.TRQ.TotRel1(T+ii-1) = (ML.TRQ.ProdMY1(T+ii-1) - ML.TRQ.ExpMY1(T+ii-1))/ML.TRQ.ExpMY1(T+ii-1)*100;
            % Q2
            ML.TRQ.ExpMY2(T+ii-1)  = nansum(data.ModIT2(floor(data.DIM) >= floor(ML.TRQ.DIM(T+ii-1)-5) & ...
                                               floor(data.DIM) <= floor(ML.TRQ.DIM(T+ii-1)+30)));
            ML.TRQ.ProdMY2(T+ii-1) = nansum(data.MYRF(floor(data.DIM) >= floor(ML.TRQ.DIM(T+ii-1)-5) & ...
                                               floor(data.DIM) <= floor(ML.TRQ.DIM(T+ii-1)+30)));
            ML.TRQ.TotRel2(T+ii-1) = (ML.TRQ.ProdMY2(T+ii-1) - ML.TRQ.ExpMY2(T+ii-1))/ML.TRQ.ExpMY2(T+ii-1)*100;
            % Q3
            ML.TRQ.ExpMY3(T+ii-1)  = nansum(data.ModIT3(floor(data.DIM) >= floor(ML.TRQ.DIM(T+ii-1)-5) & ...
                                               floor(data.DIM) <= floor(ML.TRQ.DIM(T+ii-1)+30)));
            ML.TRQ.ProdMY3(T+ii-1) = nansum(data.MYLH(floor(data.DIM) >= floor(ML.TRQ.DIM(T+ii-1)-5) & ...
                                               floor(data.DIM) <= floor(ML.TRQ.DIM(T+ii-1)+30)));
            ML.TRQ.TotRel3(T+ii-1) = (ML.TRQ.ProdMY3(T+ii-1) - ML.TRQ.ExpMY3(T+ii-1))/ML.TRQ.ExpMY3(T+ii-1)*100;
            % Q4
            ML.TRQ.ExpMY4(T+ii-1)  = nansum(data.ModIT4(floor(data.DIM) >= floor(ML.TRQ.DIM(T+ii-1)-5) & ...
                                               floor(data.DIM) <= floor(ML.TRQ.DIM(T+ii-1)+30)));
            ML.TRQ.ProdMY4(T+ii-1) = nansum(data.MYRH(floor(data.DIM) >= floor(ML.TRQ.DIM(T+ii-1)-5) & ...
                                               floor(data.DIM) <= floor(ML.TRQ.DIM(T+ii-1)+30)));
            ML.TRQ.TotRel4(T+ii-1) = (ML.TRQ.ProdMY4(T+ii-1) - ML.TRQ.ExpMY4(T+ii-1))/ML.TRQ.ExpMY4(T+ii-1)*100;
            
            % IQR EC
            indexcase = find(floor(data.DIM) >= floor(ML.TRQ.DIM(T+ii-1)) & ...
                             floor(data.DIM) <= floor(ML.TRQ.DIM(T+ii-1)+10));
            ML.TRQ.IQR1(T+ii-1) = nanmean(data.ECLF(indexcase))./nanmean(nanmean([data.ECRF(indexcase),data.ECLH(indexcase),data.ECRH(indexcase)]));
            ML.TRQ.IQR2(T+ii-1) = nanmean(data.ECRF(indexcase))./nanmean(nanmean([data.ECLF(indexcase),data.ECLH(indexcase),data.ECRH(indexcase)]));
            ML.TRQ.IQR3(T+ii-1) = nanmean(data.ECLH(indexcase))./nanmean(nanmean([data.ECLF(indexcase),data.ECRF(indexcase),data.ECRH(indexcase)]));
            ML.TRQ.IQR4(T+ii-1) = nanmean(data.ECRH(indexcase))./nanmean(nanmean([data.ECLF(indexcase),data.ECRF(indexcase),data.ECLH(indexcase)]));
            
            
            % Q1 calculate losses -5 to 30 days from detection absolute
            for iii = -5:30
                loss = data.RES1(floor(data.DIM) == floor(ML.TRQ.DIM(T+ii-1)+iii));
                if ~isempty(loss) && nansum(loss) ~= 0 && ~isnan(nansum(loss))
                    ML.TRQ.loss1(T+ii-1,iii+6) = nansum(loss);
                else
                    ML.TRQ.loss1(T+ii-1,iii+6) = NaN;
                end
            end
            % Q2 calculate losses -5 to 30 days from detection absolute
            for iii = -5:30
                loss = data.RES2(floor(data.DIM) == floor(ML.TRQ.DIM(T+ii-1)+iii));
                if ~isempty(loss) && nansum(loss) ~= 0 && ~isnan(nansum(loss))
                    ML.TRQ.loss2(T+ii-1,iii+6) = nansum(loss);
                else
                    ML.TRQ.loss2(T+ii-1,iii+6) = NaN;
                end
            end
            % Q3 calculate losses -5 to 30 days from detection absolute
            for iii = -5:30
                loss = data.RES3(floor(data.DIM) == floor(ML.TRQ.DIM(T+ii-1)+iii));
                if ~isempty(loss) && nansum(loss) ~= 0 && ~isnan(nansum(loss))
                    ML.TRQ.loss3(T+ii-1,iii+6) = nansum(loss);
                else
                    ML.TRQ.loss3(T+ii-1,iii+6) = NaN;
                end
            end
            % Q4 calculate losses -5 to 30 days from detection absolute
            for iii = -5:30
                loss = data.RES4(floor(data.DIM) == floor(ML.TRQ.DIM(T+ii-1)+iii));
                if ~isempty(loss) && nansum(loss) ~= 0 && ~isnan(nansum(loss))
                    ML.TRQ.loss4(T+ii-1,iii+6) = nansum(loss);
                else
                    ML.TRQ.loss4(T+ii-1,iii+6) = NaN;
                end
            end
                                                                                       
            
            
            
            % Q1 calculate losses -5 to 30 days from detection absolute
            for iii = -5:30
                loss = data.RELRES1(floor(data.DIM) == floor(ML.TRQ.DIM(T+ii-1)+iii));
                if ~isempty(loss) && nansum(loss) ~= 0 && ~isnan(nansum(loss))
                    ML.TRQ.lossREL1(T+ii-1,iii+6) = nanmean(loss);
                else
                    ML.TRQ.lossREL1(T+ii-1,iii+6) = NaN;
                end
            end
            
            % Q2 calculate losses -5 to 30 days from detection absolute
            for iii = -5:30
                loss = data.RELRES2(floor(data.DIM) == floor(ML.TRQ.DIM(T+ii-1)+iii));
                if ~isempty(loss)
                    ML.TRQ.lossREL2(T+ii-1,iii+6) = nanmean(loss);
                else
                    ML.TRQ.lossREL2(T+ii-1,iii+6) = NaN;
                end
            end
            % Q3 calculate losses -5 to 30 days from detection absolute
            for iii = -5:30
                loss = data.RELRES3(floor(data.DIM) == floor(ML.TRQ.DIM(T+ii-1)+iii));
                if ~isempty(loss)
                    ML.TRQ.lossREL3(T+ii-1,iii+6) = nanmean(loss);
                else
                    ML.TRQ.lossREL3(T+ii-1,iii+6) = NaN;
                end
            end
            % Q4 calculate losses -5 to 30 days from detection absolute
            for iii = -5:30
                loss = data.RELRES4(floor(data.DIM) == floor(ML.TRQ.DIM(T+ii-1)+iii));
                if ~isempty(loss)
                    ML.TRQ.lossREL4(T+ii-1,iii+6) = nanmean(loss);
                else
                    ML.TRQ.lossREL4(T+ii-1,iii+6) = NaN;
                end
            end
            
            % identify quarters with highest losses
            SUSQ = [ML.TRQ{T+ii-1,[17 20 23 26]} ; ML.TRQ{T+ii-1,[18 21 24 27]};...
                    ML.TRQ{T+ii-1,[19 22 25 28]} ;  ML.TRQ{T+ii-1,29:32}]; 
            
            M = min(SUSQ(3,:)); %  min % loss
            s = find(SUSQ(3,:)-M < 1 & SUSQ(4,:) >= 1); % Q 
            if length(s) == 1      
                ML.TRQ.susQ(T+ii-1) = s;
            else
                ML.TRQ.susQ(T+ii-1) = NaN;
            end
            
        else
            ML.TRQ.Lac(T+ii-1) = NaN;
            ML.TRQ.DIM(T+ii-1) = NaN;
        end
    end
    T = T+L;
end
clear sub T L i ii idx ind data DIM ITW ans p test h

disp(['No. of unmatched cases = ' num2str(length(find(isnan(ML.TRQ.Lac)))) ...
       ', = ' num2str(length(find(isnan(ML.TRQ.Lac)))./height(ML.TRQ)*100) '%'])
disp(['No. of matched cases = ' num2str(length(find(~isnan(ML.TRQ.Lac)))) ...
       ', = ' num2str(length(find(~isnan(ML.TRQ.Lac)))./height(ML.TRQ)*100) '%'])

% delete unmatched cases
ML.TRQ = ML.TRQ(~isnan(ML.TRQ.DIM),:);

% delete cases with no Q info
ind = find(nansum(ML.TRQ{:,9:12},2)==0);  % 158 cases zero q yield
ML.TRQ(ind,:) = [];

% 
disp(['No. of cases with no clear Q diff = ' num2str(length(find(isnan(ML.TRQ.susQ)))) ...
       ', = ' num2str(length(find(isnan(ML.TRQ.susQ)))./height(ML.TRQ)*100) '%'])
% plot histogram quarters
figure('Units','Normalized','OuterPosition',[0.2 0.2 0.25 0.45])
h = histogram(ML.TRQ.susQ);
x = h.BinEdges;
y = h.Values;
text(x(1:end-1)+0.5,y,num2str(y'),'vert','bottom','horiz','center')
h.FaceColor = [0 0.6 0.6]; h.FaceAlpha =1; h.LineWidth = 1;
ylim([0 475])
xticks([1 2 3 4]); xticklabels({'LF','RF','LH','RH'})
set(gca,'FontName','Times New Roman','FontSize',12)
xlabel('No. of cases');
title('Case detection by QMY / EC')

% correspondance
test = find(ML.TRQ.Q<5 & ML.TRQ.Q > 0 & ~isnan(ML.TRQ.Q));
disp(['% quarter not corresp = ' ...
       num2str(length(find(ML.TRQ.susQ(test) - ML.TRQ.Q(test) ~= 0))./...
       length(test)*100) '%']);


% exploratory
figure; subplot(1,3,1); histogram(ML.TRQ.Lac,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.7 0 0],...
                        'FaceColor',[0.7 0 0])
        title('Parity'); xlabel('Parity number')
        subplot(1,3,2); histogram(ML.TRQ.DIM,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0 0.7 0],...
                        'FaceColor',[0 0.7 0])
        title('Lactation stage'); xlabel('DIM (days)')
        subplot(1,3,3); histogram(ML.TRQ.Month,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.4 0 0.8],...
                        'FaceColor',[0.4 0 0.8])
        title('Month'); xlabel('Month')
disp(['No. of cases in M6/7/8 cases = ' num2str(length(find(ismember(ML.TRQ.Month,[6 7 9])))) ...
       ', = ' num2str(length(find(ismember(ML.TRQ.Month,[6 7 9])))./height(ML.TRQ)*100) '%'])
disp(['No of cases in first 10 days = ' num2str(length(find(ML.TRQ.DIM<10))) ...
       ', = ' num2str(length(find(ML.TRQ.DIM<10))./height(ML.TRQ)*100) '%'])
disp(['No of cases in first 20 days = ' num2str(length(find(ML.TRQ.DIM<21))) ...
       ', = ' num2str(length(find(ML.TRQ.DIM<21))./height(ML.TRQ)*100) '%'])
        
% plot groups
figure('Units','Normalized','OuterPosition',[0.1 0.1 0.5 0.4]);
binedges = [0 64 138 216 305]; 
subplot(1,2,1);box on; hold on; title('Proportion of cases per LS')
histogram(ML.TRQ.DIM(ML.TRQ.Lac > 1),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0 0 0.6],...
                        'FaceColor',[0 0 0.6])
histogram(ML.TRQ.DIM(ML.TRQ.Lac == 1),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0.8 1],...
                        'FaceColor',[0.2 0.8 1])
                    legend({'Lac > 1','Lac = 1'}); xlabel('DIM [days]')
subplot(1,2,2);box on; hold on; title('Raw number of cases per LS')
histogram(ML.TRQ.DIM(ML.TRQ.Lac > 1),binedges,...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0 0 0.6],...
                        'FaceColor',[0 0 0.6])
histogram(ML.TRQ.DIM(ML.TRQ.Lac == 1),binedges,...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0.8 1],...
                        'FaceColor',[0.2 0.8 1])
                    legend({'Lac > 1','Lac = 1'}); xlabel('DIM [days]')

% visualize milk losses
ind11 = find(ML.TRQ.Lac == 1 & ML.TRQ.DIM <= 63);
ind12 = find(ML.TRQ.Lac == 1 & ML.TRQ.DIM > 63 & ML.TRQ.DIM <= 138);
ind13 = find(ML.TRQ.Lac == 1 & ML.TRQ.DIM > 138 & ML.TRQ.DIM <= 216);
ind14 = find(ML.TRQ.Lac == 1 & ML.TRQ.DIM > 216);
ind21 = find(ML.TRQ.Lac > 1 & ML.TRQ.DIM <= 63);
ind22 = find(ML.TRQ.Lac > 1 & ML.TRQ.DIM > 63 & ML.TRQ.DIM <= 138);
ind23 = find(ML.TRQ.Lac > 1 & ML.TRQ.DIM > 138 & ML.TRQ.DIM <= 216);
ind24 = find(ML.TRQ.Lac > 1 & ML.TRQ.DIM > 216);


% distribution quarter 
figure;
ind = find(ML.TRQ.Q >0 & ML.TRQ.Q < 5);
scatter(1:length(ind),ML.TRQ.Q(ind)-ML.TRQ.susQ(ind));
disp(['No. correct Q = ' num2str(length(find(ML.TRQ.Q == ML.TRQ.susQ))) ' out of ' num2str(length(ind))]);

figure('Units','Normalized','OuterPosition',[0.1 0.1 0.85 0.8]);
binedges = 1:1:5;
        subplot(2,4,1); histogram(ML.TRQ.susQ(ind11),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0.8 1],...
                        'FaceColor',[0.2 0.8 1])
        subplot(2,4,2); histogram(ML.TRQ.susQ(ind12),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0.6 1],...
                        'FaceColor',[0.2 0.6 1])
        subplot(2,4,3); histogram(ML.TRQ.susQ(ind13),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0 1],...
                        'FaceColor',[0.2 0 1])                
        subplot(2,4,4); histogram(ML.TRQ.susQ(ind14),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0 0 0.6],...
                        'FaceColor',[0 0 0.6])
        subplot(2,4,5); histogram(ML.TRQ.susQ(ind21),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0.8 1],...
                        'FaceColor',[0.2 0.8 1])
        subplot(2,4,6); histogram(ML.TRQ.susQ(ind22),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0.6 1],...
                        'FaceColor',[0.2 0.6 1])
        subplot(2,4,7); histogram(ML.TRQ.susQ(ind23),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0 1],...
                        'FaceColor',[0.2 0 1])                
        subplot(2,4,8); histogram(ML.TRQ.susQ(ind24),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0 0 0.6],...
                        'FaceColor',[0 0 0.6])
    
pl.titles = {'Lac = 1, LS = 1','Lac = 1, LS = 2','Lac = 1, LS = 3','Lac = 1, LS = 4',...
             'Lac > 1, LS = 1','Lac > 1, LS = 2','Lac > 1, LS = 3','Lac > 1, LS = 4'};
pl.xlabel = {'Suspected quarter'};      
% pl.ylim = [0 0.8];
for i = 1:8
    subplot(2,4,i);
    title(pl.titles{i})
    xticks([1.5 2.5 3.5 4.5])
    xticklabels({'LF','RF','LH','RH'})
    xlabel(pl.xlabel)
%     hold on; plot([0 0], [0 1],'r-.','LineWidth',1.4);%xlim([-150 60]);
end


% Examples
i = 1500; % 985/1598 2000 199 200 585 41
    pl.farm = ML.TRQ.FarmName{i}; pl.B_ID = ML.TRQ.B_ID(i); pl.Lac = ML.TRQ.Lac(i);
    ind = find(moddata.(pl.farm).B_ID == pl.B_ID & moddata.(pl.farm).Lac == pl.Lac);
    pl.DIM = moddata.(pl.farm).DIM(ind);
    pl.QMY = moddata.(pl.farm){ind,13:16}; pl.MOD = moddata.(pl.farm){ind,22:25};
    pl.QEC = moddata.(pl.farm){ind,17:20};
    if mean(pl.QEC > 20); pl.QEC = pl.QEC./10; end
    pl.susQ = ML.TRQ.susQ(i);
    ind = find(moddataDAY.(pl.farm).B_ID == pl.B_ID & moddataDAY.(pl.farm).Lac == pl.Lac);
    pl.DIMD = moddataDAY.(pl.farm).DIM(ind);
    pl.DMY = moddataDAY.(pl.farm).MYs(ind); pl.MODD = moddataDAY.(pl.farm).ModIT(ind);
    

figure('Units','Normalized','OuterPosition',[0.1 0.1 0.65 0.8]); 
subplot(2,1,1); hold on; box on; xlabel('DIM (days)'); ylabel('Daily milk yield (kg)')
    xlim([min(pl.DIM) max(pl.DIM)])
    plot(pl.DIMD,pl.DMY,'o-','Color',[0 0.5 0.8],'LineWidth',1,'MarkerSize',3)
    plot(pl.DIMD,pl.MODD,'-','Color',[0 0.6 0.4],'LineWidth',2.5)
    plot([ML.TRQ.DIM(i) ML.TRQ.DIM(i)],[0 60],'r:','LineWidth',1.5)
    

subplot(2,1,2); hold on; box on; xlabel('DIM (days)'); ylabel('Quarter milk yield (kg)')
    xlim([min(pl.DIM) max(pl.DIM)])
    plot(pl.DIM,pl.QEC(:,pl.susQ),'s-','Color',[0.8 0.4 0.6],'LineWidth',0.8,'MarkerSize',3,'MarkerFaceColor',[0.5 0 0.5])
    plot(pl.DIM,pl.QMY(:,pl.susQ),'o-','Color',[0 0.5 0.8],'LineWidth',1,'MarkerSize',3)
    plot(pl.DIM,pl.MOD(:,pl.susQ),'-','Color',[0 0.6 0.4],'LineWidth',0.8)
    plot([ML.TRQ.DIM(i) ML.TRQ.DIM(i)],[0 15],'r:','LineWidth',1.5)



% plot milk losses in suspected vs in all quarters
nifq = []; ifq = [];
ind = find(ML.TRQ.susQ == 1);
nifq(ind,1:3) = [1*ones(length(ind),1) ML.TRQ.DIM(ind) ML.TRQ.Lac(ind)];
for j = 1:36
    ifq(ind,j) = ML.TRQ.lossREL1(ind,j);
    nifq(ind,j+3) = nanmean([ML.TRQ.lossREL2(ind,j),ML.TRQ.lossREL3(ind,j),ML.TRQ.lossREL4(ind,j)],2);
end
ind = find(ML.TRQ.susQ == 2);
nifq(ind,1:3) = [2*ones(length(ind),1) ML.TRQ.DIM(ind) ML.TRQ.Lac(ind)];
for j = 1:36
    ifq(ind,j) = ML.TRQ.lossREL2(ind,j);
    nifq(ind,j+3) = nanmean([ML.TRQ.lossREL1(ind,j),ML.TRQ.lossREL3(ind,j),ML.TRQ.lossREL4(ind,j)],2);
end
ind = find(ML.TRQ.susQ == 3);
nifq(ind,1:3) = [3*ones(length(ind),1) ML.TRQ.DIM(ind) ML.TRQ.Lac(ind)];
for j = 1:36
    ifq(ind,j) = ML.TRQ.lossREL3(ind,j);
    nifq(ind,j+3) = nanmean([ML.TRQ.lossREL1(ind,j),ML.TRQ.lossREL2(ind,j),ML.TRQ.lossREL4(ind,j)],2);
end
ind = find(ML.TRQ.susQ == 4);
nifq(ind,1:3) = [4*ones(length(ind),1) ML.TRQ.DIM(ind) ML.TRQ.Lac(ind)];
for j = 1:36
    ifq(ind,j) = ML.TRQ.lossREL4(ind,j);
    nifq(ind,j+3) = nanmean([ML.TRQ.lossREL1(ind,j),ML.TRQ.lossREL2(ind,j),ML.TRQ.lossREL3(ind,j)],2);
end
ifq(sum(ifq,2)==0,:) = [];
nifq(sum(nifq(:,4:end),2)==0,:) = [];

% plot infected vs non-infected hind vs front
pl.cols = [0.2 0.8 1; 0.2 0.6 1;0.2 0 1;0 0 0.6;0.2 0.8 1; 0.2 0.6 1;0.2 0 1;0 0 0.6];
figure('Units','Normalized','OuterPosition',[0.05 0.05 0.8 0.9]); box on;
subplot(2,2,1); hold on; box on; 
ind = find(nifq(:,1) == 1 | nifq(:,1) == 2);
h = area(-5:30,quantile(ifq(ind,:),0.75),-85);
    h.FaceColor = pl.cols(1,:);
    h.EdgeColor = pl.cols(1,:);
    h.FaceAlpha = 0.5;
h = area(-5:30,quantile(ifq(ind,:),0.25),-85);
    h.FaceColor = [1 1 1];
    h.EdgeColor = pl.cols(1,:);
plot(-5:30,nanmedian(ifq(ind,:)),'LineWidth',2.5,'Color',[0.2 0.8 1])
plot(-5:30,zeros(length(-5:30),1),'k-.','LineWidth',1)
ylim([-85 10]); title('Suspected quarter = front; Relative loss in suspected quarter')
ylabel('Relative milk loss (%)'); xlabel('Days from detection')
for i = 10:-20:-80; plot([-5 30],[i i],'LineWidth',0.01,'Color',[0.95 0.95 0.95]); end
for i = -5:5:30; plot([i i],[-90 10],'LineWidth',0.01,'Color',[0.95 0.95 0.95]); end
subplot(2,2,2); hold on; box on;
h = area(-5:30,quantile(nifq(ind,4:end),0.75),-85);
    h.FaceColor = pl.cols(1,:);
    h.EdgeColor = pl.cols(1,:);
    h.FaceAlpha = 0.5;
h = area(-5:30,quantile(nifq(ind,4:end),0.25),-85);
    h.FaceColor = [1 1 1];
    h.EdgeColor = pl.cols(1,:);
plot(-5:30,nanmedian(nifq(ind,4:end)),'LineWidth',2.5,'Color',[0.2 0.8 1])
plot(-5:30,zeros(length(-5:30),1),'k-.','LineWidth',1)
ylim([-85 10]);title('Suspected quarter = front; Average relative loss in non-suspected quarters')
ylabel('Relative milk loss (%)'); xlabel('Days from detection')
for i = 10:-20:-80; plot([-5 30],[i i],'LineWidth',0.01,'Color',[0.95 0.95 0.95]); end
for i = -5:5:30; plot([i i],[-90 10],'LineWidth',0.01,'Color',[0.95 0.95 0.95]); end
subplot(2,2,3); hold on; box on; 
ind = find(nifq(:,1) == 3 | nifq(:,1) == 4);
h = area(-5:30,quantile(ifq(ind,:),0.75),-85);
    h.FaceColor = pl.cols(1,:);
    h.EdgeColor = pl.cols(1,:);
    h.FaceAlpha = 0.5;
h = area(-5:30,quantile(ifq(ind,:),0.25),-85);
    h.FaceColor = [1 1 1];
    h.EdgeColor = pl.cols(1,:);
plot(-5:30,nanmedian(ifq(ind,:)),'LineWidth',2.5,'Color',[0.2 0.8 1])
plot(-5:30,zeros(length(-5:30),1),'k-.','LineWidth',1)
ylim([-85 10]); title('Suspected quarter = hind; Relative loss in suspected quarter')
ylabel('Relative milk loss (%)'); xlabel('Days from detection')
for i = 10:-20:-80; plot([-5 30],[i i],'LineWidth',0.01,'Color',[0.95 0.95 0.95]); end
for i = -5:5:30; plot([i i],[-90 10],'LineWidth',0.01,'Color',[0.95 0.95 0.95]); end
subplot(2,2,4); hold on; box on;
h = area(-5:30,quantile(nifq(ind,4:end),0.75),-85);
    h.FaceColor = pl.cols(1,:);
    h.EdgeColor = pl.cols(1,:);
    h.FaceAlpha = 0.5;
h = area(-5:30,quantile(nifq(ind,4:end),0.25),-85);
    h.FaceColor = [1 1 1];
    h.EdgeColor = pl.cols(1,:);
plot(-5:30,nanmedian(nifq(ind,4:end)),'LineWidth',2.5,'Color',[0.2 0.8 1])
plot(-5:30,zeros(length(-5:30),1),'k-.','LineWidth',1)
ylim([-85 10]);title('Suspected quarter = hind; Average relative loss in non-suspected quarters')
ylabel('Relative milk loss (%)'); xlabel('Days from detection')
for i = 10:-20:-80; plot([-5 30],[i i],'LineWidth',0.01,'Color',[0.95 0.95 0.95]); end
for i = -5:5:30; plot([i i],[-90 10],'LineWidth',0.01,'Color',[0.95 0.95 0.95]); end



% at quarter level:RELATIVE losses in day -5:-1 and day 0:30
figure('Units','Normalized','OuterPosition',[0.1 0.1 0.85 0.8]);
binedges = -100:10:30;
ind = find(nifq(:,2) < 63 & nifq(:,3) == 1);
        subplot(2,4,1); histogram(nanmean(ifq(ind,1:5),2),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.8 0 0],...
                        'FaceColor',[0.8 0 0])
            hold on; histogram(nanmean(nifq(ind,4:9),2),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0.8 1],...
                        'FaceColor',[0.2 0.8 1])
ind = find(nifq(:,2) >= 64 & nifq(:,2) <= 138 & nifq(:,3) == 1);
        subplot(2,4,2); histogram(nanmean(ifq(ind,1:5),2),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.8 0 0],...
                        'FaceColor',[0.8 0 0])
            hold on; histogram(nanmean(nifq(ind,4:9),2),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0.8 1],...
                        'FaceColor',[0.2 0.8 1])    
ind = find(nifq(:,2) >= 139 & nifq(:,2) < 217 & nifq(:,3) == 1);
        subplot(2,4,3); histogram(nanmean(ifq(ind,1:5),2),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.8 0 0],...
                        'FaceColor',[0.8 0 0])
            hold on; histogram(nanmean(nifq(ind,4:9),2),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0.8 1],...
                        'FaceColor',[0.2 0.8 1]) 
ind = find(nifq(:,2) >= 217 & nifq(:,3) == 1);
        subplot(2,4,4); histogram(nanmean(ifq(ind,1:5),2),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.8 0 0],...
                        'FaceColor',[0.8 0 0])
            hold on; histogram(nanmean(nifq(ind,4:9),2),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0.8 1],...
                        'FaceColor',[0.2 0.8 1]) 
ind = find(nifq(:,2) < 63 & nifq(:,3) > 1);
        subplot(2,4,5); histogram(nanmean(ifq(ind,1:5),2),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.8 0 0],...
                        'FaceColor',[0.8 0 0])
            hold on; histogram(nanmean(nifq(ind,4:9),2),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0.8 1],...
                        'FaceColor',[0.2 0.8 1])
ind = find(nifq(:,2) >= 64 & nifq(:,2) <= 138 & nifq(:,3) > 1);
        subplot(2,4,6); histogram(nanmean(ifq(ind,1:5),2),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.8 0 0],...
                        'FaceColor',[0.8 0 0])
            hold on; histogram(nanmean(nifq(ind,4:9),2),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0.8 1],...
                        'FaceColor',[0.2 0.8 1])    
ind = find(nifq(:,2) >= 139 & nifq(:,2) < 217 & nifq(:,3) > 1);
        subplot(2,4,7); histogram(nanmean(ifq(ind,1:5),2),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.8 0 0],...
                        'FaceColor',[0.8 0 0])
            hold on; histogram(nanmean(nifq(ind,4:9),2),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0.8 1],...
                        'FaceColor',[0.2 0.8 1]) 
ind = find(nifq(:,2) >= 217 & nifq(:,3) > 1);
        subplot(2,4,8); histogram(nanmean(ifq(ind,1:5),2),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.8 0 0],...
                        'FaceColor',[0.8 0 0])
            hold on; histogram(nanmean(nifq(ind,4:9),2),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0.8 1],...
                        'FaceColor',[0.2 0.8 1]) 
  
pl.titles = {'Lac = 1, LS = 1','Lac = 1, LS = 2','Lac = 1, LS = 3','Lac = 1, LS = 4',...
             'Lac > 1, LS = 1','Lac > 1, LS = 2','Lac > 1, LS = 3','Lac > 1, LS = 4'};
pl.xlabel = {'Average relative milk loss 5 days before detection (%)'};      
% pl.xlim = [-800 100];
pl.ylim = [0 0.45];
for i = 1:8
    subplot(2,4,i);
    title(pl.titles{i})
%     xlim(pl.xlim)
    ylim(pl.ylim)
    xlabel(pl.xlabel)

    hold on; plot([0 0], [0 1],'r-.','LineWidth',1.4);%xlim([-150 60]);
    legend({'susinf q','susninf q'},'AutoUpdate','off','Location','northwest')
end
        
% at quarter level:RELATIVE losses in day 0:30
figure('Units','Normalized','OuterPosition',[0.1 0.1 0.85 0.8]);
binedges = -100:10:30;
ind = find(nifq(:,2) < 63 & nifq(:,3) == 1);
        subplot(2,4,1); histogram(nanmean(ifq(ind,6:36),2),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.8 0 0],...
                        'FaceColor',[0.8 0 0])
            hold on; histogram(nanmean(nifq(ind,10:end),2),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0.8 1],...
                        'FaceColor',[0.2 0.8 1])
ind = find(nifq(:,2) >= 64 & nifq(:,2) <= 138 & nifq(:,3) == 1);
        subplot(2,4,2); histogram(nanmean(ifq(ind,6:36),2),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.8 0 0],...
                        'FaceColor',[0.8 0 0])
            hold on; histogram(nanmean(nifq(ind,10:end),2),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0.8 1],...
                        'FaceColor',[0.2 0.8 1])    
ind = find(nifq(:,2) >= 139 & nifq(:,2) < 217 & nifq(:,3) == 1);
        subplot(2,4,3); histogram(nanmean(ifq(ind,6:36),2),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.8 0 0],...
                        'FaceColor',[0.8 0 0])
            hold on; histogram(nanmean(nifq(ind,10:end),2),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0.8 1],...
                        'FaceColor',[0.2 0.8 1]) 
ind = find(nifq(:,2) >= 217 & nifq(:,3) == 1);
        subplot(2,4,4); histogram(nanmean(ifq(ind,6:36),2),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.8 0 0],...
                        'FaceColor',[0.8 0 0])
            hold on; histogram(nanmean(nifq(ind,10:end),2),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0.8 1],...
                        'FaceColor',[0.2 0.8 1]) 
ind = find(nifq(:,2) < 63 & nifq(:,3) > 1);
        subplot(2,4,5); histogram(nanmean(ifq(ind,6:36),2),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.8 0 0],...
                        'FaceColor',[0.8 0 0])
            hold on; histogram(nanmean(nifq(ind,10:end),2),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0.8 1],...
                        'FaceColor',[0.2 0.8 1])
ind = find(nifq(:,2) >= 64 & nifq(:,2) <= 138 & nifq(:,3) > 1);
        subplot(2,4,6); histogram(nanmean(ifq(ind,6:36),2),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.8 0 0],...
                        'FaceColor',[0.8 0 0])
            hold on; histogram(nanmean(nifq(ind,10:end),2),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0.8 1],...
                        'FaceColor',[0.2 0.8 1])    
ind = find(nifq(:,2) >= 139 & nifq(:,2) < 217 & nifq(:,3) > 1);
        subplot(2,4,7); histogram(nanmean(ifq(ind,6:36),2),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.8 0 0],...
                        'FaceColor',[0.8 0 0])
            hold on; histogram(nanmean(nifq(ind,10:end),2),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0.8 1],...
                        'FaceColor',[0.2 0.8 1]) 
ind = find(nifq(:,2) >= 217 & nifq(:,3) > 1);
        subplot(2,4,8); histogram(nanmean(ifq(ind,6:36),2),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.8 0 0],...
                        'FaceColor',[0.8 0 0])
            hold on; histogram(nanmean(nifq(ind,10:end),2),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0.8 1],...
                        'FaceColor',[0.2 0.8 1]) 
  
pl.titles = {'Lac = 1, LS = 1','Lac = 1, LS = 2','Lac = 1, LS = 3','Lac = 1, LS = 4',...
             'Lac > 1, LS = 1','Lac > 1, LS = 2','Lac > 1, LS = 3','Lac > 1, LS = 4'};
pl.xlabel = {'Average relative milk loss 30 days after detection (%)'};      
% pl.xlim = [-800 100];
pl.ylim = [0 0.65];
for i = 1:8
    subplot(2,4,i);
    title(pl.titles{i})
%     xlim(pl.xlim)
    ylim(pl.ylim)
    xlabel(pl.xlabel)

    hold on; plot([0 0], [0 1],'r-.','LineWidth',1.4);%xlim([-150 60]);
    legend({'susinf q','susninf q'},'AutoUpdate','off','Location','northwest')
end

% difference in slopes smoothed
for i = 1:length(ifq(:,1))
    sminf = smooth(ifq(i,:));
    smnif = smooth(nifq(i,4:end));
    slopedif(i,:) = diff(sminf)./diff(smnif);
end
medslope = [];
medslope(:,1) = nanmedian(slopedif,1);
medslope(:,2) = quantile(slopedif,0.25);
medslope(:,3) = quantile(slopedif,0.75);

figure; hold on;box on; xlim([-5 29])
h = area(-5:29,medslope(:,3),-0.5);
    h.FaceAlpha = 0.3;
h = area(-5:29,medslope(:,2),-0.5);
    h.FaceColor = 'w';
plot(-5:29,medslope(:,1),'Color',[0.2 0.6 0.6],'LineWidth',3)
plot([0 0],[-0.5 3],'r-.','LineWidth',1.6)
plot(-5:29,ones(length(-5:29),1),'k--')
xlabel('Days from treatment')
title('Difference in slopes between inf and non-inf quarters')
ylabel('% per day / % per day')


% clear variables
clear abs binedges cowlac farms h i ind ind11 ind12 ind13 ind14 medians medslope qlosrel
clear ind21 ind22 ind23 ind24 indexcase indx infected j loss m M niet ninf iii ifq ans
clear nifq inf qlosREL qrtls25 qrtls75 Qs s settings1 slopedif sminf smnif SUSQ test x y

save('D10_ML.mat','ML')

%% treatment corresponding to a perturbation in daily milk yield
% percentage overlap between perturbation daily and TR

farms = fieldnames(CLMAS);  % farms with CM data
ML.TRP = table(zeros(0,1),'VariableNames',{'B_ID'}); % prepare case table
ML.TRP.FarmName = repmat({''},0);
ML.TRP = movevars(ML.TRP, 'FarmName','Before','B_ID');
mod.Wood = @(p,t) p(1).*t.^p(2).*exp(-p(3).*t);     % wood model - prepare
T = 1;  % counter to fill table
for i = 1:length(farms)
    
    % first select cases not too close to eachother (at least 10 days
    % apart)
    data = CLMAS.(farms{i})(~isnan(CLMAS.(farms{i}).B_ID),:);
    data.diff(1,1) = 10;
    data.diff(2:end,1) = diff(datenum(data.Date));
    [~,indx] = unique(data.B_ID); % find first occurence
    data.first(indx) = 1;
    data(data.first == 0 & data.diff < 10,:) = []; % delete
    % re-add to CLMAS
    data = removevars(data,{'diff','first'});
    CLMAS.(farms{i}) = data; 
    clear data
    
    ind = find(~isnan(CLMAS.(farms{i}).B_ID));
        
    L = length(ind); % 
    ML.TRP.FarmName(T:T+L-1,1) = repmat(farms(i),L,1);
    ML.TRP.B_ID(T:T+L-1) = CLMAS.(farms{i}).B_ID(ind);
    ML.TRP.Date(T:T+L-1) = CLMAS.(farms{i}).Date(ind);
    for ii = 1:length(ind)
        % find corresponding lactation and DIM
        try
            % find the corresponding data in a lactation within 305 dim
            idxmod = find(moddataDAY.(farms{i}).B_ID == ML.TRP.B_ID(T+ii-1) & ...
                       floor(datenum(moddataDAY.(farms{i}).Date)) == floor(datenum(ML.TRP.Date(T+ii-1))));  
            ML.TRP.Lac(T+ii-1) = moddataDAY.(farms{i}).Lac(idxmod(1));
            ML.TRP.DIM(T+ii-1) = moddataDAY.(farms{i}).DIM(idxmod(1));
            idx = find(PT_zeroDAY.(farms{i}).B_ID == ML.TRP.B_ID(T+ii-1) & ...
                       PT_zeroDAY.(farms{i}).Lac == ML.TRP.Lac(T+ii-1) & ...
                       PT_zeroDAY.(farms{i}).DIMStart <= ML.TRP.DIM(T+ii-1) & ...
                       PT_zeroDAY.(farms{i}).DIMEnd >= ML.TRP.DIM(T+ii-1));
        catch
            % fill in empty
            idx = [];
        end
        % fill in
        if ~isempty(idx)
            % add Q
            ML.TRP.Q(T+ii-1) = CLMAS.(farms{i}).ms(ind(ii));
            % add length
            ML.TRP.PertLength(T+ii-1) = PT_zeroDAY.(farms{i}).PertLength(idx(1));
            % add losses
            ML.TRP.MLa(T+ii-1) = PT_zeroDAY.(farms{i}).MLa(idx(1));
            ML.TRP.MLr(T+ii-1) = PT_zeroDAY.(farms{i}).MLr(idx(1));
            ML.TRP.MLl(T+ii-1) = PT_zeroDAY.(farms{i}).MLl(idx(1));
            ML.TRP.Start(T+ii-1) = PT_zeroDAY.(farms{i}).DIMStart(idx(1));
            ML.TRP.End(T+ii-1) = PT_zeroDAY.(farms{i}).DIMEnd(idx(1));
            
        else
            ML.TRP.Lac(T+ii-1) = NaN;
            ML.TRP.DIM(T+ii-1) = NaN;
        end
    end
    T = T+L;
end
clear sub T L i ii idx ind data DIM ITW ans p test h

disp(['No. of unmatched cases = ' num2str(length(find(isnan(ML.TRP.Lac)))) ...
       ', = ' num2str(length(find(isnan(ML.TRP.Lac)))./height(ML.TRP)*100) '%'])
disp(['No. of matched cases = ' num2str(length(find(~isnan(ML.TRP.Lac)))) ...
       ', = ' num2str(length(find(~isnan(ML.TRP.Lac)))./height(ML.TRP)*100) '%'])

% delete unmatched cases
ML.TRP = ML.TRP(~isnan(ML.TRP.DIM),:);


% figure perturbation length
% visualize milk losses
ind11 = find(ML.TRP.Lac == 1 & ML.TRP.DIM <= 63);
ind12 = find(ML.TRP.Lac == 1 & ML.TRP.DIM > 63 & ML.TRP.DIM <= 138);
ind13 = find(ML.TRP.Lac == 1 & ML.TRP.DIM > 138 & ML.TRP.DIM <= 216);
ind14 = find(ML.TRP.Lac == 1 & ML.TRP.DIM > 216);
ind21 = find(ML.TRP.Lac > 1 & ML.TRP.DIM <= 63);
ind22 = find(ML.TRP.Lac > 1 & ML.TRP.DIM > 63 & ML.TRP.DIM <= 138);
ind23 = find(ML.TRP.Lac > 1 & ML.TRP.DIM > 138 & ML.TRP.DIM <= 216);
ind24 = find(ML.TRP.Lac > 1 & ML.TRP.DIM > 216);

% 
figure('Units','Normalized','OuterPosition',[0.1 0.1 0.85 0.8]);
binedges = 5:15:150;
        subplot(2,4,1);  h = histogram(ML.TRP.PertLength(ind11),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0.8 1],...
                        'FaceColor',[0.2 0.8 1]);
                    x = h.BinEdges;
                    y = round(h.Values*100);
                    text(x(1:end-1)+7.3,y/100,num2str(y'),'vert','bottom','horiz','center')
        subplot(2,4,2); h=histogram(ML.TRP.PertLength(ind12),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0.6 1],...
                        'FaceColor',[0.2 0.6 1]);
                    x = h.BinEdges;
                    y = round(h.Values*100);
                    text(x(1:end-1)+7.3,y/100,num2str(y'),'vert','bottom','horiz','center')
        subplot(2,4,3); h=histogram(ML.TRP.PertLength(ind13),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0 1],...
                        'FaceColor',[0.2 0 1]);               
                    x = h.BinEdges;
                    y = round(h.Values*100);
                    text(x(1:end-1)+7.3,y/100,num2str(y'),'vert','bottom','horiz','center')
        subplot(2,4,4); h=histogram(ML.TRP.PertLength(ind14),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0 0 0.6],...
                        'FaceColor',[0 0 0.6]);
                    x = h.BinEdges;
                    y = round(h.Values*100);
                    text(x(1:end-1)+7.3,y/100,num2str(y'),'vert','bottom','horiz','center')
        subplot(2,4,5); h=histogram(ML.TRP.PertLength(ind21),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0.8 1],...
                        'FaceColor',[0.2 0.8 1]);
                    x = h.BinEdges;
                    y = round(h.Values*100);
                    text(x(1:end-1)+7.3,y/100,num2str(y'),'vert','bottom','horiz','center')
        subplot(2,4,6); h=histogram(ML.TRP.PertLength(ind22),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0.6 1],...
                        'FaceColor',[0.2 0.6 1]);
                    x = h.BinEdges;
                    y = round(h.Values*100);
                    text(x(1:end-1)+7.3,y/100,num2str(y'),'vert','bottom','horiz','center')
        subplot(2,4,7); h=histogram(ML.TRP.PertLength(ind23),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.2 0 1],...
                        'FaceColor',[0.2 0 1]);
                    x = h.BinEdges;
                    y = round(h.Values*100);
                    text(x(1:end-1)+7.3,y/100,num2str(y'),'vert','bottom','horiz','center')
        subplot(2,4,8); h=histogram(ML.TRP.PertLength(ind24),binedges,'Normalization','probability',...
                        'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0 0 0.6],...
                        'FaceColor',[0 0 0.6]);
                    x = h.BinEdges;
                    y = round(h.Values*100);
                    text(x(1:end-1)+7.3,y/100,num2str(y'),'vert','bottom','horiz','center')

    
pl.titles = {'Lac = 1, LS = 1','Lac = 1, LS = 2','Lac = 1, LS = 3','Lac = 1, LS = 4',...
             'Lac > 1, LS = 1','Lac > 1, LS = 2','Lac > 1, LS = 3','Lac > 1, LS = 4'};
pl.xlabel = {'Perturbation length (days)'};      
pl.ylim = [0 0.75];
for i = 1:8
    subplot(2,4,i);
    title(pl.titles{i})
    xlabel(pl.xlabel)
    yticks(0:0.1:0.70)
    yticklabels({'0','10','20','30','40','50','60','70'})
    ylim(pl.ylim)
    hold on; plot([30 30], [0 1],'r-.','LineWidth',1.4);%xlim([-150 60]);
end

figure('Units','Normalized','OuterPosition',[0.1 0.1 0.35 0.48]);
binedges = 5:15:150;
h = histogram(ML.TRP.PertLength,binedges,'Normalization','probability',...
              'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.6 0 0.6],...
              'FaceColor',[0.6 0 0.6]);
               x = h.BinEdges;
               y = round(h.Values*100);
               ypos = h.Values;
               text(x(1:end-1)+7.3,ypos,num2str(y'),'vert','bottom','horiz','center')
title('Perturbation length around clinical mastitis')
xlabel('Length (days)'); ylabel('% cases')
yticks(0:0.1:0.70);yticklabels({'0','10','20','30','40','50','60','70'})


figure('Units','Normalized','OuterPosition',[0.1 0.1 0.55 0.48]);
subplot(1,2,1);
binedges = 0:5:40;
h = histogram(ML.TRP.DIM-ML.TRP.Start,binedges,'Normalization','probability',...
              'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.6 0 0.6],...
              'FaceColor',[0.6 0 0.6]);
               x = h.BinEdges;
               y = round(h.Values*100);
               ypos = h.Values;
               text(x(1:end-1)+2.3,ypos,num2str(y'),'vert','bottom','horiz','center')
title('Start perturbation before treatment DIM')
xlabel('No. of days'); ylabel('% cases')
yticks(0:0.1:0.60);yticklabels({'0','10','20','30','40','50','60'})
subplot(1,2,2);
binedges = 0:10:100;
h = histogram(ML.TRP.End-ML.TRP.DIM,binedges,'Normalization','probability',...
              'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.8 0 0.2],...
              'FaceColor',[0.8 0 0.2]);
               x = h.BinEdges;
               y = round(h.Values*100);
               ypos = h.Values;
               text(x(1:end-1)+5,ypos,num2str(y'),'vert','bottom','horiz','center')
title('End perturbation after treatment DIM')
xlabel('No. of days'); ylabel('% cases')
yticks(0:0.1:0.4);yticklabels({'0','10','20','30','40'})
ylim([0 0.40])

% average milk loss
disp(['Average milk loss = ' num2str(mean(ML.TRP.MLa)) ' +/- ' num2str(std(ML.TRP.MLa))])

% correlation length and milk loss absolute
figure('Units','Normalized','OuterPosition',[0.1 0.1 0.30 0.5]);
        box on; yyaxis left; ylabel('Milk loss (kg)');ylim([0 1500]); 
        plot(ML.TRP.PertLength,ML.TRP.MLa,'s','MarkerSize',3,'Color',[0.2 0 1])
        yyaxis right; plot(ML.TRP.PertLength,ML.TRP.MLr,'o','MarkerSize',3,'Color',[0.8 0 0.6])
        hold on; plot(ML.TRP.PertLength,ML.TRP.MLl,'o','MarkerSize',3)
        xlabel('Perturbation length (days)'); ylabel('% loss');
        legend({'Absolute losses (kg)','Relative loss % day','Relative losses % lac'},'AutoUpdate','off')
x = sortrows([ML.TRP.PertLength,ones(length(ML.TRP.MLa),1) ML.TRP.MLa],1);
        b = regress(x(:,3),x(:,[2 1]));
yyaxis left; plot(5:180,b(1)+b(2)*(5:180),'LineWidth',3.5,'Color',[0.2 0 1])
x = sortrows([ML.TRP.PertLength,ones(length(ML.TRP.MLl),1) ML.TRP.MLl],1);
        b = regress(x(:,3),x(:,[2 1]));
yyaxis right; plot(5:180,b(1)+b(2)*(5:180),'LineWidth',3.5,'Color',[1 0.4 0])



% histograms milk losses
figure('Units','Normalized','OuterPosition',[0.1 0.1 0.50 0.4]);
subplot(1,2,1); h = histogram(-ML.TRP.MLa,'Normalization','probability',...
              'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.8 0 0.2],...
              'FaceColor',[0.8 0 0.2]);
               x = h.BinEdges;
               y = round(h.Values*100);
               ypos = h.Values;
               text(x(1:end-1)+50,ypos,num2str(y'),'vert','bottom','horiz','center')
xlim([-1200 10]); title('Total absolute milk loss per perturbation (kg)')
ylim([0 0.45]); yticks(0:0.1:0.4); yticklabels({'0','10','20','30','40'}); ylabel('% cases')
xlabel('Milk loss (kg)')
          subplot(1,2,2); h = histogram(ML.TRP.MLr,'Normalization','probability',...
              'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.8 0 0.2],...
              'FaceColor',[0.8 0 0.2]);
               x = h.BinEdges;
               y = round(h.Values*100);
               ypos = h.Values;
               text(x(1:end-1)+1.5,ypos,num2str(y'),'vert','bottom','horiz','center')
        xlim([0 45]); title('Relative daily milk loss per perturbation (%)')
ylim([0 0.25]); yticks(0:0.05:0.25); yticklabels({'0','5','10','15','20','25'}); ylabel('% cases')
xlabel('Milk loss (%/d)')



%% Quarter perturbations treatments

farms = fieldnames(CLMAS);  % farms with CM data
ML.TRPQ = table(zeros(0,1),'VariableNames',{'B_ID'}); % prepare case table
ML.TRPQ.FarmName = repmat({''},0);
ML.TRPQ = movevars(ML.TRPQ, 'FarmName','Before','B_ID');
mod.Wood = @(p,t) p(1).*t.^p(2).*exp(-p(3).*t);     % wood model - prepare
T = 1;  % counter to fill table
for i = 1:length(farms)
    
    % first select cases not too close to eachother (at least 10 days
    % apart)
    data = CLMAS.(farms{i})(~isnan(CLMAS.(farms{i}).B_ID),:);
    data.diff(1,1) = 10;
    data.diff(2:end,1) = diff(datenum(data.Date));
    [~,indx] = unique(data.B_ID); % find first occurence
    data.first(indx) = 1;
    data(data.first == 0 & data.diff < 10,:) = []; % delete
    % re-add to CLMAS
    data = removevars(data,{'diff','first'});
    CLMAS.(farms{i}) = data; 
    clear data
    
    ind = find(~isnan(CLMAS.(farms{i}).B_ID));
        
    L = length(ind); % 
    ML.TRPQ.FarmName(T:T+L-1,1) = repmat(farms(i),L,1);
    ML.TRPQ.B_ID(T:T+L-1) = CLMAS.(farms{i}).B_ID(ind);
    ML.TRPQ.Date(T:T+L-1) = CLMAS.(farms{i}).Date(ind);
    
    for ii = 1:length(ind)
        % find corresponding lactation and DIM
        try
            % find the corresponding data in a lactation within 305 dim
            idxmod = find(moddata.(farms{i}).B_ID == ML.TRPQ.B_ID(T+ii-1) & ...
                       floor(datenum(moddata.(farms{i}).Date)) == floor(datenum(ML.TRPQ.Date(T+ii-1))));     
            % add lac and dim
            ML.TRPQ.Lac(T+ii-1) = moddata.(farms{i}).Lac(idxmod(1));
            ML.TRPQ.DIM(T+ii-1) = moddata.(farms{i}).DIM(idxmod(1));
            
            % find corresponding perturbations in all     
            idx = find(PT_all.(farms{i}).B_ID == ML.TRPQ.B_ID(T+ii-1) & ...
                       PT_all.(farms{i}).Lac == ML.TRPQ.Lac(T+ii-1) & ...
                       ((PT_all.(farms{i}).DIMstart1 <= ML.TRPQ.DIM(T+ii-1) & PT_all.(farms{i}).DIMend1 >= ML.TRPQ.DIM(T+ii-1))|...
                        (PT_all.(farms{i}).DIMstart2 <= ML.TRPQ.DIM(T+ii-1) & PT_all.(farms{i}).DIMend2 >= ML.TRPQ.DIM(T+ii-1))|...
                        (PT_all.(farms{i}).DIMstart3 <= ML.TRPQ.DIM(T+ii-1) & PT_all.(farms{i}).DIMend3 >= ML.TRPQ.DIM(T+ii-1))|...
                        (PT_all.(farms{i}).DIMstart4 <= ML.TRPQ.DIM(T+ii-1) & PT_all.(farms{i}).DIMend4 >= ML.TRPQ.DIM(T+ii-1))));
       
        catch
            % fill in empty
            idx = [];
        end
        
        % fill in perturbation data
        if ~isempty(idx)
            % add Q
            ML.TRPQ.Q(T+ii-1) = CLMAS.(farms{i}).ms(ind(ii));
            
            % add start and end
            ML.TRPQ.Start1(T+ii-1) = PT_all.(farms{i}).DIMstart1(idx(1));
            ML.TRPQ.Start2(T+ii-1) = PT_all.(farms{i}).DIMstart2(idx(1));
            ML.TRPQ.Start3(T+ii-1) = PT_all.(farms{i}).DIMstart3(idx(1));
            ML.TRPQ.Start4(T+ii-1) = PT_all.(farms{i}).DIMstart4(idx(1));

            ML.TRPQ.End1(T+ii-1) = PT_all.(farms{i}).DIMend1(idx(1));
            ML.TRPQ.End2(T+ii-1) = PT_all.(farms{i}).DIMend2(idx(1));
            ML.TRPQ.End3(T+ii-1) = PT_all.(farms{i}).DIMend3(idx(1));
            ML.TRPQ.End4(T+ii-1) = PT_all.(farms{i}).DIMend4(idx(1));
                        
            % add length
            ML.TRPQ.PertLength1(T+ii-1) = PT_all.(farms{i}).DIMend1(idx(1)) - PT_all.(farms{i}).DIMstart1(idx(1));
            ML.TRPQ.PertLength2(T+ii-1) = PT_all.(farms{i}).DIMend2(idx(1)) - PT_all.(farms{i}).DIMstart2(idx(1));
            ML.TRPQ.PertLength3(T+ii-1) = PT_all.(farms{i}).DIMend3(idx(1)) - PT_all.(farms{i}).DIMstart3(idx(1));
            ML.TRPQ.PertLength4(T+ii-1) = PT_all.(farms{i}).DIMend4(idx(1)) - PT_all.(farms{i}).DIMstart4(idx(1));
            
            % add losses
            ML.TRPQ.ML1a(T+ii-1) = PT_all.(farms{i}).ML1a(idx(1));
            ML.TRPQ.ML2a(T+ii-1) = PT_all.(farms{i}).ML2a(idx(1));
            ML.TRPQ.ML3a(T+ii-1) = PT_all.(farms{i}).ML3a(idx(1));
            ML.TRPQ.ML4a(T+ii-1) = PT_all.(farms{i}).ML4a(idx(1));
                  
            ML.TRPQ.ML1r(T+ii-1) = PT_all.(farms{i}).ML1r(idx(1));
            ML.TRPQ.ML2r(T+ii-1) = PT_all.(farms{i}).ML2r(idx(1));
            ML.TRPQ.ML3r(T+ii-1) = PT_all.(farms{i}).ML3r(idx(1));
            ML.TRPQ.ML4r(T+ii-1) = PT_all.(farms{i}).ML4r(idx(1));
            
            ML.TRPQ.ML1l(T+ii-1) = PT_all.(farms{i}).ML1l(idx(1));
            ML.TRPQ.ML2l(T+ii-1) = PT_all.(farms{i}).ML2l(idx(1));
            ML.TRPQ.ML3l(T+ii-1) = PT_all.(farms{i}).ML3l(idx(1));
            ML.TRPQ.ML4l(T+ii-1) = PT_all.(farms{i}).ML4l(idx(1));
            
            ML.TRPQ.smML1a(T+ii-1) = PT_all.(farms{i}).smML1a(idx(1));
            ML.TRPQ.smML2a(T+ii-1) = PT_all.(farms{i}).smML2a(idx(1));
            ML.TRPQ.smML3a(T+ii-1) = PT_all.(farms{i}).smML3a(idx(1));
            ML.TRPQ.smML4a(T+ii-1) = PT_all.(farms{i}).smML4a(idx(1));
                  
            ML.TRPQ.smML1r(T+ii-1) = PT_all.(farms{i}).smML1r(idx(1));
            ML.TRPQ.smML2r(T+ii-1) = PT_all.(farms{i}).smML2r(idx(1));
            ML.TRPQ.smML3r(T+ii-1) = PT_all.(farms{i}).smML3r(idx(1));
            ML.TRPQ.smML4r(T+ii-1) = PT_all.(farms{i}).smML4r(idx(1));
            
            ML.TRPQ.smML1l(T+ii-1) = PT_all.(farms{i}).smML1l(idx(1));
            ML.TRPQ.smML2l(T+ii-1) = PT_all.(farms{i}).smML2l(idx(1));
            ML.TRPQ.smML3l(T+ii-1) = PT_all.(farms{i}).smML3l(idx(1));
            ML.TRPQ.smML4l(T+ii-1) = PT_all.(farms{i}).smML4l(idx(1));
            
            % index of pert at day level
            idxday = find(ML.TRP.B_ID == ML.TRPQ.B_ID(T+ii-1) & ...
                          ML.TRP.Lac == ML.TRPQ.Lac(T+ii-1) & ...
                          contains(ML.TRP.FarmName,ML.TRPQ.FarmName{T+ii-1}) & ...
                          ((ML.TRP.DIM >= ML.TRPQ.Start1(T+ii-1) & ML.TRP.DIM <= ML.TRPQ.End1(T+ii-1)) | ...
                           (ML.TRP.DIM >= ML.TRPQ.Start2(T+ii-1) & ML.TRP.DIM <= ML.TRPQ.End2(T+ii-1)) | ...
                           (ML.TRP.DIM >= ML.TRPQ.Start3(T+ii-1) & ML.TRP.DIM <= ML.TRPQ.End3(T+ii-1)) | ...
                           (ML.TRP.DIM >= ML.TRPQ.Start4(T+ii-1) & ML.TRP.DIM <= ML.TRPQ.End4(T+ii-1))));
                       
            if length(idxday) == 1
                ML.TRPQ.Startday(T+ii-1) = ML.TRP.Start(idxday(1));
                ML.TRPQ.Endday(T+ii-1) = ML.TRP.End(idxday(1));
                ML.TRPQ.Lengthday(T+ii-1) = ML.TRP.PertLength(idxday(1));
                ML.TRPQ.MLaday(T+ii-1) = ML.TRP.MLa(idxday(1));
                ML.TRPQ.MLrday(T+ii-1) = ML.TRP.MLr(idxday(1));
                ML.TRPQ.MLlday(T+ii-1) = ML.TRP.MLl(idxday(1));
            elseif length(idxday) == 1
                ML.TRPQ.Startday2(T+ii-1) = ML.TRP.Start(idxday(2));
                ML.TRPQ.Endday2(T+ii-1) = ML.TRP.End(idxday(2));
                ML.TRPQ.Lengthday2(T+ii-1) = ML.TRP.PertLength(idxday(2));
                ML.TRPQ.MLaday2(T+ii-1) = ML.TRP.MLa(idxday(2));
                ML.TRPQ.MLrday2(T+ii-1) = ML.TRP.MLr(idxday(2));
                ML.TRPQ.MLlday2(T+ii-1) = ML.TRP.MLl(idxday(2));
            end
                
                
        else
            ML.TRPQ.Lac(T+ii-1) = NaN;
            ML.TRPQ.DIM(T+ii-1) = NaN;
        end
    end
    T = T+L;
end
clear sub T L i ii idx ind data DIM ITW ans p test h

disp(['No. of unmatched cases = ' num2str(length(find(isnan(ML.TRPQ.Lac)))) ...
       ', = ' num2str(length(find(isnan(ML.TRPQ.Lac)))./height(ML.TRPQ)*100) '%'])
disp(['No. of matched cases = ' num2str(length(find(~isnan(ML.TRPQ.Lac)))) ...
       ', = ' num2str(length(find(~isnan(ML.TRPQ.Lac)))./height(ML.TRPQ)*100) '%'])


% delete unmatched cases
ML.TRPQ = ML.TRPQ(~isnan(ML.TRPQ.DIM),:);

% delete cases with no Q info
ind = find(nansum(ML.TRQ{:,9:12},2)==0);  % 158 cases zero q yield
ML.TRQ(ind,:) = [];


% no. no day perturbation

disp(['No. of day perturbations = ' num2str(length(find(ML.TRPQ.Startday ~= 0))) ...
       ', = ' num2str(length(find(ML.TRPQ.Startday ~= 0))./height(ML.TRPQ)*100) '%'])





% figure no of quarters perturbation detected
for i = 1:height(ML.TRPQ)
    ML.TRPQ.NoQuart(i) = length(find(~isnan(ML.TRPQ{i,7:10})));
end

figure('Units','Normalized','OuterPosition',[0.1 0.1 0.250 0.4]);
 h = histogram(ML.TRPQ.NoQuart,'Normalization','probability',...
              'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.8 0 0.2],...
              'FaceColor',[0.8 0 0.2]);
               x = h.BinEdges;
               y = round(h.Values*100);
               ypos = h.Values;
               text(x(1:end-1)+0.50,ypos,num2str(y'),'vert','bottom','horiz','center')
title('No. of quarters with perturbation')
ylim([0 0.4]); yticks(0:0.1:0.4); yticklabels({'0','10','20','30','40'}); ylabel('% cases')
xticks(1:4); xticklabels({'1','2','3','4'}); xlabel('No. of quarters')

% assumption: longest perturbation = infected quarter
for i = 1:height(ML.TRPQ)
    ind = find(ML.TRPQ{i,15:18} == max(ML.TRPQ{i,15:18}));
    indn = find(ML.TRPQ{i,15:18} ~= max(ML.TRPQ{i,15:18}));
    if length(ind) == 1
        ML.TRPQ.susQ(i) = ind;
        ML.TRPQ.nonsusQ1(i) = indn(1);
        ML.TRPQ.nonsusQ2(i) = indn(2);
        ML.TRPQ.nonsusQ3(i) = indn(3);
    else
        ML.TRPQ.susQ(i) = NaN;
        ML.TRPQ.nonsusQ1(i) = NaN;
        ML.TRPQ.nonsusQ2(i) = NaN;
        ML.TRPQ.nonsusQ3(i) = NaN;
    end
end

figure('Units','Normalized','OuterPosition',[0.1 0.1 0.250 0.4]);
 h = histogram(ML.TRPQ.susQ,'Normalization','probability',...
              'DisplayStyle','bar','LineWidth',2,'EdgeColor',[0.8 0 0.2],...
              'FaceColor',[0.8 0 0.2]);
               x = h.BinEdges;
               y = round(h.Values*100);
               ypos = h.Values;
               text(x(1:end-1)+0.50,ypos,num2str(y'),'vert','bottom','horiz','center')
title('Suspected Quarter')
ylim([0 0.3]); yticks(0:0.1:0.4); yticklabels({'0','10','20','30','40'}); ylabel('% cases')
xticks(1:4); xticklabels({'LF','RF','LH','RH'});



% plot milk losses in suspected vs in all quarters
nifq = []; ifq = [];

ML.TRPQ.PertLength1(isnan(ML.TRPQ.PertLength1)) = 0;
ML.TRPQ.PertLength2(isnan(ML.TRPQ.PertLength2)) = 0;
ML.TRPQ.PertLength3(isnan(ML.TRPQ.PertLength3)) = 0;
ML.TRPQ.PertLength4(isnan(ML.TRPQ.PertLength4)) = 0;


% Q1 = suspected
ind = find(ML.TRPQ.susQ == 1);
ifq(ind,1:3) = [1*ones(length(ind),1) ML.TRPQ.DIM(ind) ML.TRPQ.Lac(ind)];
ifq(ind,4) = ML.TRPQ.ML1a(ind);
ifq(ind,5) = ML.TRPQ.ML1r(ind);
ifq(ind,6) = ML.TRPQ.PertLength1(ind);
nifq(ind,1:3) = [1*ones(length(ind),1) ML.TRPQ.DIM(ind) ML.TRPQ.Lac(ind)];
nifq(ind,5) = nansum([ ML.TRPQ.ML2a(ind),ML.TRPQ.ML3a(ind),ML.TRPQ.ML4a(ind)],2)./3;
nifq(ind,6) = nansum([ ML.TRPQ.ML2a(ind),ML.TRPQ.ML3a(ind),ML.TRPQ.ML4a(ind)],2);
nifq(ind,7) = nanmean([ML.TRPQ.ML2r(ind),ML.TRPQ.ML3r(ind),ML.TRPQ.ML4r(ind)],2);
nifq(ind,8) = nanmean([ML.TRPQ.PertLength2(ind),ML.TRPQ.PertLength3(ind),ML.TRPQ.PertLength4(ind)],2);
% Q2 = suspected
ind = find(ML.TRPQ.susQ == 2);
ifq(ind,1:3) = [2*ones(length(ind),1) ML.TRPQ.DIM(ind) ML.TRPQ.Lac(ind)];
ifq(ind,4) = ML.TRPQ.ML2a(ind);
ifq(ind,5) = ML.TRPQ.ML2r(ind);
ifq(ind,6) = ML.TRPQ.PertLength2(ind);
nifq(ind,1:3) = [2*ones(length(ind),1) ML.TRPQ.DIM(ind) ML.TRPQ.Lac(ind)];
nifq(ind,5) = nansum([ ML.TRPQ.ML1a(ind),ML.TRPQ.ML3a(ind),ML.TRPQ.ML4a(ind)],2)./3;
nifq(ind,6) = nansum([ ML.TRPQ.ML1a(ind),ML.TRPQ.ML3a(ind),ML.TRPQ.ML4a(ind)],2);
nifq(ind,7) = nanmean([ML.TRPQ.ML1r(ind),ML.TRPQ.ML3r(ind),ML.TRPQ.ML4r(ind)],2);
nifq(ind,8) = nanmean([ML.TRPQ.PertLength1(ind),ML.TRPQ.PertLength3(ind),ML.TRPQ.PertLength4(ind)],2);
% Q3 = suspected
ind = find(ML.TRPQ.susQ == 3);
ifq(ind,1:3) = [3*ones(length(ind),1) ML.TRPQ.DIM(ind) ML.TRPQ.Lac(ind)];
ifq(ind,4) = ML.TRPQ.ML3a(ind);
ifq(ind,5) = ML.TRPQ.ML3r(ind);
ifq(ind,6) = ML.TRPQ.PertLength3(ind);
nifq(ind,1:3) = [3*ones(length(ind),1) ML.TRPQ.DIM(ind) ML.TRPQ.Lac(ind)];
nifq(ind,5) = nansum([ ML.TRPQ.ML2a(ind),ML.TRPQ.ML1a(ind),ML.TRPQ.ML4a(ind)],2)./3;
nifq(ind,6) = nansum([ ML.TRPQ.ML2a(ind),ML.TRPQ.ML1a(ind),ML.TRPQ.ML4a(ind)],2);
nifq(ind,7) = nanmean([ML.TRPQ.ML2r(ind),ML.TRPQ.ML1r(ind),ML.TRPQ.ML4r(ind)],2);
nifq(ind,8) = nanmean([ML.TRPQ.PertLength2(ind),ML.TRPQ.PertLength1(ind),ML.TRPQ.PertLength4(ind)],2);
% Q4 = suspected
ind = find(ML.TRPQ.susQ == 4);
ifq(ind,1:3) = [4*ones(length(ind),1) ML.TRPQ.DIM(ind) ML.TRPQ.Lac(ind)];
ifq(ind,4) = ML.TRPQ.ML4a(ind);
ifq(ind,5) = ML.TRPQ.ML4r(ind);
ifq(ind,6) = ML.TRPQ.PertLength4(ind);
nifq(ind,1:3) = [4*ones(length(ind),1) ML.TRPQ.DIM(ind) ML.TRPQ.Lac(ind)];
nifq(ind,5) = nansum([ ML.TRPQ.ML2a(ind),ML.TRPQ.ML3a(ind),ML.TRPQ.ML1a(ind)],2)./3;
nifq(ind,6) = nansum([ ML.TRPQ.ML2a(ind),ML.TRPQ.ML3a(ind),ML.TRPQ.ML1a(ind)],2);
nifq(ind,7) = nanmean([ML.TRPQ.ML2r(ind),ML.TRPQ.ML3r(ind),ML.TRPQ.ML1r(ind)],2);
nifq(ind,8) = nanmean([ML.TRPQ.PertLength2(ind),ML.TRPQ.PertLength3(ind),ML.TRPQ.PertLength1(ind)],2);

nifq(:,4) = [];
ifq = array2table(ifq,'VariableNames',{'susQ','DIM','Lac','MLa','MLr','Length'});
nifq = array2table(nifq,'VariableNames',{'susQ','DIM','Lac','MLaaverage','MLatot','MLraverage','Lengthaverage'});


figure; hold on; 
subplot(1,3,1); hold on; box on; 
binedges = -600:50:50;
h = histogram(-nifq.MLaaverage,binedges,'Normalization','probability','LineWidth',2,'EdgeColor',[0.2 0.6 0.2],...
              'FaceColor',[0.2 0.6 0.2]);
               x = h.BinEdges;
               y = round(h.Values*100);
               ypos = h.Values;
               text(x(1:end-1)+25,ypos,num2str(y'),'vert','bottom','horiz','center')

h = histogram(-ifq.MLa,-600:50:0,'Normalization','probability','LineWidth',2,'EdgeColor',[0.8 0 0.2],...
              'FaceColor',[0.8 0 0.2]);
               x = h.BinEdges;
               y = round(h.Values*100);
               ypos = h.Values;
               text(x(1:end-1)+25,ypos,num2str(y'),'vert','bottom','horiz','center')
xlim([-610,60]); title('Absolute milk losses per perturbation')
legend({'Average non-infected Q','Infected Q'},'Location','northwest')
ylim([0 0.55]); yticks(0:0.1:0.5); yticklabels({'0','10','20','30','40','50'}); ylabel('% cases')

subplot(1,3,2); hold on; box on; 
binedges = -100:10:10;
h = histogram(-nifq.MLraverage,binedges,'Normalization','probability','LineWidth',2,'EdgeColor',[0.2 0.6 0.2],...
              'FaceColor',[0.2 0.6 0.2]);
               x = h.BinEdges;
               y = round(h.Values*100);
               ypos = h.Values;
               text(x(1:end-1)+5,ypos,num2str(y'),'vert','bottom','horiz','center')

h = histogram(-ifq.MLr,-100:10:0,'Normalization','probability','LineWidth',2,'EdgeColor',[0.8 0 0.2],...
              'FaceColor',[0.8 0 0.2]);
               x = h.BinEdges;
               y = round(h.Values*100);
               ypos = h.Values;
               text(x(1:end-1)+5,ypos,num2str(y'),'vert','bottom','horiz','center')
xlim([-105,15]); title('Relative milk loss per day')
legend({'Average non-infected Q','Infected Q'},'Location','northwest')
ylim([0 0.40]); yticks(0:0.1:0.5); yticklabels({'0','10','20','30','40','50'}); ylabel('% cases')

subplot(1,3,3); hold on; box on; 
binedges = 0:10:100;
h = histogram(nifq.Lengthaverage,binedges,'Normalization','probability','LineWidth',2,'EdgeColor',[0.2 0.6 0.2],...
              'FaceColor',[0.2 0.6 0.2]);
               x = h.BinEdges;
               y = round(h.Values*100);
               ypos = h.Values;
               text(x(1:end-1)+5,ypos,num2str(y'),'vert','bottom','horiz','center')

h = histogram(ifq.Length,binedges,'Normalization','probability','LineWidth',2,'EdgeColor',[0.8 0 0.2],...
              'FaceColor',[0.8 0 0.2]);
               x = h.BinEdges;
               y = round(h.Values*100);
               ypos = h.Values;
               text(x(1:end-1)+5,ypos,num2str(y'),'vert','bottom','horiz','center')
xlim([-5,105]); title('Perturbation length')
legend({'Average non-infected Q','Infected Q'},'Location','northeast')
ylim([0 0.55]); yticks(0:0.1:0.5); yticklabels({'0','10','20','30','40','50'}); ylabel('% cases')


mean(ifq.MLa)
std(ifq.MLa)
mean(ifq.MLr)
std(ifq.MLr)
mean(ifq.Length)
std(ifq.Length)

mean(nifq.MLaaverage)
std(nifq.MLaaverage)
mean(nifq.MLraverage)
std(nifq.MLraverage)
mean(nifq.Lengthaverage)
std(nifq.Lengthaverage)

% save
save('D10_ML.mat','ML')

clear ind11 ind12 ind13 ind14 ind21 ind22 ind23 ind24 i idxday idxmoc
clear ans b binedges h farms indn indx nifq smML x y ypos ifq idxmod ind
clear PT_all PT_perc PT_zero PT_zeroDAY PT_zeroDAY2 mod

%% MPR DATA
% milk yield residuals ??? 10 ??? days before and after a higher SCC
%   measurement
% proportion of high cell count sample corresponding to a perturbation
%   according to the criteria set up for the QD
% milk losses in the perturbations corresponding to the high SCC
% proportion of cases in which a difference between quarters is differences between quarters when 

load('D8_MPRcases.mat')
moddataDAY.DeSchutter = unique([moddataDAY.Deschutter(:,[1:3 5:end]);moddataDAY.DeSchutter(:,[1:3 5:end])],'rows');
moddataDAY.DeSchutter.Name(:,1) = {''};
moddataDAY.DeSchutter = movevars(moddataDAY.DeSchutter,'Name','After','CowID');
moddataDAY = rmfield(moddataDAY,'Deschutter');

moddata.DeSchutter = unique([moddata.Deschutter(:,[1:3 5:end]);moddata.DeSchutter(:,[1:3 5:end])],'rows');
moddata.DeSchutter.Name(:,1) = {''};
moddata.DeSchutter = movevars(moddata.DeSchutter,'Name','After','CowID');
moddata = rmfield(moddataDAY,'Deschutter');

% match farm names
fn = fieldnames(CaseMPR);
farms = unique(MPR.FarmName);
farms([ 15 30 36]) = []; %  delete dewulf martens and agri-vet
farms2 = fieldnames(moddataDAY);
farms2(33) = []; % delete sanders
t = array2table(farms);
t.Var2(1) = farms2(1);
t.Var2(2) = farms2(16);
t.Var2(3) = farms2(35);
t.Var2(4) = farms2(11);
t.Var2(5:10) = farms2(2:7);
t.Var2(11) = farms2(12);
t.Var2(12) = farms2(40);
t.Var2(13:14) = farms2(9:10);
t.Var2(15) = farms2(14);
t.Var2(16) = farms2(15);
t.Var2(17) = farms2(8);
t.Var2(18:19) = farms2(17:18);
t.Var2(20) = farms2(13);
t.Var2(21:23) = farms2([19 21 23]);
t.Var2(24:34) = farms2(24:34);
t.Var2(35) = farms2(42);
t.Var2(36:37) = farms2(38:39);
t.Var2(38) = farms2(22);
t.Var2(39) = farms2(36);
t.Var2(40:47) = farms2([37 41 44 45 46 43 20 47]);

% SUM gives milk production
SUM(SUM.Select == 0,:) = [];

% calculate lactation production per cow and add to sum
for i = 1:height(SUM)
    % find overlap with day data
    farm = t.Var2(contains(t.farms,SUM.FarmName{i}));
    try
        ind = find(contains(moddataDAY.(farm{1}).OffRegNo,SUM.OffRegNo{i}) == 1 & ...
            moddataDAY.(farm{1}).Lac == SUM.Lac(i));
    catch
        ind = [];
    end
      
     if ~isempty(ind)
         data = moddataDAY.(farm{1})(ind,:);
         data.Res(:,1) = data.TDMY-data.ModIT;
         data.Res(1,1) = NaN;
         data.ModIT(1,1) = NaN;
         
         SUM.LacProd(i) = nansum(data.TDMY);
         SUM.LacProdSM(i) = nansum(data.MYs);
         SUM.LacPredict(i) = nansum(data.ModIT);
         SUM.TotalResid(i) = nansum(data.Res); 
         SUM.RelativeTotal(i) = nansum(data.Res)./nansum(data.ModIT)*100;
         SUM.RelativeDay(i) = nanmean(data.Res./data.ModIT)*100;
     else
         SUM.LacProd(i) = NaN;
         SUM.LacProdSM(i) = NaN;
         SUM.LacPredict(i) = NaN;
         SUM.TotalResid(i) = NaN;
         SUM.RelativeTotal(i) = NaN;
         SUM.RelativeDay(i) = NaN;
     end        
end

t.Farm(:,1) = 1:47;
SUM2 = innerjoin(SUM,t,'LeftKeys','FarmName','RightKeys','farms');
SUM2(isnan(SUM2.LacProd),:) = [];

% number of lactations in the dataset
figure;
h = histogram(SUM2.Farm(SUM2.Lac > 1),47);
x = h.BinEdges;
    y = h.Values;
    ypos = h.Values;ind = find(y == 0); y(ind) = []; x(ind) = [];
    text(x(1:end-1)+0.4,y,num2str(y'),'vert','bottom','horiz','center')
hold on; h = histogram(SUM2.Farm(SUM2.Lac == 1),47);
x = h.BinEdges;
    y = h.Values; ind = find(y == 0); y(ind) = []; x(ind) = [];
    text(x(1:end-1)+0.4,y,num2str(y'),'vert','bottom','horiz','center')
xticks = 0:50;
xlabel('Farm No.'); ylabel('No. of lactations'); legend({'Parity > 1','Parity = 1'});

% relative production
figure('Units','normalized','OuterPosition',[0.15 0.4 0.6 0.5]); 
subplot(1,2,1); histogram(SUM2.RelativeTotal(SUM2.Lac==1 & SUM2.MaxSCC < 200),-15:0.5:2.5,'Normalization','Probability',...
                          'LineWidth',2,'EdgeColor',[0.2 0.6 0.2],'FaceColor',[0.2 0.6 0.2]);
hold on; histogram(SUM2.RelativeTotal(SUM2.Lac==1 & SUM2.MaxSCC > 600),-15:0.5:2.5,'Normalization','Probability',...
                    'LineWidth',2,'EdgeColor',[0.8 0 0.2],'FaceColor',[0.8 0 0.2])
m = median(SUM2.RelativeTotal(SUM2.Lac==1 & SUM2.MaxSCC > 600));
    plot([m m],[0 0.20],'-.','LineWidth',2,'Color',[0.8 0 0.2])
m = median(SUM2.RelativeTotal(SUM2.Lac==1 & SUM2.MaxSCC <200));
    plot([m m],[0 0.20],'-.','LineWidth',2,'Color',[0.2 0.6 0.2])
subplot(1,2,2); histogram(SUM2.RelativeTotal(SUM2.Lac>1 & SUM2.MaxSCC < 200),-15:0.5:2.5,'Normalization','Probability',...
    'LineWidth',2,'EdgeColor',[0.2 0.6 0.2],'FaceColor',[0.2 0.6 0.2]);xlim([-15 3])
hold on; histogram(SUM2.RelativeTotal(SUM2.Lac>1 & SUM2.MaxSCC > 600),-15:0.5:2.5,'Normalization','Probability',...
                    'LineWidth',2,'EdgeColor',[0.8 0 0.2],'FaceColor',[0.8 0 0.2])
m = median(SUM2.RelativeTotal(SUM2.Lac>1 & SUM2.MaxSCC > 600));
    plot([m m],[0 0.20],'-.','LineWidth',2,'Color',[0.8 0 0.2])
m = median(SUM2.RelativeTotal(SUM2.Lac>1 & SUM2.MaxSCC < 600));
    plot([m m],[0 0.20],'-.','LineWidth',2,'Color',[0.2 0.6 0.2])
[ks.h,ks.p,ks.stats] = kstest2(SUM2.RelativeTotal(SUM2.Lac==1 & SUM2.MaxSCC < 200),SUM2.RelativeTotal(SUM2.Lac==1 & SUM2.MaxSCC > 600),'Tail','smaller');
[ks.h2,ks.p2,ks.stats2] = kstest2(SUM2.RelativeTotal(SUM2.Lac>1 & SUM2.MaxSCC < 200),SUM2.RelativeTotal(SUM2.Lac>1 & SUM2.MaxSCC > 600),'Tail','smaller');    
subplot(1,2,1); title(['LAC=1, Relative milk losses, p diff = ' num2str(ks.p)]); xlabel('Relative losses (%)'); ylabel('% lactations')
legend({'< 200 c/mL','>600 c/mL','median','median'},'Location','northwest')
subplot(1,2,2); title(['LAC>1, Relative milk losses, p diff = ' num2str(ks.p2)]); xlabel('Relative losses (%)'); ylabel('% lactations')
legend({'< 200 c/mL','>600 c/mL','median','median'},'Location','northwest')

% measurements of lactations for which my data is available
MPRi = innerjoin(MPR,SUM2(:,1:3),'Keys',{'FarmName','OffRegNo','Lac'});
SUM2.RelativeDay(SUM2.RelativeDay>25) = NaN;

% summaries per farm based on SUM2
sfarm = unique(SUM2(:,[20 19 1]),'rows');
for i = 1:height(sfarm)
    ind = find(SUM2.Farm == sfarm.Farm(i)); % all cows
    ind1 = find(SUM2.Farm == sfarm.Farm(i) & SUM2.Lac == 1); % heifers
    ind2 = find(SUM2.Farm == sfarm.Farm(i) & SUM2.Lac > 1); % cows
    ind3 = find(contains(MPRi.FarmName,sfarm.FarmName{i})); % all in MPR
    ind31 = find(contains(MPRi.FarmName,sfarm.FarmName{i}) & MPRi.Lac == 1); % all in MPR
    ind32 = find(contains(MPRi.FarmName,sfarm.FarmName{i}) & MPRi.Lac > 1); % all in MPR
   
    % all cows
    sfarm.NoCows(i,1) = length(ind);
    sfarm.NoMeas(i,1) = length(ind3);
    sfarm.AvSCCall(i,1) = nanmean(MPRi.SCC(ind3));
    sfarm.AverageNoMeas(i) = mean(SUM2.NoMeas(ind));
    sfarm.AvPercHigh(i) = (mean(SUM2.No150(ind1)./SUM2.NoMeas(ind1))*length(ind1) + ...
                          mean(SUM2.No250(ind2)./SUM2.NoMeas(ind2))*length(ind2))/(length(ind1)+length(ind2))*100;
    sfarm.LacProd(i) = mean(SUM2.LacProd(ind));
    sfarm.LacPredict(i) = mean(SUM2.LacPredict(ind));
    sfarm.TotalResid(i) = mean(SUM2.TotalResid(ind));
    sfarm.RelativeTotal(i) = mean(SUM2.RelativeTotal(ind));
    sfarm.RelativeDay(i) = nanmean(SUM2.RelativeDay(ind));
    % lac = 1
    sfarm.NoHeif(i,1) = length(ind1);
    sfarm.NoMeasHeif(i,1) = length(ind31);
    sfarm.AvSCCHeif(i,1) = nanmean(MPRi.SCC(ind31));
    sfarm.AverageNoMeasHeif(i) = mean(SUM2.NoMeas(ind1));
    sfarm.AvPercHighHeif(i) = mean(SUM2.No150(ind1)./SUM2.NoMeas(ind1))*100;
    sfarm.LacProdHeif(i) = mean(SUM2.LacProd(ind1));
    sfarm.LacPredictHeif(i) = mean(SUM2.LacPredict(ind1));
    sfarm.TotalResidHeif(i) = mean(SUM2.TotalResid(ind1));
    sfarm.RelativeTotalHeif(i) = mean(SUM2.RelativeTotal(ind1));
    sfarm.RelativeDayHeif(i) = nanmean(SUM2.RelativeDay(ind1));
    % lac > 1
    sfarm.NoCowsP2(i,1) = length(ind2);
    sfarm.NoMeasCow(i,1) = length(ind32);
    sfarm.AvSCCCow(i,1) = nanmean(MPRi.SCC(ind32));
    sfarm.AverageNoMeasCow(i) = mean(SUM2.NoMeas(ind2));
    sfarm.AvPercHighCow(i) = mean(SUM2.No250(ind2)./SUM2.NoMeas(ind2))*100;
    sfarm.LacProdCow(i) = mean(SUM2.LacProd(ind2));
    sfarm.LacPredictCow(i) = mean(SUM2.LacPredict(ind2));
    sfarm.TotalResidCow(i) = mean(SUM2.TotalResid(ind2));
    sfarm.RelativeTotalCow(i) = mean(SUM2.RelativeTotal(ind2));
    sfarm.RelativeDayCow(i) = nanmean(SUM2.RelativeDay(ind2));
end
    

% plot correlations at farm level
figure('Units','normalized','OuterPosition',[-0.05 0.4 1.1 0.40]); 
for i = 1:5
    x = sfarm.AvSCCall;
    y = sfarm{:,8+i};
%     b = regress(y,[ones(length(x),1) x]);
    D = sprintf('all%i',i);
    [b,stats.(D)] = robustfit(x,y);
    p(1,i) = stats.(D).p(2);
    subplot(1,5,i); hold on; box on; xlabel('SCC [cells/mL]')
    plot(x,y,'o','LineWidth',1.5,'Color',[0.2 0.2 0.8],'MarkerSize',3)
    plot((0:450),b(1)+b(2)*(0:450),'-','LineWidth',2,'Color',[0.2 0.2 0.8])
    x = sfarm.AvSCCHeif;
    y = sfarm{:,18+i};
%     b = regress(y,[ones(length(x),1) x]);
    D = sprintf('heif%i',i);
    [b,stats.(D)] = robustfit(x,y);
    p(2,i) = stats.(D).p(2);
    subplot(1,5,i); hold on; box on; xlabel('SCC [cells/mL]')
    plot(x,y,'o','LineWidth',1.5,'Color',[0.2 0.6 0.2],'MarkerSize',3)
    plot((0:450),b(1)+b(2)*(0:450),'-','LineWidth',2,'Color',[0.2 0.6 0.2])
    x = sfarm.AvSCCCow;
    y = sfarm{:,28+i};
    D = sprintf('cow%i',i);
    [b,stats.(D)] = robustfit(x,y);
    p(3,i) = stats.(D).p(2);
    subplot(1,5,i); hold on; box on; xlabel('SCC [cells/mL]')
    plot(x,y,'o','LineWidth',1.5,'Color',[0.8 0.0 0.2],'MarkerSize',3)
    plot(0:450,b(1)+b(2)*(0:450),'-','LineWidth',2,'Color',[0.8 0.0 0.2])
    axis tight; xlim([50 450])
end

subplot(1,5,1); title('Average true production'); ylabel('MY [kg]')
subplot(1,5,2); title('Average expected production'); ylabel('MY [kg]')
subplot(1,5,3); title('Average loss per lactation'); ylabel('Loss/lac [kg]')
subplot(1,5,4); title('Relative loss per lactation'); ylabel('Loss/lac [%/lac]')
subplot(1,5,5); title('Relative loss per day'); ylabel('Loss/day [%/d]')

clear ans b D data h i ind ind1 ind2 ind3 ind31 ind32 j m p x xticks y ypos ks

% at individual level, ignore farms
figure('Units','normalized','OuterPosition',[0.05 0.2 0.7 0.80]);
for i = 1:3
    % all
    subplot(2,3,i); hold on; xlim([2.2 9.2]); box on; 
    x = log(SUM2.MaxSCC);
    y = SUM2{:,15+i};
    plot(x,y,'o','MarkerSize',3,'Color',[0.2 0.2 0.8])
    D = sprintf('all%i',i);
    [b.(D),stats.(D)] = robustfit(x,y);
    plot(0:10,b.(D)(1)+b.(D)(2)*(0:10),'-','LineWidth',2,'Color',[0 0.8 1])
    p(1,i) = stats.(D).p(2);
    % heif
    x = log(SUM2.MaxSCC(SUM2.Lac==1));
    y = SUM2{SUM2.Lac==1,15+i};
    lh = plot(x,y,'o','MarkerSize',3,'Color',[0.2 0.6 0.2]); lh.Color = [0.2 0.6 0.2 0.2];
    D = sprintf('heif%i',i);
    [b.(D),stats.(D)] = robustfit(x,y);
    plot(0:10,b.(D)(1)+b.(D)(2)*(0:10),'-','LineWidth',2,'Color',[0.4 0.8 0.6])
    p(2,i) = stats.(D).p(2);
    % cows
    x = log(SUM2.MaxSCC(SUM2.Lac>1));
    y = SUM2{SUM2.Lac>1,15+i};
    lh = plot(x,y,'o','MarkerSize',3,'Color',[0.8 0 0.2]); lh.Color = [0.8 0 0.2 0.8];
    D = sprintf('cow%i',i);
    [b.(D),stats.(D)] = robustfit(x,y);
    plot(0:10,b.(D)(1)+b.(D)(2)*(0:10),'-','LineWidth',2,'Color',[1 0.4 0.4])
    p(3,i) = stats.(D).p(2);
   
end
subplot(2,3,1); ylim([-2000 100]); title('Max SCC vs. milk loss/lac'); xlabel('log(SCC)'); ylabel('Milk loss (kg)')
subplot(2,3,2); ylim([-20 2]); title('Max SCC vs. milk loss/lac %'); xlabel('log(SCC)'); ylabel('Milk loss (%/lac)')
subplot(2,3,3); ylim([-20 2]); title('Max SCC vs. milk loss/day'); xlabel('log(SCC)'); ylabel('Milk loss (%/d)')

for i = 1:3
    subplot(2,3,3+i); hold on; box on

    D = sprintf('all%i',i);
    plot(exp(3:10)-exp(2:9),b.(D)(1)+b.(D)(2)*(2:9),'LineWidth',2,'Color',[0 0.8 1])

    D = sprintf('heif%i',i);
    plot(exp(3:10)-exp(2:9),b.(D)(1)+b.(D)(2)*(2:9),'LineWidth',2,'Color',[0.4 0.8 0.6])

    D = sprintf('cow%i',i);
    plot(exp(3:10)-exp(2:9),b.(D)(1)+b.(D)(2)*(2:9),'LineWidth',2,'Color',[1 0.4 0.4])
    grid on; xlim([0 10000])
    legend({'all','Lac=1','Lac>1'})
end
subplot(2,3,4); ylim([-325 -180]); title('Max SCC vs. milk loss/lac'); xlabel('SCC [cells/mL]'); ylabel('Milk loss (kg)')
subplot(2,3,5); ylim([-3.25 -2]); title('Max SCC vs. milk loss/lac %'); xlabel('SCC [cells/mL]'); ylabel('Milk loss (%/lac)')
subplot(2,3,6); ylim([-3.25 -2]); title('Max SCC vs. milk loss/day'); xlabel('SCC [cells/mL]'); ylabel('Milk loss (%/d)')


clear b ans D i ks lh p pl x y p 


% milk losses surrounding each MPR measurement
tic
for i = 1:height(MPRi)
    farm = t.Var2(contains(t.farms,MPRi.FarmName{i}));
    
    ind = find(contains(moddataDAY.(farm{1}).OffRegNo,MPRi.OffRegNo{i}) == 1 & ...
               moddataDAY.(farm{1}).Lac == MPRi.Lac(i) & ...
               moddataDAY.(farm{1}).DIM >= MPRi.DIM(i)-15 & ...
               moddataDAY.(farm{1}).DIM <= MPRi.DIM(i)+15);
             
    data = moddataDAY.(farm{1})(ind,:);
    ind = find(data.DIM == MPRi.DIM(i)); % 
    
    data.idx = data.DIM - MPRi.DIM(i)+16;
    
    
    data.Res(:,1) = data.MYs-data.ModIT;
    data.RelRes(:,1) = data.Res(:,1)./data.ModIT(:,1)*100;
    
    MPRi.NoNegRes(i) = length(find(data.Res < 0));
    MPRi.SumRes(i) = sum(data.Res); 
    MPRi.AvRelRes(i) = mean(data.RelRes);
    MPRi.RelRes(i) = sum(data.MYs)/sum(data.ModIT)*100;
    MPRi.Loss(i,data.idx) = data.Res';
    MPRi.RelLoss(i,data.idx) = data.RelRes';
end
toc
save('D11_MPRloss.mat','MPRi')

% delete MPRi
MPRi(MPRi.DIM>305,:) = [];   % 57245



% 
for i = 1:31
    MPRi.Loss(MPRi.Loss(:,i) == 0,i) = NaN;
    MPRi.RelLoss(MPRi.RelLoss(:,i) == 0,i) = NaN;
end

% Sumamrize losses
th = [0 100 200 300 500 700 900 1500 10000];
for i = 1:length(th)-1
    
    ind = find(MPRi.SCC > th(i) & MPRi.SCC <= th(i+1));
    for j = 1:31
        medians(i,j) = nanmedian(MPRi.Loss(ind,j));
        quart25(i,j) = quantile(MPRi.Loss(ind,j),0.25);
        quart75(i,j) = quantile(MPRi.Loss(ind,j),0.75);
        
        mediansREL(i,j) = nanmedian(MPRi.RelLoss(ind,j));
        quart25REL(i,j) = quantile(MPRi.RelLoss(ind,j),0.25);
        quart75REL(i,j) = quantile(MPRi.RelLoss(ind,j),0.75);
    end
end

% visualize 'losses'
figure('Units','normalized','OuterPosition',[0.1 0.1 0.6 0.5]);hold on; box on;
col = [zeros(8,1) (1:-1/7:0)' 0.6*ones(8,1)];
ylim([-5.5 1.5]); 
for i = 1:8
    hold on
    plot(-15:15, medians(i,:),'LineWidth',3.5,'Color',col(i,:))
end
plot([0 0],[-5.5 1.5],'-','Color',[0.8 0 0.2],'LineWidth',3.5)
plot([-15 15],[0 0],'k-.','LineWidth',1)
legend({'SCC<100','SCC<200','SCC<300','SCC<500','SCC<700','SCC<900','SCC<1500','SCC>1500'},'Location','Southwest','AutoUpdate','off')
for i = 1:8
    hold on
    plot(-15:15, quart25(i,:),':','LineWidth',1,'Color',col(i,:))
    plot(-15:15, quart75(i,:),':','LineWidth',1,'Color',col(i,:))
end
xlabel('Days from SCC measurement')
ylabel('Milk loss [kg]')
title('Milk losses around SCC measurements')

% visualize 'losses'
figure('Units','normalized','OuterPosition',[0.1 0.1 0.6 0.5]);hold on; box on;
col = [zeros(8,1) (1:-1/7:0)' 0.6*ones(8,1)];
ylim([-16.5 4]); 
for i = 1:8
    hold on
    plot(-15:15, mediansREL(i,:),'LineWidth',3.5,'Color',col(i,:))
end
plot([0 0],[-16.5 4],'-','Color',[0.8 0 0.2],'LineWidth',3.5)
plot([-15 15],[0 0],'k-.','LineWidth',1)
legend({'SCC<100','SCC<200','SCC<300','SCC<500','SCC<700','SCC<900','SCC<1500','SCC>1500'},'Location','Southwest','AutoUpdate','off')
for i = 1:8
    hold on
    plot(-15:15, quart25REL(i,:),':','LineWidth',1,'Color',col(i,:))
    plot(-15:15, quart75REL(i,:),':','LineWidth',1,'Color',col(i,:))
end
xlabel('Days from SCC measurement')
ylabel('Milk loss [%/d]')
title('Milk losses around SCC measurements')


%
figure; hold on
for i = 1:length(th)-1
    ind = find(MPRi.SCC > th(i) & MPRi.SCC <= th(i+1));

    histogram(MPRi.RelLoss(ind,:),'Normalization','probability')
    xlim([-100 35])
end



% clear
clear ans quart25 quart25REL quart75 quart75REL th i ind j mean medians mediansREL col
clear data means farms farms2 farm fn


% example figure
ans = find(MPRi.SCC > 4000);
ind = find(contains(moddataDAY.Almey.OffRegNo,'BE 410970530') & moddataDAY.Almey.Lac == 6);

figure;
yyaxis left; hold on;
plot(moddataDAY.Almey.DIM(ind),moddataDAY.Almey.MYs(ind),'o-','MarkerSize',3)
plot(moddataDAY.Almey.DIM(ind),moddataDAY.Almey.ModIT(ind),'-','LineWidth',2.5)

ind = find(contains(MPRi.OffRegNo,'BE 410970530') & MPRi.Lac == 6);
yyaxis right; plot(MPRi.DIM(ind),MPRi.SCC(ind),'s','LineWidth',1.5)
xlim([0 250])
ylabel('SCC (*1000 cells/mL)')
xlabel('DIM (days)')
%% summarize daily and quarter per farm
length(fieldnames(moddata))
length(fieldnames(moddataDAY))

fn  = fieldnames(moddata);
for i = 1:length(fieldnames(moddata))
%     sumday(i,1) = height(moddataDAY.(fn{i})); % no meas
%     sumday(i,2) = length(unique(moddataDAY.(fn{i}).B_ID));
%     sumday(i,3) = height(unique(moddataDAY.(fn{i})(:,[2 8]),'rows'));
%     sumday(i,4) = max(datenum(moddataDAY.(fn{i}).Date))-min(datenum(moddataDAY.(fn{i}).Date));
    
%     summilk(i,1) = height(moddata.(fn{i})); % no meas
    summilk(i,2) = length(unique(moddata.(fn{i}).B_ID));
    summilk(i,3) = height(unique(moddata.(fn{i})(:,[2 6]),'rows'));
    summilk(i,4) = max(datenum(moddata.(fn{i}).Date))-min(datenum(moddata.(fn{i}).Date));

end

sum(sumday)
sum(summilk)




%% summarize pert quarter elvel losses

figure; subplot(1,2,1); hold on; box on
h = histogram(-sumary.MLaF);
    x = h.BinEdges;
    y = h.Values;
    ypos = h.Values;ind = find(y == 0); y(ind) = []; x(ind) = [];
    text(x(1:end-1)+4.5,y,num2str(y'),'vert','bottom','horiz','center')

h.EdgeColor = [0 0.8 0.6]; h.FaceColor = [0 0.8 0.6]; h.LineWidth = 2;
h = histogram(-sumary.MLaH);
    x = h.BinEdges;
    y = h.Values;
    ypos = h.Values;ind = find(y == 0); y(ind) = []; x(ind) = [];
    text(x(1:end-1)+4.5,y,num2str(y'),'vert','bottom','horiz','center')
h.EdgeColor = [0 0.4 0.6]; h.FaceColor = [0 0.4 0.6]; h.LineWidth = 2;
xlabel('Absolute milk loss (kg)'); ylabel('No. of farms'); ylim([0 21])
legend({'Front Q','Hind Q'},'Location','northwest')
title('Absolute milk loss per perturbation')


subplot(1,2,2); hold on; box on
binedges = -24:2:0;
h = histogram(-sumary.MLrF,binedges);
    x = h.BinEdges;
    y = h.Values;
    ypos = h.Values;ind = find(y == 0); y(ind) = []; x(ind) = [];
    text(x(1:end-1)+0.9,y,num2str(y'),'vert','bottom','horiz','center')
h.EdgeColor = [0 0.8 0.6]; h.FaceColor = [0 0.8 0.6]; h.LineWidth = 2;
h = histogram(-sumary.MLrH,binedges);
    x = h.BinEdges;
    y = h.Values;
    ypos = h.Values;ind = find(y == 0); y(ind) = []; x(ind) = [];
    text(x(1:end-1)+0.9,y,num2str(y'),'vert','bottom','horiz','center')
h.EdgeColor = [0 0.4 0.6]; h.FaceColor = [0 0.4 0.6]; h.LineWidth = 2;
xlabel('Relative milk loss (%)'); ylabel('No. of farms'); ylim([0 18])
legend({'Front Q','Hind Q'},'Location','northeast')
title('Relative milk loss per milking')





























































