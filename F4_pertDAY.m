function [PT_zero,h] = F4_pertDAY(DIM,MY,MOD,minProcLoss,minLengthNeg,makefig)
% this function returns the different perturbations of a quarter lactation,
% based on the DIM and the residuals of a model given by RES.
% INPUTS:
%       1   DIM days in milk of measurements
%       2   QMY = matrix of 4 cols, each = QMY of 1 quarter
%       3   MOD = matrix of 4 cols, each = model of QMY of 1 quarter
%       4   minProcLoss = threshold minimal loss
%       5   minLengthNeg = minimal length below zero (neg pert) = DIM
%       6   makefig = whether or not to make a figure
% OUTPUTS:
%       1   PT_zero  is the overview of overlapping perturbations
%       2   h       is the figure handle if makefig is 1
% 
% STEPS:
%   STEP 1:  calculate the residuals of TMY against MOD 
%   STEP 2:  express residuals in percentages
%   STEP 3:  smooth with median smoother
%   STEP 4:  identify indices of smoother below zero
%   STEP 5:  identify indices of smoother below minProcLoss
%   STEP 6:  detect start and end dates of the perturbations
%   STEP 7:  detect overlap and merge to single perturbations
%   STEP 8:  detect perturbations that comply with both criteria
%   STEP 9:  assign perturbation numbers + construct 'overlap' dataset 
%   STEP 10: prepare outputs



% ================= % for development purposes only % ================= %
minProcLoss = 80;   % at least 20% loss
minLengthProc = 5;  % minimum 5 days negative 
minLengthNeg = 5;
data = moddataDAY.Vanhullebus;
i = 263;
cowlac = unique(data(:,[2 8]),'rows');
DIM = data.DIM(data.B_ID == cowlac.B_ID(i) & data.Lac == cowlac.Lac(i));
MY = data{data.B_ID == cowlac.B_ID(i) & data.Lac == cowlac.Lac(i),13};
MOD = data{data.B_ID == cowlac.B_ID(i) & data.Lac == cowlac.Lac(i),14};
makefig = 1;
% ================= % for development purposes only % ================= %

warning('off','all')

%% FUNCTION

% create NUM = order of measurements
NUM = (1:length(DIM))';

if makefig == 1
    % prepare plot
    h = figure('Units','normalized','Outerposition',[0 0.5 1 0.5]); 
    subplot(2,1,1); hold on; axis tight; box on; ylim([0 max(MY)+5]);xlabel('DIM (days)'); ylabel('Daily MY (kg)')
    plot(DIM, MY, 'o:','LineWidth',1.2,'MarkerSize',3,'Color',[0.8 0.6 0.6])
    plot(DIM,MOD,'-','LineWidth',2,'Color',[0.2 0 0.6])    
    subplot(2,1,2); hold on; axis tight; box on; xlabel('DIM (days)'); ylabel('Residual MY (%)')
    plot(DIM,zeros(length(DIM),1),'k--','LineWidth',1.5)
    plot(DIM,-(100-minProcLoss)*ones(length(DIM),1),':','LineWidth',1.5,'Color',[0 0.4 0.8])
else
    h = 0;
end

% calculate residuals
RES = MY-MOD;

% express in percentages
RES = (RES./MOD)*100;

if makefig == 1
    subplot(2,1,2); plot(DIM,RES,'o:','LineWidth',1.2,'MarkerSize',3,'Color',[0.8 0.6 0.6])
    ylim([-100 100])
end


% find measurements below zero
ind = [];
ind = find(RES < 0);

% find measurements below minProcLoss
idx = [];
idx = find(RES < -(100-minProcLoss));

% select measurements below 0
ptb = [];
ptb = [DIM(ind) NUM(ind) RES(ind)];

if makefig == 1
    subplot(2,1,2);
    plot(ptb(:,1),ptb(:,3),'kx','MarkerSize',4,'LineWidth',1)
end

% select measurements below threshold
ptb2 = [];
ptb2 = [DIM(idx) NUM(idx) RES(idx)];


if makefig == 1
    subplot(2,1,2);
    plot(ptb2(:,1),ptb2(:,3),'rx','MarkerSize',4,'LineWidth',1)
end


% start and end of each perturbation, DIM and index
PT_zero = array2table(ind, 'VariableNames',{'Index'});
for i = 1:size(ind,1)
    idx1 = find(RES > 0 & DIM < DIM(ind(i)), 1, 'last');% index of last above 0
    idx1 = idx1+1;                                      % index of first below 0
    if isempty(idx1) ==1; idx1 = 1; end                 % no index then is first
    
    PT_zero.IndexStart(i,1) = idx1;                     % fill in index
    PT_zero.DIMStart(i,1) = DIM(idx1);                  % fill in DIM
    
    idx2 = find(RES > 0 & DIM > DIM(ind(i)), 1,'first');% index first above 0
    idx2 = min(idx2-1,length(DIM));                     % index of last below 0
    if isempty(idx2) ==1; idx2 = length(DIM); end       % no index then is last
    PT_zero.IndexEnd(i,1) = idx2;                       % fill in index
    PT_zero.DIMEnd(i,1) = DIM(idx2);                    % fill in DIM
end

% select unique perturbations
if ~isempty(PT_zero)
    [~,uniPidx] = unique(PT_zero(:,2:5),'rows');
    PT_zero = PT_zero(uniPidx,:);
    PT_zero = sortrows(PT_zero,1);

    % calculate length and only keep > minLengthNeg
    PT_zero.PertLength(:,1) = PT_zero.DIMEnd - PT_zero.DIMStart+1;
    PT_zero(PT_zero.PertLength < minLengthNeg,:) = [];

    if makefig == 1
        for i = 1:height(PT_zero)
            subplot(2,1,2)
            plot([PT_zero.DIMStart(i) PT_zero.DIMStart(i)],[-100 100],'-.','LineWidth',1.6,'Color',[0 0.6 00])
            plot([PT_zero.DIMEnd(i) PT_zero.DIMEnd(i)],[-100 100],'-.','LineWidth',1.6,'Color',[0.8 0 0.2])
        end
    end
else
    PT_zero = array2table(zeros(0,10),'VariableNames',{'NoPert','threshold',...
                        'IndexStart','DIMStart','IndexEnd','DIMEnd','PertLength',...
                        'MLa','MLr','MLl'});
end

% find whether there is at least one measurement below the threshold
for i = 1:height(PT_zero)
    idx = find(DIM >= PT_zero.DIMStart(i) & ...
               DIM <= PT_zero.DIMEnd(i) & ...
               RES <= -(100-minProcLoss)); % at least ONCE below threshold
    
    PT_zero.threshold(i,1) = length(idx); 
end

try
    PT_zero(PT_zero.threshold == 0,:) = [];
    PT_zero = removevars(PT_zero,{'threshold','Index'});
catch
    PT_zero = removevars(PT_zero,{'Index'});    
end

% add number of perturbaions or prepare table output without
if isempty(PT_zero)
    % prepare PT_zero
    PT_zero = array2table(zeros(0,9),'VariableNames',{'NoPert',...
                        'IndexStart','DIMStart','IndexEnd','DIMEnd','PertLength',...
                        'MLa','MLr','MLl'});
else
    PT_zero.NoPert(:,1) = 1:height(PT_zero);
    PT_zero = movevars(PT_zero,'NoPert','Before','IndexStart');
end

if ~isempty(PT_zero) && makefig == 1
    subplot(2,1,2); 
    for i = 1:height(PT_zero)
        plot([PT_zero.DIMStart(i) PT_zero.DIMEnd(i)] ,[45 45],'LineWidth',6,'Color',[0.6 0 0.6])
    end
end




% clear variables
clear A ans cowsel i fileName idx idx1 idx2 ind j minLengthNeg 
clear minLengthProc minProcLoss NUM opts ptb ptb2 ptb_sel T uniPidx yesno

% add milk losses to PT_zero
if ~isempty(PT_zero)
    PT_zero.MLa(:,1) = NaN;       % prepare absolute loss Q1
    PT_zero.MLr(:,1) = NaN;       % prepare relative loss Q1
    PT_zero.MLl(:,1) = NaN;       % prepare lactation loss Q1

    for i = 1:height(PT_zero) % all perturbations
        idx = [0 0];
        idx(1) = PT_zero.IndexStart(i);     % idx of DIM start
        idx(2) = PT_zero.IndexEnd(i);      % idx of DIM end
        PT_zero.MLa(i) = -sum(MY(idx(1):idx(2))-MOD(idx(1):idx(2))); % absolute losses = total produced - total expected in kg
        PT_zero.MLr(i) = 100-(sum(MY(idx(1):idx(2)))/sum(MOD(idx(1):idx(2)))*100); % relative losses = 100%- total produced / total expected
        PT_zero.MLl(i) = (PT_zero.MLa(i)/sum(MOD))*100; % (total produced-total expected)/(lactation expected)
    end
end

% put warnings on again
warning('on','all')