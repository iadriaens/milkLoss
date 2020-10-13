function [PT_zero,PT_perc,PT_all,h] = F2_pertQuarter(DIM,QMY,MOD,minProcLoss,minLengthNeg,minLengthProc,makefig)
% this function returns the different perturbations of a quarter lactation,
% based on the DIM and the residuals of a model given by RES.
% INPUTS:
%       1   DIM days in milk of measurements
%       2   QMY = matrix of 4 cols, each = QMY of 1 quarter
%       3   MOD = matrix of 4 cols, each = model of QMY of 1 quarter
%       4   minProcLoss = threshold minimal loss
%       5   minLengthNeg = minimal length below zero (neg pert) = DIM
%       6   minLengthProc = minimal length below thresholds = MEASUREMENTS
% OUTPUTS:
%       1   PT_zero is an overview of the pertubations below zero longer 
%                   than minLengthNeg in DIM
%       2   PT_perc is an overview of the perturbation below minProcLoss
%                   and longer than minLengthProc in MEASUREMENTS
%       3   PT_all  is the overview of overlapping perturbations at QUARTER 
%       3           level
%       4   h       is the figure handle if makefig is 1
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
% % minProcLoss = 80;   % at least 20% loss
% % minLengthNeg = 10;  % minimum 10 days negative 
% % minLengthProc = 5;  % minimum 5 days below threshold
% % % % 
% % data = moddata.Wildemauwe;
% % i = 2;
% % cowlac = unique(data(:,[2 6]),'rows');
% % DIM = data.DIM(data.B_ID == cowlac.B_ID(i) & data.Lac == cowlac.Lac(i));
% % QMY = data{data.B_ID == cowlac.B_ID(i) & data.Lac == cowlac.Lac(i),[13:16]};
% % MOD = data{data.B_ID == cowlac.B_ID(i) & data.Lac == cowlac.Lac(i),[22:25]};
% % makefig = 1;
% ================= % for development purposes only % ================= %

warning('off','all')

%% FUNCTION

% create NUM = order of measurements
NUM = (1:length(DIM))';

if makefig == 1
    % prepare plot
    h = figure('Units','normalized','Outerposition',[0 0 1 1]); 
    for j =1:4
        subplot(2,2,j); hold on; axis tight; box on; xlabel('DIM (days)'); ylabel('Residual MY (%)')
    end
else
    h = 0;
end

% calculate residuals
RES = QMY-MOD;

% express in percentages
RES = (RES./MOD)*100;

if makefig == 1
    for j = 1:4
        subplot(2,2,j); plot(DIM,RES(:,j),'o:','LineWidth',1.2,'MarkerSize',3,'Color',[0.8 0.6 0.6])
        ylim([-100 100])
    end
end

% smooth with median smoother, window of 7 measurements
RES = medfilt1(RES,7,[],'omitnan','truncate'); % smoothed res

% alternative for calculations
RES2 = medfilt1(QMY-MOD,7,[],'omitnan','truncate');

if makefig == 1
    for j = 1:4
        subplot(2,2,j)
        plot(DIM,RES(:,j),'-','LineWidth',1.4,'Color',[0 0 0.6])
        plot(DIM,-(100-minProcLoss)*ones(length(DIM),1),'-','LineWidth',1.7,'Color',[0 0.6 0.6])
        plot(DIM,zeros(length(DIM),1),':','LineWidth',1.2,'Color',[0.6 0 0])
    end
end

% find measurements below zero
ind = [];
[ind(:,1),ind(:,2)] = find(RES < 0);
ind = sortrows(ind,2);

% find measurements below minProcLoss
idx = [];
[idx(:,1),idx(:,2)] = find(RES < -(100-minProcLoss));
idx = sortrows(idx,2);

% select measurements below 0
ptb = [];
for j = 1:4
    ptb = [ptb; DIM(ind(ind(:,2)==j,1)) NUM(ind(ind(:,2)==j,1),1) RES(ind(ind(:,2)==j,1),j) j*ones(length(find(ind(:,2)==j)),1)];
end

if makefig == 1
    for j = 1:4
        subplot(2,2,j);
        plot(ptb(ptb(:,4) == j,1),ptb(ptb(:,4) == j,3),'kx','MarkerSize',3,'LineWidth',1)
    end
end

% select measurements below threshold
ptb2 = [];
for j = 1:4
    ptb2 = [ptb2; DIM(idx(idx(:,2)==j,1)) NUM(idx(idx(:,2)==j,1),1) RES(idx(idx(:,2)==j,1),j) j*ones(length(find(idx(:,2)==j)),1)];
end

if makefig == 1
    for j = 1:4
        subplot(2,2,j);
        plot(ptb2(ptb2(:,4) == j,1),ptb2(ptb2(:,4) == j,3),'rx','MarkerSize',3,'LineWidth',1)
    end
end

% start and end of each perturbation, DIM and index
PT_zero = array2table(ind, 'VariableNames',{'Index','Quarter'});
for i = 1:size(ind,1)
    idx1 = find(RES(:,ind(i,2)) > 0 & DIM < DIM(ind(i,1)), 1, 'last');% index of last above 0
    idx1 = idx1+1;                                      % index of first below 0
    if isempty(idx1) ==1; idx1 = 1; end                 % no index then is first
    
    PT_zero.IndexStart(i,1) = idx1;                     % fill in index
    PT_zero.DIMStart(i,1) = DIM(idx1);                  % fill in DIM
    
    idx2 = find(RES(:,ind(i,2)) > 0 & DIM > DIM(ind(i,1)), 1,'first');% index first above 0
    idx2 = min(idx2-1,length(DIM));                     % index of last below 0
    if isempty(idx2) ==1; idx2 = length(DIM); end       % no index then is last
    PT_zero.IndexEnd(i,1) = idx2;                       % fill in index
    PT_zero.DIMEnd(i,1) = DIM(idx2);                    % fill in DIM
end

% select unique perturbations
if ~isempty(PT_zero)
    [~,uniPidx] = unique(PT_zero(:,[2 3 4 5 6]),'rows');
    PT_zero = PT_zero(uniPidx,:);
    PT_zero = sortrows(PT_zero,[2 1]);

    % calculate length and only keep > minLengthNeg
    PT_zero.PertLength(:,1) = PT_zero.DIMEnd - PT_zero.DIMStart;
    PT_zero(PT_zero.PertLength < minLengthNeg,:) = [];

    if makefig == 1
        for i = 1:height(PT_zero)
            subplot(2,2,PT_zero.Quarter(i))
            plot([PT_zero.DIMStart(i) PT_zero.DIMStart(i)],[-100 100],'-.','LineWidth',1.6,'Color',[0 0.6 00])
            plot([PT_zero.DIMEnd(i) PT_zero.DIMEnd(i)],[-100 100],'-.','LineWidth',1.6,'Color',[0.8 0 0.2])
        end
    end
end

% detect / set minimum of the perturbation
if ~isempty(PT_zero)
    for k = 1:height(PT_zero)
        [~,sortedindx] = sortrows(RES(PT_zero.IndexStart(k):PT_zero.IndexEnd(k),PT_zero.Quarter(k)));
        sst = mean(DIM(PT_zero.IndexStart(k)+sortedindx(1:3)-1));
        win = 0.5*(PT_zero.DIMEnd(k)-PT_zero.DIMStart(k)); % window 30%
        PT_zero.DIMwinMin(k) = max(sst-win,PT_zero.DIMStart(k));       % minimum of window for overlap
        PT_zero.DIMwinMax(k) = min(sst+win,PT_zero.DIMEnd(k));       % max window for overlap
        PT_zero.DIMwin(k) = sst;            % minimum
    end
end


% start and end of each perturbation, DIM and index percent
if ~isempty(idx)
    PT_perc = array2table(idx, 'VariableNames',{'Index','Quarter'});
    for i = 1:size(idx,1)
        idx1 = find(RES(:,idx(i,2)) > -(100-minProcLoss) & DIM < DIM(idx(i,1)), 1, 'last');% index of last above 0
        idx1 = idx1+1;                                      % index of first below 0
        if isempty(idx1) ==1; idx1 = 1; end                 % no index then is first

        PT_perc.IndexStart(i,1) = idx1;                     % fill in index
        PT_perc.DIMStart(i,1) = DIM(idx1);                  % fill in DIM

        idx2 = find(RES(:,idx(i,2)) > -(100-minProcLoss) & DIM > DIM(idx(i,1)), 1,'first');% index first above 0
        idx2 = min(idx2-1,length(DIM));                     % index of last below 0
        if isempty(idx2) ==1; idx2 = length(DIM); end       % no index then is last
        PT_perc.IndexEnd(i,1) = idx2;                       % fill in index
        PT_perc.DIMEnd(i,1) = DIM(idx2);                    % fill in DIM
    end
else
    PT_perc = array2table(zeros(0,11),'VariableNames',{'NoPert','Quarter',...
                        'IndexStart','DIMStart','IndexEnd','DIMEnd','PertLength','OverlapCritDIM',...
                        'MLa','MLr','MLl'}); 
end

% select unique perturbations
if ~isempty(PT_perc)
    [~,uniPidx] = unique(PT_perc(:,[2 3 4 5 6]),'rows');
    PT_perc = PT_perc(uniPidx,:);
    PT_perc = sortrows(PT_perc,[2 1]);

    % calculate length and only keep > minMeasNeg
    PT_perc.PertLength(:,1) = PT_perc.IndexEnd - PT_perc.IndexStart;
    PT_perc(PT_perc.PertLength < minLengthProc,:) = [];
    
    if makefig == 1
        for i = 1:height(PT_perc)
            subplot(2,2,PT_perc.Quarter(i))
            plot([PT_perc.DIMStart(i) PT_perc.DIMEnd(i)],[45 45],'-','LineWidth',6,'Color',[1 0.4 0])
        end
    end
end

% if all deleted
if isempty(PT_perc)
    PT_perc = array2table(zeros(0,14),'VariableNames',{'NoPert','Quarter',...
              'IndexStart','DIMStart','IndexEnd','DIMEnd','PertLength','OverlapCritDIM',...
              'MLa','MLr','MLl','smMLa','smMLr','smMLl'}); 
end


% overlap between PT_zero & PT_perc
if ~isempty(PT_perc) && ~isempty(PT_zero)
    for i = 1:height(PT_zero)
        yesno = find(PT_perc.Quarter == PT_zero.Quarter(i) & ...
                     PT_perc.IndexStart >= PT_zero.IndexStart(i) &...
                     PT_perc.IndexEnd <= PT_zero.IndexEnd(i));
        PT_zero.yesno(i) = length(yesno);
    end
    PT_zero(PT_zero.yesno == 0,:) = [];
end

if ~isempty(PT_perc) && ~isempty(PT_zero)
    PT_zero.Properties.VariableNames{11} = 'NoNegEp';
elseif isempty(PT_perc) && ~isempty(PT_zero)
    PT_zero.NoNegEp(:,1) = 0;
end
    

if ~isempty(PT_zero)
    % select / remove if not both criteria
    PT_zero = removevars(PT_zero,{'Index'});
end

% arrange in table that combines overlapping perturbations of different
% quarters
PT_all = array2table(zeros(0,9),'VariableNames',{'PertNo','DIMstart1','DIMend1',...
                            'DIMstart2','DIMend2','DIMstart3','DIMend3',...
                            'DIMstart4','DIMend4'});
if ~isempty(PT_zero) && ~isempty(PT_perc)
    T = 0;
    PT_del = sortrows(PT_zero,'DIMStart');
    while ~isempty(PT_del)
        T =T+1;
        PT_all.PertNo(T) = T;
        PT_all{T,2:end} = NaN;
        
        % using just overlap
% % %         ind = find(~(PT_del.IndexStart > PT_del.IndexEnd(1)) & ... % find overlap, this means
% % %             ~(PT_del.IndexEnd < PT_del.IndexStart(1))); % find end not before start and start not after end
        
        % using the window defined
        ind = find(PT_del.DIMwin >= PT_del.DIMwinMin(1) & ...
                   PT_del.DIMwin <= PT_del.DIMwinMax(1));

        % fill in
        if length(ind)==1  % if no overlapping perturbations in other quarters
            PT_all{T,1+2*PT_del.Quarter(ind)-1} = PT_del.DIMStart(ind);
            PT_all{T,1+2*PT_del.Quarter(ind)} = PT_del.DIMEnd(ind);
        else
            for j = 1:length(ind)
                PT_all{T,1+2*PT_del.Quarter(ind(j))-1} = PT_del.DIMStart(ind(j));
                PT_all{T,1+2*PT_del.Quarter(ind(j))} = PT_del.DIMEnd(ind(j));
            end
        end
        PT_del(ind,:) = [];
    end
else
    PT_all = array2table(zeros(0,33),'VariableNames',{'PertNo','DIMstart1',...
                        'DIMend1','DIMstart2','DIMend2','DIMstart3','DIMend3','DIMstart4','DIMend4',...
                        'ML1a','ML2a','ML3a','ML4a','ML1r','ML2r','ML3r','ML4r','ML1l','ML2l','ML3l','ML4l',...
                        'smML1a','smML2a','smML3a','smML4a','smML1r','smML2r','smML3r','smML4r',...
                        'smML1l','smML2l','smML3l','smML4l'});
end
clear PT_del



if isempty(PT_zero)
    % prepare PT_zero
    PT_zero = array2table(zeros(0,14),'VariableNames',{'NoPert','Quarter',...
                        'IndexStart','DIMStart','IndexEnd','DIMEnd','PertLength','NoNegEp',...
                        'MLa','MLr','MLl','smMLa','smMLr','smMLl'});
else
    PT_zero = removevars(PT_zero,{'DIMwinMin','DIMwinMax','DIMwin'});
end

% clear variables
clear A ans cowsel i fileName idx idx1 idx2 ind j minLengthNeg 
clear minLengthProc minProcLoss NUM opts ptb ptb2 ptb_sel T uniPidx yesno

% add pert number to PT_perc and PT_zero
if ~isempty(PT_zero) && ~isempty(PT_all)
    A = PT_all{:,2:end};
    B = reshape(A',1,8*size(A,1));
    C = reshape(B',2,size(B,2)/2)';
    D = [sortrows(repmat(PT_all.PertNo,4,1)) C];
    D = D(~isnan(D(:,2)),:);
    D = array2table(D,'VariableNames',{'NoPert','DIMStart','DIMEnd'});

    PT_zero = innerjoin(PT_zero,D,'Keys',{'DIMStart','DIMEnd'});
    PT_zero = movevars(PT_zero,'NoPert','Before','Quarter');
end
clear A B C D PT_del test


% overlap between PT_perc & PT_zero & assign NoPert
if ~isempty(PT_perc)
    PT_perc.yesno(:,1) = 0;
    PT_perc.NoPert(:,1) = 0;
end

if ~isempty(PT_perc) && ~isempty(PT_zero)
    for i = 1:height(PT_perc)
        yesno = find(PT_zero.Quarter == PT_perc.Quarter(i) & ...
                     PT_zero.IndexStart <= PT_perc.IndexStart(i) &...
                     PT_zero.IndexEnd >= PT_perc.IndexEnd(i));
        PT_perc.yesno(i) = length(yesno);
        if ~isempty(yesno)
            PT_perc.NoPert(i) = mean(PT_zero.NoPert(yesno));
        end
    end
end

if ~isempty(PT_perc)
    % select / remove if not both criteria
    PT_perc = removevars(PT_perc,{'Index'});
    PT_perc.Properties.VariableNames{7} = 'OverlapCritDIM';
    PT_perc = movevars(PT_perc,'NoPert','Before','Quarter');
end

if ~isempty(PT_zero)
    % plot last = overlapping of all 
    if makefig == 1
        for i = 1:height(PT_zero)
            subplot(2,2,PT_zero.Quarter(i))
            plot([PT_zero.DIMStart(i) PT_zero.DIMEnd(i)],[60 60],'-','LineWidth',6,'Color',[0.8 0 0.8])
        end
    end
end

% add milk losses to PT_all
if ~isempty(PT_all)
    PT_all.ML1a(:,1) = NaN;       % prepare absolute loss Q1
    PT_all.ML2a(:,1) = NaN;       % prepare absolute loss Q2
    PT_all.ML3a(:,1) = NaN;       % prepare absolute loss Q3
    PT_all.ML4a(:,1) = NaN;       % prepare absolute loss Q4
    PT_all.ML1r(:,1) = NaN;       % prepare relative loss Q1
    PT_all.ML2r(:,1) = NaN;       % prepare relative loss Q2
    PT_all.ML3r(:,1) = NaN;       % prepare relative loss Q3
    PT_all.ML4r(:,1) = NaN;       % prepare relative loss Q4
    PT_all.ML1l(:,1) = NaN;       % prepare lactation loss Q1
    PT_all.ML2l(:,1) = NaN;       % prepare lactation loss Q2
    PT_all.ML3l(:,1) = NaN;       % prepare lactation loss Q3
    PT_all.ML4l(:,1) = NaN;       % prepare lactation loss Q4
    PT_all.smML1a(:,1) = NaN;     % prepare absolute loss Q1 ~smooth
    PT_all.smML2a(:,1) = NaN;     % prepare absolute loss Q2 ~smooth
    PT_all.smML3a(:,1) = NaN;     % prepare absolute loss Q3 ~smooth
    PT_all.smML4a(:,1) = NaN;     % prepare absolute loss Q4 ~smooth
    PT_all.smML1r(:,1) = NaN;     % prepare relative loss Q1 ~smooth
    PT_all.smML2r(:,1) = NaN;     % prepare relative loss Q2 ~smooth
    PT_all.smML3r(:,1) = NaN;     % prepare relative loss Q3 ~smooth
    PT_all.smML4r(:,1) = NaN;     % prepare relative loss Q4 ~smooth
    PT_all.smML1l(:,1) = NaN;     % prepare lactation loss Q1 ~smooth
    PT_all.smML2l(:,1) = NaN;     % prepare lactation loss Q2 ~smooth
    PT_all.smML3l(:,1) = NaN;     % prepare lactation loss Q3 ~smooth
    PT_all.smML4l(:,1) = NaN;     % prepare lactation loss Q4 ~smooth

    for i = 1:height(PT_all) % all perturbations
        for j = 1:4     % all quarters
            if ~isnan(PT_all{i,1+2*j-1})  % if perturbation detected
                idx = [0 0];
                idx(1) = find(DIM == PT_all{i,1+2*j-1});    % idx of DIM start  
                idx(2) = find(DIM == PT_all{i,1+2*j});      % idx of DIM end
                PT_all{i,9+j} = -sum(QMY(idx(1):idx(2),j)-MOD(idx(1):idx(2),j)); % absolute losses = total produced - total expected in kg
                PT_all{i,13+j} = 100-(sum(QMY(idx(1):idx(2),j))/sum(MOD(idx(1):idx(2),j))*100); % relative losses = 100%- total produced / total expected
                PT_all{i,17+j} = (PT_all{i,9+j}/sum(MOD(:,j)))*100; % (total produced-total expected)/(lactation expected)
                PT_all{i,21+j} = -sum(RES2(idx(1):idx(2),j)); % absolute losses = total produced - total expected in kg
                PT_all{i,25+j} = -sum(RES2(idx(1):idx(2),j))/sum(MOD(idx(1):idx(2),j))*100; % relative losses = total produced / total expected
                PT_all{i,29+j} = (PT_all{i,21+j}/sum(MOD(:,j)))*100; % (total produced-total expected)/(lactation expected)
            else
                idx = [0 0];
                idx(1) = find(DIM == min(PT_all{i,2:9}));    % idx of DIM start  
                idx(2) = find(DIM == max(PT_all{i,2:9}));    % idx of DIM end

                PT_all{i,9+j} = -sum(QMY(idx(1):idx(2),j)-MOD(idx(1):idx(2),j)); % absolute losses = total produced - total expected in kg
                PT_all{i,13+j} = 100-(sum(QMY(idx(1):idx(2),j))/sum(MOD(idx(1):idx(2),j))*100); % relative losses = 100%- total produced / total expected
                PT_all{i,17+j} = length(find(RES2(idx(1):idx(2),j)<0))./diff(idx)*100; % percentage of negative residuals
            end
        end       
    end
end

% add milk losses to PT_perc
if ~isempty(PT_perc)
    PT_perc.MLa(:,1) = NaN;     % prepare absolute loss
    PT_perc.MLr(:,1) = NaN;     % prepare relative loss
    PT_perc.MLl(:,1) = NaN;     % prepare lactation loss
    PT_perc.smMLa(:,1) = NaN;     % prepare absolute loss
    PT_perc.smMLr(:,1) = NaN;     % prepare relative loss
    PT_perc.smMLl(:,1) = NaN;     % prepare lactation loss

    for i = 1:height(PT_perc) % all perturbations
        idx = [0 0];
        idx(1) = PT_perc.IndexStart(i);    % idx of DIM start
        idx(2) = PT_perc.IndexEnd(i);      % idx of DIM end
        j = PT_perc.Quarter(i);            % Quarter 
        PT_perc.MLa(i) = -sum(QMY(idx(1):idx(2),j)-MOD(idx(1):idx(2),j)); % absolute losses = total produced - total expected in kg
        PT_perc.MLr(i) = 100-(sum(QMY(idx(1):idx(2),j))/sum(MOD(idx(1):idx(2),j))*100); % relative losses = total produced / total expected
        PT_perc.MLl(i) = -(PT_perc.MLa(i)/sum(MOD(:,j)))*100; % (total produced-total expected)/(lactation expected)
        PT_perc.smMLa(i) = -sum(RES2(idx(1):idx(2),j)); % absolute losses = total produced - total expected in kg
        PT_perc.smMLr(i) = -sum(RES2(idx(1):idx(2),j))/sum(MOD(idx(1):idx(2),j))*100; % relative losses = total produced / total expected
        PT_perc.smMLl(i) = (PT_perc.smMLa(i)/sum(MOD(:,j)))*100; % (total produced-total expected)/(lactation expected)
    end
else
    PT_zero.NoPert(:,1) = 0;
    PT_zero = movevars(PT_zero,'NoPert','Before','Quarter'); % add that no pert number could be assigned
end

% add milk losses to PT_zero
if ~isempty(PT_zero)
    PT_zero.MLa(:,1) = NaN;     % prepare absolute loss
    PT_zero.MLr(:,1) = NaN;     % prepare relative loss
    PT_zero.MLl(:,1) = NaN;     % prepare lactation loss

    for i = 1:height(PT_zero) % all perturbations
        idx = [0 0];
        idx(1) = PT_zero.IndexStart(i);    % idx of DIM start
        idx(2) = PT_zero.IndexEnd(i);      % idx of DIM end
        j = PT_zero.Quarter(i);            % Quarter 
        PT_zero.MLa(i) = sum(QMY(idx(1):idx(2),j)-MOD(idx(1):idx(2),j)); % absolute losses = total produced - total expected in kg
        PT_zero.MLr(i) = 100-(sum(QMY(idx(1):idx(2),j))/sum(MOD(idx(1):idx(2),j))*100); % relative losses = total produced / total expected
        PT_zero.MLl(i) = (PT_zero.MLa(i)/sum(MOD(:,j)))*100; % (total produced-total expected)/(lactation expected)
        PT_zero.smMLa(i) = -sum(RES2(idx(1):idx(2),j)); % absolute losses = total produced - total expected in kg
        PT_zero.smMLr(i) = -sum(RES2(idx(1):idx(2),j))/sum(MOD(idx(1):idx(2),j))*100; % relative losses = total produced / total expected
        PT_zero.smMLl(i) = (PT_zero.smMLa(i)/sum(MOD(:,j)))*100; % (total produced-total expected)/(lactation expected)
    end
end

% put warnings on again
warning('on','all')

















































