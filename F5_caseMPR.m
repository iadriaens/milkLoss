function [Cases] = F5_caseMPR(MPR,T1,T2,NoMeas, THapart)
% this function will select cases from DHI/MR data based on threee
% different methods, and returns the identifiers plus dates
%
%
% INPUTS:   1   MPR     data containing ID, date, Lac and SCC
%           2   M       which method to use to select the cases
%           3   T1      threshold for 1st lactation cows
%           4   T2      threshold for 2+ lactation cows
%           5   NoMeas  number of successive measurements for a case
%
% OUTPUTS:  1   Cases   CowID / Date / method
%
%
% METHOD 1: M1  simple thresholds based on T1 and T2
%   STEP 1:     find the measurements above T
%   STEP 2:     list cases
%
% METHOD 2: M2  at least 2 successive measurments within 2 months above T1
%               and T2
%   STEP 1:     find the measurements above T
%   STEP 2:     check which measurements are successive and select
%   STEP 3:     list cases
%

%%%%----------- set up script ------------%%%%
% % % % T1 = 150;   % threshold for lactation 1
% % % % T2 = 250;   % threshold for higher lactations
% % % % NoMeas = 1; % number of successive above the threshold
% % % % THapart = 1*8*7;% number of days successives can be apart


%%%%----------- set up script ------------%%%%

%% STEP 1: find measurmeents per cows higher than thresholds

% find the indices of FarmName, OffRegNo and Lac
T(1) = find(startsWith(MPR.Properties.VariableNames,'FarmName')==1);
T(2) = find(startsWith(MPR.Properties.VariableNames,'OffRegNo')==1);
T(3) = find(endsWith(MPR.Properties.VariableNames,'Lac')==1);

% select unique farm/cow/lac in the dataset and add unique number to MPR
cowlac = unique(MPR(:,T),'rows');
cowlac.Uniques(:,1) = 1:height(cowlac);
MPR = innerjoin(MPR,cowlac,'Keys',MPR.Properties.VariableNames(T));

% sort the rows to ensure sort at farm / cow / lac level and add index
MPR = sortrows(MPR,T);
MPR.Index(:,1) = 1:height(MPR);

% find all measurements above threshold
MPR.aboveT(:,1) = 0;
MPR.aboveT((MPR.Lac == 1 & MPR.SCC > T1) | (MPR.Lac >1 & MPR.SCC > T2),1) = 1;


% find last measurement above threshold of each lactation
ind = find(MPR.aboveT == 1);
[~,uni] = unique(MPR(ind,T),'rows','last');
MPR.isLastAbove(:,1) = 0;
MPR.isLastAbove(ind(uni),1) = 1;

% find successive measurements above threshold in measurments
MPR.diffT(:,1) = 0;
dif = [diff(MPR.Index(MPR.aboveT == 1)); 0];  % find diff in index with next
MPR.diffT(MPR.aboveT == 1,1) = dif;         % add difference in index with next
MPR.diffT(MPR.isLastAbove == 1) = 0;        % if last delete

% find time between successive measurements above threshold
MPR.diffDays(:,1) = 0;
dif = [diff(datenum(MPR.Date(MPR.aboveT == 1))); 0];  % find diff in days with next
MPR.diffDays(MPR.aboveT == 1,1) = dif;         % add difference in index with next
MPR.diffDays(MPR.isLastAbove == 1) = 0;        % if last delete


% if NoMeas more than 1
idx(1) = find(startsWith(MPR.Properties.VariableNames,'CowID'));
idx(2) = find(startsWith(MPR.Properties.VariableNames,'Date'));
idx(3) = find(startsWith(MPR.Properties.VariableNames,'SCC'));
idx(4) = find(startsWith(MPR.Properties.VariableNames,'Index'));

if NoMeas > 1
    % select cases by detecting patterns
    pat = ones(NoMeas-1,1);
    ind = strfind(MPR.diffT',pat');
    indx = find(MPR.diffDays(ind) < THapart);
    Cases = MPR(ind(indx),[T(1:2) idx(1) T(3) idx(2:4)]);  % select cases
else
    Cases = MPR(MPR.aboveT == 1,[T(1:2) idx(1) T(3) idx(2:4)]);  % select cases
end 
 
% find number of successive measurements above threshold    
Cases.diff(:,1) = 0;
Cases.diff(1,1) = 10;
Cases.diff(2:end) = diff(Cases.Index);
Cases.New(Cases.diff ~= 1,1) = 1;

% add end of cases
ind = find(Cases.New == 0);
dif = diff(Cases.Index(ind));  % on if same case, higher if not same
% if Cases.New(end) == 0 % the last case belongs to other meas
    dif = [dif; 10];
% end
Cases.End(:,1) = 0;
Cases.End(ind(dif ~= 1)) = 1;

% find all with length > 1
ind = find(Cases.New(1:end-1) == 1 & Cases.New(2:end) == 0);
Cases.NoHigh(:,1) = 0;
Cases.NoHigh(Cases.New == 1,1) = NoMeas;
Cases.NoHigh(ind,1) = Cases.Index(Cases.End==1)-Cases.Index(ind)+NoMeas;

% add start and end date
Cases.StartDate(:,1) = NaT;
Cases.StartDate(Cases.New == 1,1) = Cases.Date(Cases.New == 1,1);
Cases.EndDate(:,1) = NaT;
ind = find(Cases.New == 1)+Cases.NoHigh(Cases.New ==1)-NoMeas;
Cases.EndDate(Cases.New == 1,1) = Cases.Date(ind,1);


% removevars
Cases = removevars(Cases, {'Index','New','diff','End'});












































