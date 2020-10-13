%% S1_CaseSelection / data
% This script will load the data, select cases,  check for overlap with DAY
% and MILK data and uniformize IDs
% 
% STEP 1: load datasets
%           dir = MILK LOSSES scripts
%           D1 = data treatment registers
%           D2 = data MPR
%           D3 = quarter milk yield
%           D4 = daily milk yield
%
% STEP 2A: Select treatment ~ overlap DAY / MILK
% STEP 2B: Select MPR
% STEP 2C: Select QMY/QEC cases

clear variables
close all
clc


%% STEP 1: load data
% directory
%dir = 'C:\Users\u0084712\Documents\Box Sync\Documents\MastiMan\Research\Data mining\MILK LOSSES scripts\'
dir = 'C:\Users\u0084712\Documents\MastiMan\Research\dataMining\MILK LOSSES scripts\A1_results\';

% files
load([dir 'D1_CM.mat'])  % Data treatment registers
load([dir 'D2_MPR.mat'])



% =============== Not sure DAY will be used ============= %
% % % % day data
% % % FN = ls([dir 'DAY_*']);
% % % for i = 1:length(FN)
% % %     numLoc = regexp(FN(i,:),'_');       % this functions finds the unique positions of 2 successive numbers in a filename
% % %     farmsDay{i} = FN(i,numLoc(1)+1:numLoc(2)-1);
% % %     disp(['Farm = ' farms{i}])
% % %  
% % %     opts = detectImportOptions([dir deblank(FN(i,:))]); % get Import options
% % %     day.(farms{i}) = readtable([dir deblank(FN(i,:))],opts); % readtable 
% % % end
% =============== Not sure DAY will be used ============= %






%% Cases Quarter Milk Yield
% load data
% model
% detect perturbations
% save outcomes
clc
clear variables
close all

% milk data
dir = 'C:\Users\u0084712\Documents\MastiMan\Research\dataMining\MILK LOSSES scripts\A1_results\';
FN = ls([dir 'MILK_*']);
for i = 91:length(FN)
    numLoc = regexp(FN(i,:),'_');       % this functions finds the unique positions of 2 successive numbers in a filename
    farms{i} = FN(i,numLoc(1)+1:numLoc(2)-1);
    disp(['Load data, Farm = ' farms{i}])
 
    opts = detectImportOptions([dir deblank(FN(i,:))]); % get Import options
    milk.(farms{i}) = readtable([dir deblank(FN(i,:))],opts); % readtable 
    
    % select
%     cowsel = 1:10;%[23 24 25 26 28 70 71 104 105 111 112];
    cowlac = unique(milk.(farms{i})(:,[2 6]),'rows');
%     cowlac = cowlac(cowsel,:); % select a number interesting lactations to train
    % select data of selected cows
    milk.(farms{i}) = innerjoin(milk.(farms{i}),cowlac,'Keys',{'B_ID','Lac'}); % select
    milk.(farms{i}) = milk.(farms{i})(milk.(farms{i}).DIM < 306,:);  % trim lactations up to 305 days
    [~,ind] = unique(milk.(farms{i})(:,[2 6 7]),'rows');
    milk.(farms{i}) = milk.(farms{i})(ind,:);               % unique milkings
    milk.(farms{i}).MI(milk.(farms{i}).MI > 50) = NaN;      % MI > 50 = NaN
    clear opts cowlac cowsel fileName    
end
clear i FN dir numLoc opts 


% run model function
farms = fieldnames(milk);
for j = 1:length(farms)
    disp(['Model fit, farm name = ' farms{j}])
    
    % prepare savename
    savename = ['C:\Users\u0084712\Documents\MastiMan\Research\mastitisMilkLoss\A1_results\' farms{j}];
    makefig = 0;    % overruled inside the function
    
    % run function to model quarters
    tic 
    [moddata.(farms{j}), model.(farms{j})] = F1_modelQuarter(milk.(farms{j}),makefig,savename); % with fig, time for 10 = 48s, without = 0.86s
    toc
end

% run perturbation function
farms = fieldnames(moddata);
clear PT_zero PT_all PT_perc PT_zero2 PT_all2 PT_perc2
for j =9 1:length(farms)
    disp(['Perturbations, farm name = ' farms{j}])
    savename = ['C:\Users\u0084712\Documents\MastiMan\Research\mastitisMilkLoss\A1_results\' farms{j}];

    % summarize cowlac
    cowlac = array2table(unique([moddata.(farms{j}).B_ID,moddata.(farms{j}).Lac],'rows'),'VariableNames',{'B_ID','Lac'});
    
    % prepare arrays minProcLoss = 80,minLengthNeg = 10,minLengthProc = 5
    PT_zero.(farms{j}) = array2table(zeros(0,16),'VariableNames',{'B_ID','Lac','NoPert','Quarter',...
                        'IndexStart','DIMStart','IndexEnd','DIMEnd','PertLength','NoNegEp',...
                        'MLa','MLr','MLl','smMLa','smMLr','smMLl'});
    PT_perc.(farms{j}) = array2table(zeros(0,16),'VariableNames',{'B_ID','Lac','NoPert','Quarter',...
                        'IndexStart','DIMStart','IndexEnd','DIMEnd','PertLength','OverlapCritDIM',...
                        'MLa','MLr','MLl','smMLa','smMLr','smMLl'});
    PT_all.(farms{j}) = array2table(zeros(0,35),'VariableNames',{'B_ID','Lac','PertNo','DIMstart1',...
                        'DIMend1','DIMstart2','DIMend2','DIMstart3','DIMend3','DIMstart4','DIMend4',...
                        'ML1a','ML2a','ML3a','ML4a','ML1r','ML2r','ML3r','ML4r','ML1l','ML2l','ML3l','ML4l',...
                        'smML1a','smML2a','smML3a','smML4a','smML1r','smML2r',...
                        'smML3r','smML4r','smML1l','smML2l','smML3l','smML4l'});
    
    % prepare arrays minProcLoss = 80,minLengthNeg = 5,minLengthProc = 2
    PT_zero2.(farms{j}) = array2table(zeros(0,16),'VariableNames',{'B_ID','Lac','NoPert','Quarter',...
                        'IndexStart','DIMStart','IndexEnd','DIMEnd','PertLength','NoNegEp',...
                        'MLa','MLr','MLl','smMLa','smMLr','smMLl'});
    PT_perc2.(farms{j}) = array2table(zeros(0,16),'VariableNames',{'B_ID','Lac','NoPert','Quarter',...
                        'IndexStart','DIMStart','IndexEnd','DIMEnd','PertLength','OverlapCritDIM',...
                        'MLa','MLr','MLl','smMLa','smMLr','smMLl'});
    PT_all2.(farms{j}) = array2table(zeros(0,35),'VariableNames',{'B_ID','Lac','PertNo','DIMstart1',...
                        'DIMend1','DIMstart2','DIMend2','DIMstart3','DIMend3','DIMstart4','DIMend4',...
                        'ML1a','ML2a','ML3a','ML4a','ML1r','ML2r','ML3r','ML4r','ML1l','ML2l','ML3l','ML4l',...
                        'smML1a','smML2a','smML3a','smML4a','smML1r','smML2r',...
                        'smML3r','smML4r','smML1l','smML2l','smML3l','smML4l'});
    
    % estimate perturbations
    count = -1; % for display purposes
    for i = 1:height(cowlac)
        % find data of cow lac
        ind = find(moddata.(farms{j}).B_ID == cowlac.B_ID(i) & moddata.(farms{j}).Lac == cowlac.Lac(i));
        
        % select data for function
        set.DIM = moddata.(farms{j}).DIM(ind);
        set.QMY = [moddata.(farms{j}).MYLF(ind) moddata.(farms{j}).MYRF(ind) ...
                   moddata.(farms{j}).MYLH(ind) moddata.(farms{j}).MYRH(ind)];
        set.MOD = [moddata.(farms{j}).ModIT1(ind) moddata.(farms{j}).ModIT2(ind) ...
                   moddata.(farms{j}).ModIT3(ind) moddata.(farms{j}).ModIT4(ind)];

        set.perc = floor(i/height(cowlac)*100);
        if ismember(set.perc,0:10:100) && set.perc > count
            toc
            disp(['Calculate losses = ' num2str(set.perc) '% completed, i = ' num2str(i) ' out of ' num2str(height(cowlac))])
            count = set.perc;
            tic
        end
               
        % only run the loop if model data available
        if sum(sum(~isnan(set.MOD))) > 0

            % set makefig
            if i==200000000 %sum(ismember(i,[22 44 66]))>0
                makefig = 1;
                % detect perturbations
                [Pzero,Pperc,Pall,h] = F2_pertQuarter(set.DIM,set.QMY,set.MOD,80,10,5,makefig);
                [Pzero2,Pperc2,Pall2,h2] = F2_pertQuarter(set.DIM,set.QMY,set.MOD,70,7,3,makefig);
            else 
                makefig = 0;
                % detect perturbations
                [Pzero,Pperc,Pall] = F2_pertQuarter(set.DIM,set.QMY,set.MOD,80,10,5,makefig);
                [Pzero2,Pperc2,Pall2] = F2_pertQuarter(set.DIM,set.QMY,set.MOD,70,7,3,makefig);
            end

            % store perturbations in arrays 
            PT_zero.(farms{j}) = [PT_zero.(farms{j});repmat(cowlac(i,:),height(Pzero),1) Pzero];
            PT_perc.(farms{j}) = [PT_perc.(farms{j}); repmat(cowlac(i,:),height(Pperc),1) Pperc];
            PT_all.(farms{j}) = [PT_all.(farms{j}); repmat(cowlac(i,:),height(Pall),1) Pall];
            % store perturbations in arrays 
            PT_zero2.(farms{j}) = [PT_zero2.(farms{j});repmat(cowlac(i,:),height(Pzero2),1) Pzero2];
            PT_perc2.(farms{j}) = [PT_perc2.(farms{j}); repmat(cowlac(i,:),height(Pperc2),1) Pperc2];
            PT_all2.(farms{j}) = [PT_all2.(farms{j}); repmat(cowlac(i,:),height(Pall2),1) Pall2];

            if makefig == 1
    %             saveas(h,[savename 'pert_Cow' num2str(cowlac.B_ID(i)) '_Lac_' num2str(cowlac.Lac(i)) '.fig'])
                saveas(h,[savename '_pert_Cow' num2str(cowlac.B_ID(i)) '_Lac' num2str(cowlac.Lac(i)) '.tif'])
    %             saveas(h2,[savename 'pert_Cow' num2str(cowlac.B_ID(i)) '_Lac_' num2str(cowlac.Lac(i)) '_2' '.fig'])
                saveas(h2,[savename '_pert_Cow' num2str(cowlac.B_ID(i)) '_Lac' num2str(cowlac.Lac(i)) '_2' '.tif'])
                close all
            end
        end
    end
end

% clear variables
clear cowlac h h2 ind j makefig Pall Pal2 Pperc Pperc2
clear Pzero Pzero2 savename set i Pall2

% save settings
settings1 = [80 10 5];  % settings1
settings2 = [70 7 3];   % settings2

% save statements
save('D1_PT.mat','PT_all','PT_perc','PT_zero','settings1')
save('D2_PT2.mat','PT_all2','PT_perc2','PT_zero2','settings2')
save('D3_Mod.mat','model')

% clear variables
clear count PT_perc PT_perc2 PT_zero PT_zero2 PT_all PT_all2 settings1 settings2 model 
clear DIM idx1 idx2 ii k L minLengthNeg minLengthProc minProcLoss MOD NUM ptb ptb2 QMY RES RES2 sortedindx sst uniPidx win


%% MPR data case selection
% clear variables
close all
clc

% directory
dir = 'C:\Users\u0084712\Documents\MastiMan\Research\dataMining\MILK LOSSES scripts\';

% load MPR file
load([dir 'D2_MPR.mat'])
clear dir

% select the variables needed
MPR = MPR(:,[2 3 4 6 7 9 10 12 15 16 17 18]);

% summary
SUM = unique(MPR(:,[1 2 5]),'rows','first');

% summarize high SCC
tic
for i = 1:height(SUM)
    ind = find(startsWith(MPR.FarmName,SUM.FarmName(i)) & startsWith(MPR.OffRegNo,SUM.OffRegNo(i)) & MPR.Lac == SUM.Lac(i));
    SUM.NoMeas(i,1) = length(ind);
    SUM.FirstMeas(i,1) = MPR.DIM(ind(1));
    SUM.LastMeas(i,1) = MPR.DIM(ind(end));
    SUM.No150(i,1) = length(find(MPR.SCC(ind)>150));
    SUM.No250(i,1) = length(find(MPR.SCC(ind)>250));
    SUM.StartDate(i,1) = MPR.Date(ind(1));
    SUM.EndDate(i,1) = MPR.Date(ind(end));
    SUM.MaxSCC(i,1) = max(MPR.SCC(ind));
end
toc

% number of unique cows
COW = unique(MPR(:,[1 2]),'rows');

% sum
SUM.Select(:,1) = 0;
SUM.Select(SUM.NoMeas >= 5 & SUM.FirstMeas <= 42 & SUM.LastMeas >= 150,1) = 1; 
disp(['No Lact. selected = ' num2str(sum(SUM.Select))])


% select MPR data
MPRall = MPR;
MPR = innerjoin(MPRall,SUM(SUM.Select == 1,[1 2 3]),'Keys',{'FarmName','OffRegNo','Lac'});

load('D8_MPRcases.mat');
% select cases
CaseMPR.T150_T250_1 = F5_caseMPR(MPR,150,250,1,1);
CaseMPR.T150_T250_2 = F5_caseMPR(MPR,150,250,2,8*7);
CaseMPR.T250_T350_1 = F5_caseMPR(MPR,250,350,1,1);
CaseMPR.T350_T450_1 = F5_caseMPR(MPR,350,450,1,1);
CaseMPR.T450_T550_1 = F5_caseMPR(MPR,450,550,1,1);
CaseMPR.T550_T650_1 = F5_caseMPR(MPR,550,650,1,1);
CaseMPR.T650_T750_1 = F5_caseMPR(MPR,650,750,1,1);
CaseMPR.T200_T200_1 = F5_caseMPR(MPR,200,200,1,1);
CaseMPR.T500_T500_1 = F5_caseMPR(MPR,500,500,1,1);
CaseMPR.T800_T800_1 = F5_caseMPR(MPR,800,800,1,1);
CaseMPR.T1000_T1000_1 = F5_caseMPR(MPR,1000,1000,1,1);
CaseMPR.T1500_T1500_1 = F5_caseMPR(MPR,1500,1500,1,1);
CaseMPR.T2000_T2000_1 = F5_caseMPR(MPR,2000,2000,1,1);


% save dataset
save('D8_MPRcases.mat','CaseMPR', 'MPR','MPRall','SUM');

clear ans Cases COW cowlac dif i idx ind indx T T2 T1


%% Add animal identification

% load all ID information
load('D4_IDs.mat')

% files
dir = 'C:\Users\u0084712\Documents\MastiMan\Research\dataMining\MILK LOSSES scripts\';
load([dir 'D1_CM.mat'])  % Data treatment registers
clear dir

% change names and adjust
id.Debrabander = id.DeBrabander;
id.VandenMeijenberg = id.Vandenmeijdenberg;
id.Vereenooghe = id.Vanhullebus;
CLMAS = rmfield(CLMAS,'Dewulf');
CLMAS = rmfield(CLMAS,'ZBentley');
CLMAS = rmfield(CLMAS,'ZHunt');

% sum cases
farms = fieldnames(CLMAS);
SUM = array2table((1:length(farms))','VariableNames',{'FarmID'});
SUM.FarmName(:,1) = farms';

% try to find the right IDs and merge
for i = 1:length(farms)
  
    % prepare adding Backup ID
    CLMAS.(farms{i}).B_ID(:,1) = NaN; % set NaN
    if sum(contains(CLMAS.(farms{i}).Properties.VariableNames,'Name')) == 0
        CLMAS.(farms{i}).Name(:,1) = {'NoName'};
    end
    
    % find corresponding IDs and check date
    for j = 1:height(CLMAS.(farms{i}))
        % try matching with OffRegNo
        ind = find(contains(id.(farms{i}).OfficialRegNo,CLMAS.(farms{i}).OffRegNo(j)) == 1 & ...
                   datenum(id.(farms{i}).BDate) < datenum(CLMAS.(farms{i}).Date(j)));
        
        % try matching with Name
        if isempty(ind)
            ind = find(contains(id.(farms{i}).Name,CLMAS.(farms{i}).Name(j)) == 1 & ...
                   datenum(id.(farms{i}).BDate) < datenum(CLMAS.(farms{i}).Date(j)));
        end
        
        % try matching with CowID
        if isempty(ind)
            ind = find(id.(farms{i}).CowID == CLMAS.(farms{i}).CowID(j) & ...
                   datenum(id.(farms{i}).BDate) < datenum(CLMAS.(farms{i}).Date(j)));
        end   
        
        % match cases
        if ~isempty(ind)
            CLMAS.(farms{i}).B_ID(j) = id.(farms{i}).B_ID(ind(1));
        end
    end
    
    % summarize the number and percentage of cases matched
    SUM.NoCases(i,1) = height(CLMAS.(farms{i}));
    SUM.CasesMatched(i,1) = sum(~isnan(CLMAS.(farms{i}).B_ID));
    SUM.PercMatched(i,1) = (SUM.CasesMatched(i,1)/SUM.NoCases(i,1))*100;
      
    % take out the unmatched cases
    Unmatched.(farms{i}) = CLMAS.(farms{i})(isnan(CLMAS.(farms{i}).B_ID),:);
    if isempty(Unmatched.(farms{i}))
        Unmatched = rmfield(Unmatched,farms{i});
    end
end

SUM.CasesUnmatched(:,1) =  SUM.NoCases - SUM.CasesMatched;

disp(['Total number cases in treatment registers = ' num2str(sum(SUM.NoCases))])
disp(['Total number cases matched with backups = ' num2str(sum(SUM.CasesMatched))...
      ' of ' num2str(length(SUM.NoCases(SUM.CasesMatched>0))) ' farms'])
disp(['Total number unmatched = ' num2str(sum(SUM.CasesUnmatched))])

save('D9_TR_cases.mat','CLMAS')

% We have now matched the backup IDs with the clinical mastitis cases
% Next, we want to check whether we have day or milk information to
% calculate the milk losses


%% DAILY DATA

% >>>> ========== DAILY DATA - LOAD DATA ========== <<<<
dir = 'C:\Users\u0084712\Documents\MastiMan\Research\dataMining\MILK LOSSES scripts\A1_results\';
FN = ls([dir 'DAY_*']);
for i = 1:length(FN)
    numLoc = regexp(FN(i,:),'_');       % this functions finds the unique positions of 2 successive numbers in a filename
    farms{i} = FN(i,numLoc(1)+1:numLoc(2)-1);
    disp(['Load data, Farm = ' farms{i}])
 
    opts = detectImportOptions([dir deblank(FN(i,:))]); % get Import options
    day.(farms{i}) = readtable([dir deblank(FN(i,:))],opts); % readtable 
    
    % select
%     cowsel = 1:10;%[23 24 25 26 28 70 71 104 105 111 112];
    cowlac = unique(day.(farms{i})(:,[2 8]),'rows');
%     cowlac = cowlac(cowsel,:); % select a number interesting lactations to train
    % select data of selected cows
    day.(farms{i}) = innerjoin(day.(farms{i}),cowlac,'Keys',{'B_ID','Lac'}); % select
    day.(farms{i}) = day.(farms{i})(day.(farms{i}).DIM < 306,:);  % trim lactations up to 305 days
    [~,idx] = unique(day.(farms{i})(:,[2 9 10]),'rows');
    day.(farms{i}) = day.(farms{i})(idx,:);

    clear opts cowlac cowsel fileName    
end
clear i FN dir numLoc opts idx


% >>>> ========== DAILY DATA - MODEL ITW ========== <<<<
farms = fieldnames(day);
for j = 1:length(farms)
    disp(['Model fit DAY, farm name = ' farms{j}])
    
    % prepare savename
    savename = ['C:\Users\u0084712\Documents\MastiMan\Research\mastitisMilkLoss\A1_results\' farms{j}];
    
    % prepare make fig or not >> overruled in function
    makefig = 0;

    % run function to model quarters
    tic 
    [moddataDAY.(farms{j}), modelDAY.(farms{j})] = F3_modelDay(day.(farms{j}),makefig,savename); % with fig, time for 10 = 48s, without = 0.86s    
    toc
    close all
end
clear farms id j makefig NoMeas uni THapart


% >>>> ========== DAILY DATA - DETECT PERTURBATIONS ========== <<<<
farms = fieldnames(moddataDAY);
clear PT_zero PT_all PT_perc PT_zero2 PT_all2 PT_perc2
for j = 1:length(farms)
    disp(['Perturbations, farm name = ' farms{j}])
    savename = ['C:\Users\u0084712\Documents\MastiMan\Research\mastitisMilkLoss\A1_results\' farms{j}];

    % summarize cowlac
    cowlac = array2table(unique([moddata.(farms{j}).B_ID,moddata.(farms{j}).Lac],'rows'),'VariableNames',{'B_ID','Lac'});
    
    % prepare arrays minProcLoss = 80,minLengthNeg = 10,minLengthProc = 5
    PT_zeroDAY.(farms{j}) = array2table(zeros(0,11),'VariableNames',{'B_ID','Lac','NoPert',...
                        'IndexStart','DIMStart','IndexEnd','DIMEnd','PertLength',...
                        'MLa','MLr','MLl'});
    
    % prepare arrays minProcLoss = 80,minLengthNeg = 5,minLengthProc = 2
    PT_zeroDAY2.(farms{j}) = array2table(zeros(0,11),'VariableNames',{'B_ID','Lac','NoPert',...
                        'IndexStart','DIMStart','IndexEnd','DIMEnd','PertLength',...
                        'MLa','MLr','MLl'});
    
    % estimate perturbations
    count = -1; % for display purposes
    for i = 1:height(cowlac)
        % find data of cow lac
        ind = find(moddataDAY.(farms{j}).B_ID == cowlac.B_ID(i) & moddataDAY.(farms{j}).Lac == cowlac.Lac(i));
        
        % select data for function
        set.DIM = moddataDAY.(farms{j}).DIM(ind);
        set.MY = moddataDAY.(farms{j}).MYs(ind);
        set.MOD = moddataDAY.(farms{j}).ModIT(ind);

        set.perc = floor(i/height(cowlac)*100);
        if ismember(set.perc,0:10:100) && set.perc > count
            toc
            disp(['Calculate losses = ' num2str(set.perc) '% completed, i = ' num2str(i) ' out of ' num2str(height(cowlac))])
            count = set.perc;
            tic
        end
               
        % only run the loop if model data available
        if sum(sum(~isnan(set.MOD))) > 0

            % set makefig
            if sum(ismember(i,[22 44 66]))>0
                makefig = 0;
                % detect perturbations
                [Pzero,h] = F4_pertDAY(set.DIM,set.MY,set.MOD,80,5,makefig);
                [Pzero2,h2] = F4_pertDAY(set.DIM,set.MY,set.MOD,70,5,makefig);
            else 
                makefig = 0;
                % detect perturbations
                Pzero = F4_pertDAY(set.DIM,set.MY,set.MOD,80,5,makefig);
                Pzero2 = F4_pertDAY(set.DIM,set.MY,set.MOD,70,5,makefig);
            end

            % store perturbations in arrays 
            PT_zeroDAY.(farms{j}) = [PT_zeroDAY.(farms{j});repmat(cowlac(i,:),height(Pzero),1) Pzero];
            % store perturbations in arrays 
            PT_zeroDAY2.(farms{j}) = [PT_zeroDAY2.(farms{j});repmat(cowlac(i,:),height(Pzero2),1) Pzero2];

            if makefig == 1
    %             saveas(h,[savename 'pert_Cow' num2str(cowlac.B_ID(i)) '_Lac_' num2str(cowlac.Lac(i)) '.fig'])
                saveas(h,[savename '_pert_Cow' num2str(cowlac.B_ID(i)) '_Lac' num2str(cowlac.Lac(i)) '.tif'])
    %             saveas(h2,[savename 'pert_Cow' num2str(cowlac.B_ID(i)) '_Lac_' num2str(cowlac.Lac(i)) '_2' '.fig'])
                saveas(h2,[savename '_pert_Cow' num2str(cowlac.B_ID(i)) '_Lac' num2str(cowlac.Lac(i)) '_2' '.tif'])
                close all
            end
        end
    end
end


% >>>> ========== DAILY DATA - SAVE PERTURBATIONS ========== <<<<
settings1 = [80,5];
settings2 = [70,5];
save('D5_PTday.mat','PT_zeroDAY','settings1')
save('D6_PTday2.mat','PT_zeroDAY2','settings2')
save('D7_Mod.mat','modelDAY')


clear ans count farms h h2 j ind makefig savename set settings1 settings2 Pzero Pzero2 i







%% MPR - Overview of cows with  high scc

ans =find(MPRi.SCC > 1000);
ind = find(contains(moddataDAY.Almey.OffRegNo,'BE 111940663') & moddataDAY.Almey.Lac == 3);
ind2 = find(contains(MPRi.OffRegNo,'BE 111940663') & MPRi.Lac == 3);
figure();
plot(moddataDAY.Almey.DIM(ind),moddataDAY.Almey.TDMY(ind));
hold on; 
plot(moddataDAY.Almey.DIM(ind),moddataDAY.Almey.ModIT(ind));
yyaxis right; plot(MPRi.DIM(ind2),MPRi.SCC(ind2),'s')










%% milk losses per case
PT_all.DeSchutter = PT_all.Deschutter;
PT_all = rmfield(PT_all,'Deschutter');
fn = fieldnames(PT_all);
for i = 1:length(fn)
    farm = fn{i};
    % no of perturbations
    sumary(i,1) = height(PT_all.(farm));
    
    % no pertu/lactation
    sumary(i,2) = height(PT_all.(farm))./height(unique(moddata.(farm)(:,[2 6]),'rows'));
    
    % average duration / milk loss q1
    ind= find(isnan(PT_all.(farm).DIMstart1)==0);
    sumary(i,3) = nanmean(PT_all.(farm).DIMend1(ind)-PT_all.(farm).DIMstart1(ind));
    sumary(i,7) = nanmean(PT_all.(farm).ML1a(ind));
    sumary(i,11) = nanmean(PT_all.(farm).ML1r(ind));
%     q2
    ind= find(isnan(PT_all.(farm).DIMstart2)==0);
    sumary(i,4) = nanmean(PT_all.(farm).DIMend2(ind)-PT_all.(farm).DIMstart2(ind));
    sumary(i,8) = nanmean(PT_all.(farm).ML2a(ind));
    sumary(i,12) = nanmean(PT_all.(farm).ML2r(ind));

%     q3
    ind= find(isnan(PT_all.(farm).DIMstart3)==0);
    sumary(i,5) = nanmean(PT_all.(farm).DIMend3(ind)-PT_all.(farm).DIMstart3(ind));
    sumary(i,9) = nanmean(PT_all.(farm).ML3a(ind));
    sumary(i,13) = nanmean(PT_all.(farm).ML3r(ind));

%     q4
    ind= find(isnan(PT_all.(farm).DIMstart4)==0);
    sumary(i,6) = nanmean(PT_all.(farm).DIMend4(ind)-PT_all.(farm).DIMstart4(ind));
    sumary(i,10) = nanmean(PT_all.(farm).ML4a(ind));
    sumary(i,14) = nanmean(PT_all.(farm).ML4r(ind));

    % all front
    ind= find(isnan(PT_all.(farm).DIMstart1)==0| isnan(PT_all.(farm).DIMstart2)==0);
    sumary(i,15) = nanmean([PT_all.(farm).DIMend1(ind)-PT_all.(farm).DIMstart1(ind); PT_all.(farm).DIMend2(ind)-PT_all.(farm).DIMstart2(ind)]);
    sumary(i,17) = nanmean([PT_all.(farm).ML1a(ind);PT_all.(farm).ML2a(ind)]);
    sumary(i,19) = nanmean([PT_all.(farm).ML1r(ind);PT_all.(farm).ML2r(ind)]);
    
    
    % all hind
    ind= find(isnan(PT_all.(farm).DIMstart3)==0| isnan(PT_all.(farm).DIMstart4)==0);
    sumary(i,16) = nanmean([PT_all.(farm).DIMend3(ind)-PT_all.(farm).DIMstart3(ind); PT_all.(farm).DIMend4(ind)-PT_all.(farm).DIMstart4(ind)]);
    sumary(i,18) = nanmean([PT_all.(farm).ML3a(ind);PT_all.(farm).ML4a(ind)]);
    sumary(i,20) = nanmean([PT_all.(farm).ML3r(ind);PT_all.(farm).ML4r(ind)]);

end

sumary(20,:) = [];

sumary = array2table(sumary,'VariableNames',{'Nopert','Nopertlac','Length1','Length2','Length3','Length4','MLa1','MLa2','MLa3','MLa4','MLr1','MLr2','MLr3','MLr4','LengthF','LengthH','MLaF','MLaH','MLrF','MLrH'});
    
    
sum(sumary.Nopert)

save('D12_pertQsum.mat','sumary')

























