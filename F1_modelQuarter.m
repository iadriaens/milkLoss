function [data,mod] = F1_modelQuarter(data,makefig,savename)
% this function models the quarter milk yields and estimates the
% unperturbed curve at milking level
% INPUTS:   1   data:    contains at least MY at Q level, DIM, MI, B_ID, Lac
%                        order 4 consecutive QMY = LF / RF / LH / RH
%           2   makefig: whether or not figures are to be made / output
%           3   savedir: directory to save
%
% OUTPUTS:  1   data:    original data with additional columns added
%                        cols added = ModIT1 / ModIT12 / ModIT3 / ModIT4
%           2   mod:     information about the modelling   
%
% STEPS:
%
%   STEP 1:  prepare plotting and modelling wood
%   STEP 2:  detect unique cows in the dataset
%   STEP 3:  transform QMY for modeling to QMY/h*24
%   STEP 4:  detect outliers upper boundary using moving median
%   STEP 5:  initiate modelling iterations
%   STEP 6:  iterative modelling Wood & 1.6*std residuals
%   STEP 7:  store parameters and prepare outputs
%   STEP 8:  prepare figures and save in savedir
%   STEP 9:  save
%
% 
% ================= % for development purposes only % ================= %
% % % % % set file name to read data
% % % % fileName = ['C:\Users\u0084712\Documents\MastiMan\Research\dataMining\' ...
% % % %              'MILK LOSSES scripts\A1_results\MILK_Huzen_20060330_20191223.txt'];
% % % % %            'MILK LOSSES scripts\A1_results\MILK_Sanders_20161106_20200409.txt'];
% % % % cowsel = 730:750;%[23 24 25 26 28 70 71 104 105 111 112];
% % % %      
% % % % % read file
% % % % opts = detectImportOptions(fileName);
% % % % data = readtable(fileName,opts);
% % % % cowlac = unique(data(:,[2 6]),'rows');
% % % % cowlac = cowlac(cowsel,:); % select a number interesting lactations to train
% % % % 
% % % % % select data of selected cows
% % % % data = innerjoin(data,cowlac,'Keys',{'B_ID','Lac'}); % select
% % % % data = data(data.DIM < 306,:);  % trim lactations up to 305 days
% % % % clear opts cowlac cowsel fileName
% % % % 
% % % % % makefig prepare + savename
% % % % makefig = 1;
% % % % savename = 'C:\Users\u0084712\Documents\MastiMan\Research\mastitisMilkLoss\A1_results\Huzen_';
% ================= % for development purposes only % ================= %


%% STEP 0: Model and plot preparation

% prepare plotting
pl.LW = 'LineWidth';
pl.MS = 'MarkerSize';
pl.C = 'Color';
pl.MFC = 'MarkerFaceColor';

% prepare modeling Wood
mod.Wood = @(p,t) p(1).*t.^p(2).*exp(-p(3).*t);     % wood model - prepare
mod.opts = optimset('Display','off', 'MaxIter',1000);% set display to off
mod.UB = [25 1 0.01];                               % set upper boundary
mod.LB = [0 0 0];                                   % set lower boundary

% if longer than 50 hours, MI is set to NaN (missing milkings or new lac) 
data.MI(data.MI>30) = NaN;                          % set MI > 30 to NaN


%% STEP 1: INDIVIDUAL MODEL FIT/RESIDUALS -ITERATIVE
% summarize cow lactations based on B_ID and Lac vars
cowlac = array2table(unique([data.B_ID,data.Lac],'rows'),'VariableNames',{'B_ID','Lac'});

% modelling loop
for i = 1:height(cowlac)        % for all cow lactations entered
    % find data of the specific lactation
    T = find(startsWith(data.Properties.VariableNames,'MYLF'));  % index MY
    ind = find(data.B_ID == cowlac.B_ID(i) & data.Lac == cowlac.Lac(i) & ~isnan(sum(data{:,T:T+3},2)));
    
    
    % prepare make fig or not >>> OVERRULE INPUT
    if sum(ismember(i,[22 44 66 8 10]))>0
        makefig = 0;
    else
        makefig = 0;
    end

    
    % prepare figure of model
    if makefig == 1
        h=figure('Units','Normalized','Outerposition',[0 0 1 1]);
    end
    
    for j = 1:4     % for each quarter
        if nansum(data{ind,T+j-1}) > 10 && length(find(data{ind,T+j-1} >0))> 30 % at least 30 measurements and >10kg over whole lac
            % data selection
            DIM = data.DIM(ind);           % days in milk
            MI = min(20,data.MI(ind)); MI = max(MI,3);  % adjust extreme MI
            MY = (data{ind,T+j-1}./MI)*24; % milk yield / hour * 24
            
%           =============> OL detect via moving window <============
            [OL.LT, OL.low, OL.up,OL.cent] = isoutlier(MY(:,1),'movmedian',36,'thresholdfactor',4); % 1 week MM = 18meas
            OL.mad = movmad(MY(:,1)-OL.cent,36);
            OL.up = OL.cent + OL.mad*4;
%           =============> OL detect via moving window <============
            
            % select measurements excluding outliers
            MY(MY<OL.up,2) = 1;            % all measurement below the upper threshold
            
            if sum(MY<OL.up) > 1
                mod.my = MY(MY(:,2) ==1,1);    % select measurments
                mod.dim = DIM(MY(:,2) == 1);   % select dim
            else
                mod.my = 0;
                mod.dim = 0;
            end
            
            % visualisation of point selection
            if makefig == 1
                subplot(2,2,j); hold on; box on; axis tight;
                xlabel('DIM (days)'); ylabel('MY/h*24 (kg)')
                title(['i = ' num2str(i) ...
                    ', CowID = ' num2str(cowlac.B_ID(i)) ...
                    ', Lac = ' num2str(cowlac.Lac(i)) ...
                    ', Q = ' data.Properties.VariableNames{T+j-1}])
                plot(DIM,MY(:,1),'o:',pl.LW,1,pl.MS,3,pl.C,[0 0.80 0.8 1])
                plot(DIM, OL.up,pl.LW,1.2,pl.C, [1 0.4 0.2])
                plot(DIM, OL.cent ,pl.LW,1.7,pl.C, [1 0.6 0.2])
                plot(DIM(MY(:,2)==0), MY(MY(:,2)==0,1),'ks',pl.LW,2,pl.MS,3.5)
            end
            
            if length(mod.my) > 10 
                % model preparation
                mod.p0 = [nanmean(mod.my) 0.20 0.006];   % set initial parameters
                mod.p = lsqcurvefit(mod.Wood,mod.p0,mod.dim,mod.my,mod.LB,mod.UB,mod.opts); % estimate p
                mod.yval = mod.Wood(mod.p,mod.dim);      % y value initial wood fit

                % calculate std using all measurements after DIM=5
                mod.std = std(mod.my(mod.dim>=5) - mod.yval(mod.dim>=5)); % std of the residuals
                idx = find(mod.my-(mod.yval-1.6.*mod.std)<0 & ...
                           mod.dim > 5);  % this are the INITIAL measurements outside the range of 1.6 STD
                mod.rmse = sqrt(mean((mod.my-mod.yval).^2));   % calculate rmse
                mod.pIni(i,:) = [cowlac{i,:} mod.p mod.rmse];  % store initial p and rmse

                % visualisation initalisation modeling
                if makefig == 1
                    plot(mod.dim,mod.Wood(mod.p0,mod.dim),'b--',pl.LW,1.2);
                    plot(mod.dim,mod.Wood(mod.p,mod.dim),'r--',pl.LW,2);
                    plot(mod.dim(idx),mod.my(idx),'s','LineWidth',2,'MarkerSize',3.5, 'Color',[177/255 0 0 ]);
                end

                % preparation iterations
                mod.dimIT = mod.dim;            % prepare DIM
                mod.myIT = mod.my;              % prepare MY
                mod.stdIT = mod.std;            % prepare std
                mod.yvalIT = mod.yval;          % prepare y wood
                mod.pIT = mod.p;                % prepare p
                mod.rmseIT = 0;                 % prepare rmse
                mod.n = 0;                      % prepare no iter counter

                % iterations with while loop
                while abs(mod.rmse - mod.rmseIT) > 0.1 && mod.n <= 20
                    mod.n=mod.n+1;              % iteration n

                    % select measurements to delete
                    idx = find(mod.myIT > (mod.yvalIT-1.6*mod.stdIT) | mod.dimIT <= 5); % select measurements
                    mod.myIT = mod.myIT(idx);   % TMY of all milkings above yval+1.6*std or dim < 5
                    mod.dimIT = mod.dimIT(idx); % DIM of all milkings above yval+1.6*std or dim < 5

                    % refit model only when > 30 measurements are kept
                    if length(idx) > 30
                        mod.pIT = lsqcurvefit(mod.Wood,mod.pIT,mod.dimIT,mod.myIT,mod.LB,mod.UB,mod.opts);  % refit model
                        mod.yvalIT = mod.Wood(mod.pIT,mod.dimIT);      % calculate model values
                        mod.stdIT = std(mod.myIT(mod.dimIT>=5) - mod.yvalIT(mod.dimIT>=5)); % std of the residuals

                        mod.rmseIT = mod.rmse;   % Set the RMSE to the previous ones
                        mod.rmse = sqrt(mean((mod.myIT-mod.yvalIT).^2)); % new rmse
                    else
                        mod.n = 100;             % stop iteration by setting n=100
                    end

                    % visualisation iterations
                    if makefig == 1
                        plot(DIM,mod.Wood(mod.pIT,DIM),'--',pl.LW,1.5,pl.C,[0.6 0 0.8])
                        plot(mod.dimIT(mod.myIT <= (mod.yvalIT-1.6*mod.stdIT) & mod.dimIT > 5),...
                             mod.myIT(mod.myIT <= (mod.yvalIT-1.6*mod.stdIT) & mod.dimIT>5),...
                             's',pl.LW,2,pl.MS,3.5,pl.C,[0.6 0 0.8]);
                    end
                end

                % plot final model
                if makefig == 1
                    plot(DIM, mod.Wood(mod.pIT,DIM),'-',pl.LW,2.5,pl.C,[1 0 1])
                    ylim([0 max(MY(:,1)+0.5)]);
                end

                % final model
                mod.finMod = (mod.Wood(mod.pIT,data.DIM(ind))./24).*MI; % retransform for final model
                mod.rmseFin(i,2+j) = mod.rmse;      % rmse
                mod.nFin(i,2+j) = mod.n;            % n
                
            else
                mod.pIT = [NaN NaN NaN];            % no values
                mod.finMod = NaN*zeros(length(ind),1); % no values
                mod.rmseFin(i,2+j) = NaN;           % store final rmse
                mod.nFin(i,2+j) = NaN;              % store n iterations
            end
        else
            mod.pIT = [NaN NaN NaN];            % no values
            mod.finMod = NaN*zeros(length(ind),1); % no values
            mod.rmseFin(i,2+j) = NaN;           % store final rmse
            mod.nFin(i,2+j) = NaN;              % store n iterations
        end
        
        % final pars
        mod.pFin(i,[1 2]) = cowlac{i,:};        % store cowlac
        mod.pFin(i,2+3*j-2:2+3*j) = mod.pIT;    % store final p
        mod.rmseFin(i,[1 2]) = cowlac{i,:};     % store cowlac rmse
        mod.nFin(i,[1 2]) = cowlac{i,:};        % store cowlac n
        
        % prepare smoothing
        N = 7; % window of % residual smoothing
        
        % store vars
        if j == 1
            data.ModIT1(ind) = mod.finMod; % retransform
            data.ResIT1(ind) = data.MYLF(ind) - data.ModIT1(ind);        % calculate residuals
            data.ResPIT1(ind) = (data.ResIT1(ind)./data.ModIT1(ind)).*100; % calculate % residuals
            data.SmoRESPIT1(ind) = medfilt1(data.ResPIT1(ind),N,[],'omitnan','truncate'); % smoothed res
        elseif j == 2
            data.ModIT2(ind) = mod.finMod; % retransform
            data.ResIT2(ind) = data.MYRF(ind) - data.ModIT2(ind);        % calculate residuals
            data.ResPIT2(ind) = (data.ResIT2(ind)./data.ModIT2(ind)).*100; % calculate % residuals
            data.SmoRESPIT2(ind) = medfilt1(data.ResPIT2(ind),N,[],'omitnan','truncate'); % smoothed res
        elseif j == 3
            data.ModIT3(ind) = mod.finMod; % retransform
            data.ResIT3(ind) = data.MYLH(ind) - data.ModIT3(ind);        % calculate residuals
            data.ResPIT3(ind) = (data.ResIT3(ind)./data.ModIT3(ind)).*100; % calculate % residuals
            data.SmoRESPIT3(ind) = medfilt1(data.ResPIT3(ind),N,[],'omitnan','truncate'); % smoothed res
        else
            data.ModIT4(ind) = mod.finMod; % retransform
            data.ResIT4(ind) = data.MYRH(ind) - data.ModIT4(ind);        % calculate residuals
            data.ResPIT4(ind) = (data.ResIT4(ind)./data.ModIT4(ind)).*100; % calculate % residuals
            data.SmoRESPIT4(ind) = medfilt1(data.ResPIT4(ind),N,[],'omitnan','truncate'); % smoothed res
        end
    end

    if makefig == 1
        saveas(h,[savename 'mod_ex_' num2str(i) '.fig'])
        saveas(h,[savename 'mod_ex_' num2str(i) '.tif'])
    end
end
clear i DIM ind idx MY OL prevMI T j MI 


% plot percentual
if makefig ==1
    close all
    for i = 1:height(cowlac)
        % find data of cow
        ind = find(data.B_ID == cowlac.B_ID(i) & data.Lac == cowlac.Lac(i));
        h = figure(i+120); set(gcf,'units','normalized','outerposition',[0 0 1 1]);
        T1 = find(startsWith(data.Properties.VariableNames,'ModIT1'));  % index of first
        T2 = find(startsWith(data.Properties.VariableNames,'MYLF'));  % index of first
        GO = 0;
        for j = 1:4
            if sum(data{ind,T2+j-1}) > 10 && length(find(data{ind,T2+j-1} >0))> 30 
                % prepare plot data
                subplot(4,2,2*j-1); hold on; box on; grid on; 
                xlabel('DIM (days)'); ylabel('MY (kg)');
                title(['i = ' num2str(i) ...
                       ', CowID = ' num2str(cowlac.B_ID(i)), ...
                       ', Lac = ' num2str(cowlac.Lac(i))])

                % plot the data and the model
                plot(data.DIM(ind),data{ind,T2+(j-1)},'o',pl.LW,1.2,pl.MS,3,pl.C,[0.8 0.4 0.6]);
                plot(data.DIM(ind),data{ind,T1+(4*j)-4},'-',pl.LW,1,pl.MS,3,pl.C,[0.2 0 0.6]);
                axis tight; ylim([0 max(data{ind,T2+(j-1)})+1])

                % prepare plot residuals
                subplot(4,2,2*j); hold on; box on; grid on
                xlabel('DIM (days)'); ylabel('residual MY (%)');
                title(['i = ' num2str(i) ...
                       ', CowID = ' num2str(cowlac.B_ID(i)), ...
                       ', Lac = ' num2str(cowlac.Lac(i))])

                % plot the residuals procentual + smoothed
                plot(data.DIM(ind),0*ones(length(ind),1),'k-',pl.LW,1.4)
                plot(data.DIM(ind),data{ind,T1+(4*j)-2},'o',...
                     pl.LW,1.2,pl.MS,3,pl.C,[0.8 0.4 0.6]);
                plot(data.DIM(ind),data{ind,T1+(4*j)-1},'-',...
                     pl.LW,1.5,pl.MS,3,pl.C,[0.2 0 0.6]);
                plot(data.DIM(ind),-20*ones(length(ind),1),'-',pl.LW,1.9,pl.C,[0 0.8 0.6])
                axis tight; ylim([-100 50])
                GO =1;
            end
        end
        if GO == 1
            saveas(h,[savename 'pert_ex_' num2str(i) '.fig'])
            saveas(h,[savename 'pert_ex_' num2str(i) '.tif'])
        end
    end
end

clear ans F i ind j mSel N T1 T2 h GO pl makefig savename
close all

% shape outputs
data = removevars(data,{'ResIT1','ResPIT1','SmoRESPIT1',...
                        'ResIT2','ResPIT2','SmoRESPIT2',...
                        'ResIT3','ResPIT3','SmoRESPIT3',...
                        'ResIT4','ResPIT4','SmoRESPIT4'});
T = find(startsWith(data.Properties.VariableNames,'MYLF'));  % index MY
data = data(~isnan(sum(data{:,T:T+3},2)),:);  % delete NaN
% prepare out tables
mod.pFin = array2table(mod.pFin,'VariableNames',{'B_ID','Lac',...
                'Q1a','Q1b','Q1c','Q2a','Q2b','Q2c',...
                'Q3a','Q3b','Q3c','Q4a','Q4b','Q4c'});
try
    mod = rmfield(mod,{'opts','UB','LB','my','dim','p0','p','yval','std',...
                   'rmse','pIni','dimIT','myIT','stdIT','yvalIT','pIT'...
                   'rmseIT','n','finMod','Wood'});
catch   % if no QMY data available
    mod = rmfield(mod,{'opts','UB','LB',...
                   'pIT',...
                   'finMod','Wood'});
    
end
mod.nFin = array2table(mod.nFin,'VariableNames',{'B_ID','Lac','Q1n','Q2n','Q3n','Q4n'});
mod.rmseFin = array2table(mod.rmseFin,'VariableNames',{'B_ID','Lac','Q1rmse','Q2rmse','Q3rmse','Q4rmse'});