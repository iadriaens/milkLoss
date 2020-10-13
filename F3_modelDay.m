function [data,mod] = F3_modelDay(data,makefig,savename)
% this function models the quarter milk yields and estimates the
% unperturbed curve at milking level
% INPUTS:   1   data:    contains MY at DAY level, DIM, B_ID, Lac
%           2   makefig: whether or not figures are to be made / output
%           3   savedir: directory to save
%
% OUTPUTS:  1   data:    original data with additional columns added
%                        cols added = DAY & MOD
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
% set file name to read data
% % % fileName = ['C:\Users\u0084712\Documents\MastiMan\Research\dataMining\MILK LOSSES scripts\A1_results\' ...
% % %              'DAY_Kestier_20150307_20190225.txt'];
% % % %            'MILK LOSSES scripts\A1_results\MILK_Sanders_20161106_20200409.txt'];
% % % cowsel = 1:10;%[23 24 25 26 28 70 71 104 105 111 112];
% % %      
% % % % read file
% % % opts = detectImportOptions(fileName);
% % % data = readtable(fileName,opts);
% % % cowlac = unique(data(:,[2 6]),'rows');
% % % cowlac = cowlac(cowsel,:); % select a number interesting lactations to train
% % % 
% % % % select data of selected cows
% % % data = innerjoin(data,cowlac,'Keys',{'B_ID','Lac'}); % select
% % % data = data(data.DIM < 306,:);  % trim lactations up to 305 days
% % % clear opts cowlac cowsel fileName
% % % 
% % % % makefig prepare + savename
% % % makefig = 1;
% % % savename = 'C:\Users\u0084712\Documents\MastiMan\Research\mastitisMilkLoss\A1_results\Kestier_';
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
mod.UB = [70 5 0.01];                               % set upper boundary
mod.LB = [0 0 0];                                   % set lower boundary



%% STEP 1: INDIVIDUAL MODEL FIT/RESIDUALS -ITERATIVE
% summarize cow lactations based on B_ID and Lac vars
cowlac = array2table(unique([data.B_ID,data.Lac],'rows'),'VariableNames',{'B_ID','Lac'});

% modelling loop
for i = 1:height(cowlac)        % for all cow lactations entered
    % find data of the specific lactation
    T = find(startsWith(data.Properties.VariableNames,'TDMY'));  % index MY
    ind = find(data.B_ID == cowlac.B_ID(i) & data.Lac == cowlac.Lac(i) & ~isnan(sum(data{:,T},2)));
    
    % determine figure making >> overrul input
    if sum(ismember(i,[22 44 66 8 10]))>0
        makefig = 0;
    else
        makefig = 0;
    end
    
    % prepare figure of model
    if makefig == 1
        h=figure('Units','Normalized','Outerposition',[0 0.5 1 0.5]);
    end
    
    if nansum(data{ind,T}) > 10 && length(find(data{ind,T} >0))> 30 % at least 30 measurements and >10kg over whole lac
        % data selection
        DIM = data.DIM(ind);           % days in milk
        MY = (data{ind,T}); % milk yield / hour * 24
        
        % check the variability in the data
        if data.RobotType == 1%var(MY) > 17 % this is DELAVAL DATA
            mod.my = smooth(MY,3);
        else
            mod.my = MY;
        end
        mod.dim = DIM;
        
        % visualisation of point selection
        if makefig == 1
            hold on; box on; axis tight;
            xlabel('DIM (days)'); ylabel('MY (kg)')
            title(['i = ' num2str(i) ...
                ', CowID = ' num2str(cowlac.B_ID(i)) ...
                ', Lac = ' num2str(cowlac.Lac(i)) ...
                ])
            plot(mod.dim,mod.my(:,1),'o:',pl.LW,1,pl.MS,3,pl.C,[0 0.80 0.8 1])
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
                ylim([0 max(MY(:,1)+5)]);
            end
            
            % final model
            mod.finMod = mod.Wood(mod.pIT,data.DIM(ind)); % retransform for final model
            mod.rmseFin(i,3) = mod.rmse;        % store final rmse
            mod.nFin(i,3) = mod.n;                % store final n
        
        else
            mod.pIT = [NaN NaN NaN];            % no values
            mod.finMod = NaN*zeros(length(ind),1); % no values
            mod.rmseFin(i,3) = NaN;           % store final rmse
            mod.nFin(i,3) = NaN;              % store n iterations
        end
    else
        mod.pIT = [NaN NaN NaN];            % no values
        mod.finMod = NaN*zeros(length(ind),1); % no values
        mod.rmseFin(i,3) = NaN;           % store final rmse
        mod.nFin(i,3) = NaN;              % store n iterations
    end
    
    % final pars
    mod.pFin(i,[1 2]) = cowlac{i,:};        % store cowlac
    mod.pFin(i,3:5) = mod.pIT;    % store final p
    mod.rmseFin(i,[1 2]) = cowlac{i,:};     % store cowlac rmse
    mod.nFin(i,[1 2]) = cowlac{i,:};        % store cowlac n

    
% % %     % prepare smoothing
% % %     N = 7; % window of % residual smoothing
    
    % store vars
    data.MYs(ind) = mod.my;
    data.ModIT(ind) = mod.finMod; % final model
        
    if makefig == 1
        saveas(h,[savename 'DAY_mod_ex_' num2str(i) '.fig'])
        saveas(h,[savename 'DAY_mod_ex_' num2str(i) '.tif'])
    end
end
clear i DIM ind idx MY OL prevMI T j MI 


% prepare out tables
mod.pFin = array2table(mod.pFin,'VariableNames',{'B_ID','Lac',...
                'a','b','c'});
try
    mod = rmfield(mod,{'opts','UB','LB','my','dim','p0','p','yval','std',...
                   'rmse','pIni','dimIT','myIT','stdIT','yvalIT','pIT'...
                   'rmseIT','n','finMod','Wood'});
catch   % if no QMY data available
    mod = rmfield(mod,{'opts','UB','LB',...
                   'pIT',...
                   'finMod','Wood'});
    
end
mod.nFin = array2table(mod.nFin,'VariableNames',{'B_ID','Lac','n'});
mod.rmseFin = array2table(mod.rmseFin,'VariableNames',{'B_ID','Lac','rmse'});