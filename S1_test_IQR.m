%% S1_test_IQR

clear variables 
close all
clc

% set file name to read data
fileName = ['C:\Users\u0084712\Documents\MastiMan\Research\dataMining\' ...
           'MILK LOSSES scripts\A1_results\MILK_Sanders_20161106_20200409.txt'];
cowsel = [23 24 25 26 28 70 71 104 105 111 112];
     
% read file
opts = detectImportOptions(fileName);
data = readtable(fileName,opts);
cowlac = unique(data(:,[2 6]),'rows');
cowlac = cowlac(cowsel,:); % select a number interesting lactations to train

% select data of selected cows
data = innerjoin(data,cowlac,'Keys',{'B_ID','Lac'}); % select
data = data(data.DIM < 306,:);  % trim lactations up to 305 days

% prepare plotting
pl.LW = 'LineWidth';
pl.MS = 'MarkerSize';
pl.C = 'Color';
pl.MFC = 'MarkerFaceColor';



%% INDIVIDUAL MODEL FIT/RESIDUALS -ITERATIVE
% prepare Wood modelling
mod.Wood = @(p,t) p(1).*t.^p(2).*exp(-p(3).*t);     % wood model - prepare
mod.opts = optimset('Display','off', 'MaxIter',1000);% set display to off
mod.UB = [25 1 0.01];                               % set upper boundary
mod.LB = [0 0 0];                                   % set lower boundary
F =0;
data.MI(data.MI>30) = NaN;                          % set MI > 30 to NaN
close all
for i = 1:height(cowlac)
    prevMI = [5; data.MI(1:end-1)]; % array with previous milking intervals
    ind = find(data.B_ID == cowlac.B_ID(i) & data.Lac == cowlac.Lac(i));
    T = find(startsWith(data.Properties.VariableNames,'MYLF'));  % index
    
    h=figure('Units','Normalized','Outerposition',[0 0 1 1]);
    for j = 1:4
        if length( find(data{ind,T+j-1} >0))> 30 % at least 30 measurements
            % data selection
            DIM = data.DIM(ind);           % days in milk
            MI = min(20,data.MI(ind)); MI = max(MI,3);  % adjust extreme MI
            MY = (data{ind,T+j-1}./MI)*24; % milk yield / hour * 24
            
%           =============> OL detect via moving window <============
            [OL.LT, OL.low, OL.up,OL.cent] = isoutlier(MY(:,1),'movmedian',36,'thresholdfactor',4); % 1 week MM = 18meas
            OL.mad = movmad(MY(:,1)-OL.cent,36);
            OL.up = OL.cent + OL.mad*4;
%           =============> OL detect via moving window <============
            
            
            MY(MY<OL.up,2) = 1;            % all measurement below the upper threshold
            mod.my = MY(MY(:,2) ==1,1);    % select measurments
            mod.dim = DIM(MY(:,2) == 1);   % select dim
            
            % visualisation of point selection
            subplot(2,2,j); hold on; box on; axis tight;
            xlabel('DIM (days)'); ylabel('MY/h*24 (kg)')
            title(['i = ' num2str(i) ...
                ', CowID = ' num2str(cowlac.B_ID(i)) ...
                ', Lac = ' num2str(cowlac.Lac(i)) ...
                ', Q = ' data.Properties.VariableNames{T+j-1}])
            plot(DIM,MY(:,1),'o:',pl.LW,1,pl.MS,3,pl.C,[0 0.80 0.8 1])
            plot(DIM, OL.up,pl.LW,1.2,pl.C, [1 0.4 0.2])
            plot(DIM, OL.cent ,pl.LW,1.2,pl.C, [0.8 0 0.6])
            plot(DIM(MY(:,2)==0), MY(MY(:,2)==0,1),'ks',pl.LW,2,pl.MS,3.5)
            
            % model preparation
            mod.p0 = [nanmean(mod.my) 0.20 0.006];   % set initial parameters
            mod.p = lsqcurvefit(mod.Wood,mod.p0,mod.dim,mod.my,mod.LB,mod.UB,mod.opts); % estimate p
            mod.yval = mod.Wood(mod.p,mod.dim);      % y value initial wood fit
            
            % calculate std using all measurements after DIM=5
            mod.std = std(mod.my(mod.dim>=5) - mod.yval(mod.dim>=5)); % std of the residuals
            idx = find(mod.my-(mod.yval-1.6.*mod.std)<0);  % this are the INITIAL measurements outside the range of 1.6 STD
            mod.rmse = sqrt(mean((mod.my-mod.yval).^2));   % calculate rmse
            mod.pIni(i,:) = [cowlac{i,:} mod.p mod.rmse];  % store initial p and rmse
            
            % visualisation initalisation modeling
            plot(mod.dim,mod.Wood(mod.p0,mod.dim),'b--',pl.LW,1.2);
            plot(mod.dim,mod.Wood(mod.p,mod.dim),'r--',pl.LW,2);
            plot(mod.dim(idx),mod.my(idx),'s','LineWidth',2,'MarkerSize',3.5, 'Color',[177/255 0 0 ]);
                        
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
                plot(DIM,mod.Wood(mod.pIT,DIM),'--',pl.LW,1.5,pl.C,[0.6 0 0.8])
                plot(mod.dimIT(mod.myIT <= (mod.yvalIT-1.6*mod.stdIT) & mod.dimIT > 5),...
                     mod.myIT(mod.myIT <= (mod.yvalIT-1.6*mod.stdIT) & mod.dimIT>5),...
                     'o',pl.LW,2,pl.MS,3.5,pl.C,[0.8 0 0.6],...
                     pl.MFC,[0.8 0 0.6]);
            end
            
            % plot final model
            plot(DIM, mod.Wood(mod.pIT,DIM),'-',pl.LW,2.5,pl.C,[1 0 1])
            ylim([0 max(MY(:,1)+0.5)]);
            
            % final model
            mod.finMod = (mod.Wood(mod.pIT,data.DIM(ind))./24).*MI; % retransform for final model
        else
            mod.pIT = [NaN NaN NaN];  % no values
            mod.finMod = NaN*zeros(length(ind),1); % no values
        end
        
        % final pars
        mod.pFin(i,[1 2]) = cowlac{i,:};  % store initial p and rmse
        mod.pFin(i,2+3*j-2:2+3*j) = mod.pIT;          % store final p
        
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

%     saveas(h,['C:\Users\u0084712\Documents\MastiMan\Research\mastitisMilkLoss\A1_results\Mod_ex_' num2str(i) '_Vandevloet.fig'])
%     saveas(h,['C:\Users\u0084712\Documents\MastiMan\Research\mastitisMilkLoss\A1_results\Mod_ex_' num2str(i) '_Vandevloet.tif'])
end
clear i DIM ind idx MY OL prevMI T j MI 

% plot percentual 
% close all
for i = 1:height(cowlac)
    % find data of cow
    ind = find(data.B_ID == cowlac.B_ID(i) & data.Lac == cowlac.Lac(i));
    h = figure(i+120); set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    T1 = find(startsWith(data.Properties.VariableNames,'ModIT1'));  % index of first
    T2 = find(startsWith(data.Properties.VariableNames,'MYLF'));  % index of first

    for j = 1:4
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
    end
%     saveas(h,['C:\Users\u0084712\Documents\MastiMan\Research\mastitisMilkLoss\A1_results\MPert_ex' num2str(i) '_Vandevloet.fig'])
%     saveas(h,['C:\Users\u0084712\Documents\MastiMan\Research\mastitisMilkLoss\A1_results\MPert_ex' num2str(i) '_Vandevloet.tif'])
end

% visualize pars
figure('Units','normalized','OuterPosition',[0 0.5 0.8 0.5]); 
subplot(1,3,1); 
    histogram([mod.pFin(:,3); mod.pFin(:,6); mod.pFin(:,9); mod.pFin(:,12)])
    
    title('Wood: a')
subplot(1,3,2); 
    histogram([mod.pFin(:,4); mod.pFin(:,7); mod.pFin(:,10); mod.pFin(:,13)])
    title('Wood: b')
subplot(1,3,3); 
    histogram([mod.pFin(:,5); mod.pFin(:,8); mod.pFin(:,11); mod.pFin(:,14)])
    title('Wood: c')

clear ans F i ind j mSel N T1 T2 h 
% close all

%% Test Case selection
% criteria: 


% 


















%% SMOOTH with LOWESS

for i = 1:height(cowlac)
    ind = find(data.B_ID == cowlac.B_ID(i) & data.Lac == cowlac.Lac(i));
    
    % calculate the iqr (not per hour)
    %     figure; hold on;
    data.IQsmooth1(ind) = smooth(data.MYLF(ind),'lowess')./smooth(nanmean([data.MYRF(ind),data.MYLH(ind),data.MYRH(ind)],2),'lowess');
    %      plot(data.DIM(1:800), data.iqMY1(1:800),'o:',pl.LW,1.5,pl.MS,4);
    data.IQsmooth2(ind) = smooth(data.MYRF(ind),'lowess')./smooth(nanmean([data.MYLF(ind),data.MYLH(ind),data.MYRH(ind)],2),'lowess');
    data.IQsmooth3(ind) = smooth(data.MYLH(ind),'lowess')./smooth(nanmean([data.MYLF(ind),data.MYRF(ind),data.MYRH(ind)],2),'lowess');
    data.IQsmooth4(ind) = smooth(data.MYRH(ind),'lowess')./smooth(nanmean([data.MYLF(ind),data.MYRF(ind),data.MYLH(ind)],2),'lowess');
end
% plot
close all
for i = 1:height(cowlac)
    ind = find(data.B_ID == cowlac.B_ID(i) & data.Lac == cowlac.Lac(i));
    h = figure(floor((i-0.01)/4)+1); set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    subplot(2,2,ceil(rem((i-0.01),4))); hold on; box on;
    xlabel('DIM (days)'); ylabel('data');
    title(['i = ' num2str(i) ...
        ', CowID = ' num2str(cowlac.B_ID(i)), ...
        ', Lac = ' num2str(cowlac.Lac(i))])
%     plot(data.DIM(ind),data.MYLF(ind)./data.MI(ind),'o:',...
%         pl.LW,1.2,pl.MS,3,pl.C,[0.6 0.4 0.6]);
%     plot(data.DIM(ind),data.MYRF(ind)./data.MI(ind),'o:',...
%         pl.LW,1.2,pl.MS,3,pl.C,[0.6 0.6 0.6]);
%     plot(data.DIM(ind),data.MYLH(ind)./data.MI(ind),'o:',...
%         pl.LW,1.2,pl.MS,3,pl.C,[0.6 0.8 0.6]);
%     plot(data.DIM(ind),data.MYRH(ind)./data.MI(ind),'o:',...
%         pl.LW,1.2,pl.MS,3,pl.C,[0.6 1 0.6]);
%     yyaxis right
    plot(data.DIM(ind),data.IQsmooth1(ind),'o:',...
        pl.LW,1.2,pl.MS,3,pl.C,[0.6 0 0.6]);
    plot(data.DIM(ind),data.IQsmooth2(ind),'o:',...
        pl.LW,1.2,pl.MS,3,pl.C,[0.6 0.1 0.4]);
    plot(data.DIM(ind),data.IQsmooth3(ind),'o:',...
        pl.LW,1.2,pl.MS,3,pl.C,[0.2 0.6 0.2]);
    plot(data.DIM(ind),data.IQsmooth4(ind),'o:',...
        pl.LW,1.2,pl.MS,3,pl.C,[0.4 0.8 0.2]);
    axis tight; ylim([0 2.5])
end
    
%% MEDIAN SMOOTH

for i = 1:height(cowlac)
    ind = find(data.B_ID == cowlac.B_ID(i) & data.Lac == cowlac.Lac(i));
    N = 15; % number of measurments
    % calculate the iqr (not per hour)
    %     figure; hold on;
    data.IQfilt1(ind) = medfilt1(data.MYLF(ind),N,[],'omitnan','truncate')./medfilt1(nanmean([data.MYRF(ind),data.MYLH(ind),data.MYRH(ind)],2),N,[],'omitnan','truncate');
    %      plot(data.DIM(1:800), data.iqMY1(1:800),'o:',pl.LW,1.5,pl.MS,4);
    data.IQfilt2(ind) = medfilt1(data.MYRF(ind),N,[],'omitnan','truncate')./medfilt1(nanmean([data.MYLF(ind),data.MYLH(ind),data.MYRH(ind)],2),N,[],'omitnan','truncate');
    data.IQfilt3(ind) = medfilt1(data.MYLH(ind),N,[],'omitnan','truncate')./medfilt1(nanmean([data.MYLF(ind),data.MYRF(ind),data.MYRH(ind)],2),N,[],'omitnan','truncate');
    data.IQfilt4(ind) = medfilt1(data.MYRH(ind),N,[],'omitnan','truncate')./medfilt1(nanmean([data.MYLF(ind),data.MYRF(ind),data.MYLH(ind)],2),N,[],'omitnan','truncate');
end
% plot
close all
for i = 1:height(cowlac)
    ind = find(data.B_ID == cowlac.B_ID(i) & data.Lac == cowlac.Lac(i));
    h = figure(floor((i-0.01)/4)+4); set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    subplot(2,2,ceil(rem((i-0.01),4))); hold on; box on;
    xlabel('DIM (days)'); ylabel('data');
    title(['i = ' num2str(i) ...
        ', CowID = ' num2str(cowlac.B_ID(i)), ...
        ', Lac = ' num2str(cowlac.Lac(i))])
%     plot(data.DIM(ind),data.MYLF(ind)./data.MI(ind),'o:',...
%         pl.LW,1.2,pl.MS,3,pl.C,[0.6 0 0.6]);
%     plot(data.DIM(ind),data.MYRF(ind)./data.MI(ind),'o:',...
%         pl.LW,1.2,pl.MS,3,pl.C,[0.6 0.1 0.4]);
%     plot(data.DIM(ind),data.MYLH(ind)./data.MI(ind),'o:',...
%         pl.LW,1.2,pl.MS,3,pl.C,[0.2 0.6 0.2]);
%     plot(data.DIM(ind),data.MYRH(ind)./data.MI(ind),'o:',...
%         pl.LW,1.2,pl.MS,3,pl.C,[0.4 0.8 0.2]);
%      axis tight; ylim([0 1.2])
%     yyaxis right;
    plot(data.DIM(ind),data.IQfilt1(ind),'o:',...
        pl.LW,1.2,pl.MS,3,pl.C,[0.6 0 0.6]);
    plot(data.DIM(ind),data.IQfilt2(ind),'o:',...
        pl.LW,1.2,pl.MS,3,pl.C,[0.6 0.1 0.4]);
    plot(data.DIM(ind),data.IQfilt3(ind),'o:',...
        pl.LW,1.2,pl.MS,3,pl.C,[0.2 0.6 0.2]);
    plot(data.DIM(ind),data.IQfilt4(ind),'o:',...
        pl.LW,1.2,pl.MS,3,pl.C,[0.4 0.8 0.2]);
    axis tight; %ylim([0 2.5])
end
































































