% script to run simulations of control model with and without vaccination,
% produces Figure 2 and Figure S5 (strategies without/with vaccination) from
% the manuscript

% NB: This file is analagous to Simulations.m with a difference in subpanel
% formatting. To save thresholds to mat file, run Simulations.m
clear; close all

%Plotting preferences
set(0,'defaultlinelinewidth',2)
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultTextInterpreter','latex')
set(0,'defaultaxesfontsize',16)

%%%%%%%%%%%%%%%%%  INPUT FOR SUBPANELS  %%%%%%%%%%%%%%%%%

%Fig 2:  vstarts = 2160 (or anything > 1080)
%Fig S5: vstarts = 360
%FigXa:  which_strat = 1 (Cautious Easing)
%FigXb:  which_strat = 2 (Suppression)
%FigXc:  which_strat = 3 (Slow Control)
%FigXd:  which_strat = 4 (Rapid Control)

vstarts = 2160;
which_strat = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% MAIN SCRIPT

% load default parameters
if not(isfolder('mats'))
    disp('No parameters saved: Running define_params.m')
    define_params
end

para0 = load('./mats/Parameters.mat');

% Define time to run model for
t_init = 30;    % preliminary run
maxtime = 720;  % main simulation

% define strategy numbers and switching thresholds
if not(exist('mats/Thresholds.mat','file'))
    disp('No strategy thresholds saved: Running Simulations.m')
    Simulations
    close all;
end

load('mats/Thresholds.mat')
strategies = 1:length(thresholds);

% which strategy to reproduce
SELECTED_STRAT = strategies(which_strat);

% plotting preperation
figure('Position',[200 400 500 250])
stIDs = {'S1', 'S2', 'S3', 'S4'};
stnames = {'(Cautious easing)', '(Suppression)', '(Slow control)', '(Rapid control)'};
stpos = [25 175 325 500; 25 175 325 475; 75 225 375 525; 100 250 400 550];

dt = t_init/3;
markedtimes = dt:dt:(maxtime-dt);
mygreen = [0.1 0.5 0.1];
myblue = [0 0.2 0.8];

tic
for strat = SELECTED_STRAT
    % Preliminary run - no control, 30 day build-up
    para0.vstart = vstarts;
    [Prelim, Prelim_ICs] = Get_ICs(para0);
        
    % add control thresholds defined by strategy
    para = para0;
    para.maxtime = maxtime;
    para.T10 = thresholds(strat,1);
    para.T01 = thresholds(strat,2);
    para.T21 = thresholds(strat,3);
    para.T12 = thresholds(strat,4);

    % starting control state
    if sum(Prelim.IH(end,:)) < para.T12
        para.init = 1;
    else
        para.init = 2;
    end

    % Run model
    [Classes, burden, stringency, peak_hospital] = ODEmodel(para,Prelim_ICs);     
    peak_hospital
    round(sum(burden))
    sum(stringency)

    % Post-Processing ixs, SD Class (only needed for plotting)
    % Find times where control measures are enforced
    ix1 = find(Classes.SD(:,2)==1);
    ix2 = find(Classes.SD(:,2)==2);

    % append SD for lockdown computation if we end in a restriction
    if Classes.SD(end,2) ~= 0
        Classes.SD(end+1:end+2,:) = [Classes.t(end) 0; Classes.t(end) 0];
    end


    % plotting active hospitalisations
    %yyaxis left
    ax1 = gca;
    ax1.YColor = 'k';
    ax1.FontSize = 16;
    ax1.FontSizeMode = 'manual';

    for i = ix1'
        patch([Classes.SD(i,1) Classes.SD(i,1) Classes.SD(i+2,1) Classes.SD(i+2,1)], [0 20000 20000 0], 'y', 'Facealpha',0.25, 'EdgeAlpha',0)
        hold on
    end
    for i = ix2'
        patch([Classes.SD(i,1) Classes.SD(i,1) Classes.SD(i+2,1) Classes.SD(i+2,1)], [0 20000 20000 0], 'r', 'Facealpha',0.25, 'EdgeAlpha',0)
        hold on
    end

    if vstarts < maxtime
        if strat == 1
            xline(para.vstart,'-',{'Vaccine','Arrival'},'Color',mygreen,'Linewidth',2,'FontSize',18,'Interpreter','latex','LabelOrientation','horizontal','LabelHorizontalAlignment','left')
        else
            xline(para.vstart,'-','Color',mygreen,'Linewidth',2)
        end
    end

%     plot(markedtimes, para.T01.*ones(size(markedtimes)), 'Marker', '_', 'MarkerSize', 2, 'MarkerEdgeColor', 'b', 'LineStyle', 'none')
%     hold on
%     plot(markedtimes, para.T12.*ones(size(markedtimes)), 'Marker', '_', 'MarkerSize', 2, 'MarkerEdgeColor', 'b', 'LineStyle', 'none')
%     hold on
%     plot(markedtimes, para.T10.*ones(size(markedtimes)), 'Marker', 'o', 'MarkerSize', 2, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'LineStyle', 'none')
%     hold on
%     plot(markedtimes, para.T21.*ones(size(markedtimes)), 'Marker', 'o', 'MarkerSize', 2, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'LineStyle', 'none')
%     hold on
    T1 = plot(markedtimes, para.T01.*ones(size(markedtimes)), '--', 'Color', myblue, 'LineWidth', 1.5);
    hold on
    T2 = plot(markedtimes, para.T12.*ones(size(markedtimes)), '--', 'Color', myblue, 'LineWidth', 1.5);
    hold on
    T3 = plot(markedtimes, para.T10.*ones(size(markedtimes)), '.', 'Color', myblue, 'LineWidth', 1.5);
    hold on
    T4 = plot(markedtimes, para.T21.*ones(size(markedtimes)), '.', 'Color', myblue, 'LineWidth', 1.5);
    hold on

    plot(Classes.t, sum(Classes.IH,2), 'k', 'LineWidth', 3.5)

    if strat > 1/2
        xlabel('Time (days)')
    end

    ax1.Position(1) = ax1.Position(1) + 0.07;
    ax1.Position(3) = ax1.Position(3) - 0.05;

    axis([0 maxtime 0 1100])

    ylab = ylabel('$I^H(t)$','Rotation',0);
    ylab.Position(1) = ylab.Position(1) - 25;


%     yyaxis right
%     ax2 = gca;
%     ax2.YColor = 'k';
%     ax2.TickDir = 'none';
%     ax2Y = ax2.YAxis(2,1);
%     ax2Y.FontSize = 14;
%     ax2.FontSizeMode = 'manual';
% 
% 
% 
%     plot(markedtimes, para.T01.*ones(size(markedtimes)), 'Marker', '_', 'MarkerSize', 5, 'MarkerEdgeColor', mygrey, 'LineStyle', 'none')
%     hold on
%     plot(markedtimes, para.T12.*ones(size(markedtimes)), 'Marker', 'o', 'MarkerSize', 5, 'MarkerEdgeColor', mygrey, 'LineStyle', 'none')
%     hold on
%     plot(markedtimes, para.T10.*ones(size(markedtimes)), 'Marker', 'x', 'MarkerSize', 5, 'MarkerEdgeColor', mygrey, 'LineStyle', 'none')
%     hold on
%     plot(markedtimes, para.T21.*ones(size(markedtimes)), 'Marker', 's', 'MarkerSize', 5, 'MarkerEdgeColor', mygrey, 'LineStyle', 'none')
%     yticks(stpos(strat,:))
%     
%     if para.T21 > para.T01
%         yticklabels({'$T_{10}$','$T_{01}$','$T_{21}$','$T_{12}$'})
%     else
%         yticklabels({'$T_{10}$','$T_{21}$','$T_{01}$','$T_{12}$'})
%     end

    xticks(0:180:maxtime)
    xtickangle(0)
    axis([0 maxtime 0 1100])
            
    grid on
 
end
toc

%save figure
if not(isfolder('./figs/sim_images'))
    mkdir('./figs/sim_images')
end

saveas(gcf,strcat('./figs/sim_images/','simulation_',num2str(para.vstart),'_strat',num2str(strat),'.png'))

