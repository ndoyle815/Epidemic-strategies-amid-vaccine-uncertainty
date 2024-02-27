% script to run simulations of control model with and without vaccination,
% produces Figure 2 and Figures S4-5 (strategies without/with vaccination) from
% the manuscript
clear

%Plotting preferences
set(0,'defaultlinelinewidth',2)
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultTextInterpreter','latex')
set(0,'defaultaxesfontsize',16)

% load default parameters
if not(exist('mats/Parameters.mat','file'))
    disp('No parameters saved: Running define_params.m')
    define_params
end

para0 = load('./mats/Parameters.mat');

% vaccination
vstarts = 2160;

% Define time to run model for
t_init = 30;    % preliminary run
maxtime = 720;  % main simulation

% define strategy numbers and switching thresholds
thresholds = [50 100 150 500; 50 100 150 200; 100 200 400 500; 200 300 350 450];
%load('./mats/Thresholds.mat')
strategies = 1:length(thresholds);

% plotting preperation
plot_cumulative = 1;
figlength = [450 900];
figure('Position',[200 400 figlength(plot_cumulative) 900])
stIDs = {'S1', 'S2', 'S3', 'S4'};
stnames = {'(Cautious easing)', '(Suppression)', '(Slow control)', '(Rapid control)'};
stpos = [25 175 325 500; 25 175 325 475; 75 225 375 525; 100 250 400 550];

tic
for strat = strategies
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
    %sum(Classes.Cases(end,:),2)/sum(Classes.Hosp(end,:),2)

    % Post-Processing ixs, SD Class (only needed for plotting)
    % Find times where control measures are enforced
    ix1 = find(Classes.SD(:,2)==1);
    ix2 = find(Classes.SD(:,2)==2);

    % append SD for lockdown computation if we end in a restriction
    if Classes.SD(end,2) ~= 0
        Classes.SD(end+1:end+2,:) = [Classes.t(end) 0; Classes.t(end) 0];
    end

    plotidx = plot_cumulative*(strat-1) + 1;  % index for subfigure

    % plotting active hospitalisations
    subplot(length(strategies),plot_cumulative,plotidx)

    yyaxis left
    ax1 = gca;
    ax1.YColor = 'k';
    ax1.FontSize = 16;
    ax1.FontSizeMode = 'manual';

    for i = ix1'
        patch([Classes.SD(i,1) Classes.SD(i,1) Classes.SD(i+2,1) Classes.SD(i+2,1)], [0 20000 20000 0], 'y', 'Facealpha',0.3, 'EdgeAlpha',0)
        hold on
    end
    for i = ix2'
        patch([Classes.SD(i,1) Classes.SD(i,1) Classes.SD(i+2,1) Classes.SD(i+2,1)], [0 20000 20000 0], 'r', 'Facealpha',0.3, 'EdgeAlpha',0)
        hold on
    end
    plot(Classes.t, sum(Classes.IH,2), 'k', 'LineWidth', 2.5)
    
    if vstarts < maxtime
        if strat == 1
            xline(para.vstart,'-',{'Vaccine','Arrival'},'Color','b','Linewidth',2,'FontSize',16,'Interpreter','latex','LabelOrientation','horizontal','LabelHorizontalAlignment','left')
        else
            xline(para.vstart,'-','Color','b','Linewidth',2)
        end
    end

    if strat == strategies(end)
        xlabel('Time (days)')
    end

    if plot_cumulative == 2
        %ylab = ylabel({char(stIDs(strat)) ; char(stnames(strat))});
    end

    if strat == 1
        title('Active $I^H(t)$');
    end

    axis([0 maxtime 0 1100])

    yyaxis right
    ax2 = gca;
    ax2.YColor = 'k';
    ax2.TickDir = 'none';
    ax2Y = ax2.YAxis(2,1);
    ax2Y.FontSize = 14;
    ax2.FontSizeMode = 'manual';

    plot(Classes.t, para.T01.*ones(size(Classes.t)), 'k--', 'LineWidth',0.5)
    hold on
    plot(Classes.t, para.T12.*ones(size(Classes.t)), 'k--', 'LineWidth',0.5)
    hold on
    plot(Classes.t, para.T10.*ones(size(Classes.t)), 'k--', 'LineWidth',0.5)
    hold on
    plot(Classes.t, para.T21.*ones(size(Classes.t)), 'k--', 'LineWidth',0.5)
    yticks(stpos(strat,:))
    
    if para.T21 > para.T01
        yticklabels({'$T_{10}$','$T_{01}$','$T_{21}$','$T_{12}$'})
    else
        yticklabels({'$T_{10}$','$T_{21}$','$T_{01}$','$T_{12}$'})
    end

    xticks(0:180:maxtime)
    xtickangle(0)
    axis([0 maxtime 0 1100])
            
    grid on
 

    % plotting cumulative hospitalisations
    if plot_cumulative == 2
        plotidx = plotidx + 1;
        subplot(length(strategies),plot_cumulative,plotidx)
    
        for i = ix1'
            patch([Classes.SD(i,1) Classes.SD(i,1) Classes.SD(i+2,1) Classes.SD(i+2,1)], [0 40000 40000 0], 'y', 'Facealpha',0.3, 'EdgeAlpha',0)
            hold on
        end
        for i = ix2'
            patch([Classes.SD(i,1) Classes.SD(i,1) Classes.SD(i+2,1) Classes.SD(i+2,1)], [0 40000 40000 0], 'r', 'Facealpha',0.3, 'EdgeAlpha',0)
            hold on
        end
        Dhosp = [Classes.Hosp(1,:) - Prelim.Hosp(end-1,:) ; Classes.Hosp(2:end,:) - Classes.Hosp(1:end-1,:)];
        It = sum(Classes.IS1,2) + sum(Classes.IS2,2) + sum(Classes.IS3,2) + sum(Classes.IPH1,2) + sum(Classes.IPH2,2) + sum(Classes.IPH3,2) + sum(Classes.IH,2);
        %plot(Classes.t, sum(Dhosp,2), 'k', 'LineWidth', 2.5)
        %plot(Classes.t, It, 'k', 'LineWidth', 2.5)
        plot(Classes.t, sum(Classes.Hosp,2)/1000, 'k', 'LineWidth', 2.5)
        
        if vstarts < maxtime
            if strat == 1
                xline(para.vstart,'-',{'Vaccine','Arrival'},'Color','b','Linewidth',2,'FontSize',16,'Interpreter','latex','LabelOrientation','horizontal','LabelHorizontalAlignment','left')
            else
                xline(para.vstart,'-','Color','b','Linewidth',2)
            end
        end
        
        %axis([0 maxtime 0 max([max(sum(Dhosp,2)),30])])
        axis([0 maxtime 0 max([max(It),30])])
        axis([0 maxtime 0 max([1.01*sum(Classes.Hosp(end,:))/1000,30])])
        
        if strat == strategies(end)
            xlabel('Time (days)')
        end
    
        if strat == 1
            title('Cumulative $I^H(t)$ (thousands)');
            %title('Cumulative Cases $I^S(t)$');
        end
    
        set(gca,'FontSize',16)
        grid on

    end
    xticks(0:180:maxtime)
    xtickangle(0)
end
toc

%save figure
if not(isfolder('sim_images'))
    mkdir('sim_images')
end

saveas(gcf,strcat('./sim_images/','simulation_',num2str(para.vstart),'.png'))

%save thresholds used to mat file for reproducibility
save("./mats/Thresholds.mat","thresholds",'-mat')

%% Plotting cumulative vaccinations
% Uncomment to reproduce Figure S4 from manuscript

% f2 = figure(2);
% f2.Position = [1250 400 450 300];
% plot(Classes.t, Classes.V(:,1)./1000, 'LineWidth', 2.5, 'DisplayName', '0-19')
% hold on
% plot(Classes.t, Classes.V(:,2)./1000, 'LineWidth', 2.5, 'DisplayName', '20-64')
% hold on
% plot(Classes.t, Classes.V(:,3)./1000, 'LineWidth', 2.5, 'DisplayName', '65+')
% hold on
% plot(Classes.t, sum(Classes.V,2)./1000, 'k', 'LineWidth', 2.5, 'DisplayName', 'Total')
% 
% axis([0 maxtime 0 sum(para.N)./1000])
% xline(para.vstart,'--','DisplayName','Arrival','Color','k')
% xlabel('Time (days)')
% ylabel('Population (thousands)')
% title('Cumulative Vaccinations')
% legend('Interpreter','latex','Location','west')
% grid on

% saveas(f2,'./vacc_images/Tvacc.png')
