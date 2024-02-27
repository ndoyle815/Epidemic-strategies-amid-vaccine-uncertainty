% script to rank the control strategies under each vaccine-related joint
% distribution according to whether it minimises the mean, mode or a given
% percentile cost (eg. 95th). Reproduces Figure S8 from the manuscript
clear

%Plotting preferences
set(0,'defaultlinelinewidth',2)
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultTextInterpreter','latex')
set(0,'defaultaxesfontsize',18)

%for saving figures
if not(isfolder('vacc_images'))
    mkdir('vacc_images')
end

% load cost outputs
if not(exist('mats/cost_tensor.mat','file'))
    disp('No cost outputs saved: Running dists_joint.m')
    dists_joint
end

load('./mats/cost_tensor.mat');

% load joint distributions
if not(exist('mats/jointdists.mat','file'))
    disp('No vaccine pmfs saved: Running FigS1S2_jointdists.m')
    FigS1S2_jointdists
end

load('./mats/jointdists.mat');

% strategy markers
markers = {'o','^','s','d'};
cols = [0.9290 0.6940 0.1250; 0.3290, 0.6940, 0.1250; 0.4940 0.1840 0.5560; 0 0.5470 0.9410];

%% Acquire costs based on each metric

strategies = [1:size(fs,2)];
weights = [0.0:0.01:1.0];
etas = [0.0:0.05:1.0];
Ts = [60:60:1080];
dweights = [1:5:size(fs,1)];  % only plot some weights with markers

ndists = size(P,1);
percentile = 0.95;
dc = 0.0001;

Meancosts     = zeros(length(ndists), length(weights), length(strategies));
Modecosts     = zeros(length(ndists), length(weights), length(strategies));
Quantilecosts = zeros(length(ndists), length(weights), length(strategies));
Meanranks     = zeros(length(ndists), length(weights), length(strategies));
Moderanks     = zeros(length(ndists), length(weights), length(strategies));
Quantileranks = zeros(length(ndists), length(weights), length(strategies));

Plowest        = zeros(length(ndists), length(weights), length(strategies));
Plowestranks   = zeros(length(ndists), length(weights), length(strategies));
Phighest       = zeros(length(ndists), length(weights), length(strategies));
Phighestranks  = zeros(length(ndists), length(weights), length(strategies));
Worstcosts     = zeros(length(ndists), length(weights), length(strategies));
Worstranks     = zeros(length(ndists), length(weights), length(strategies));

for dist = 1:ndists

    % obtain joint distribution
    PD = reshape(P(dist,:,:), [length(Ts) length(etas)]);

    % and find its mode for later (finds multiple in case of floating point
    % error)
    maxval = max(PD(:));
    idx = find(maxval - PD(:)  < 1e-10);

    for w = 1:length(weights)

        % find which strategy is lowest per (T,eta) combination
        [~,minix] = min([fs(w,1,:,:);fs(w,2,:,:);fs(w,3,:,:);fs(w,4,:,:)]);
        minix = reshape(minix, [length(Ts) length(etas)]);

        % find which strategy is highest per (T,eta) combination
        [~,maxix] = max([fs(w,1,:,:);fs(w,2,:,:);fs(w,3,:,:);fs(w,4,:,:)]);
        maxix = reshape(maxix, [length(Ts) length(etas)]);

        for strat = strategies

            % obtain function values
            FWS = reshape(fs(w,strat,:,:), [length(Ts) length(etas)]);

            % mean cost
            Meancosts(dist,w,strat) = sum(PD.*FWS,'all');

            % mode cost (average if multiple modes)
            Modecosts(dist,w,strat) = mean(FWS(idx));

            % percentile cost
            C = 1.0;
            while sum(PD(FWS < C)) > percentile
                C = C - dc;
            end
            Quantilecosts(dist,w,strat) = C;
            
            % probability cost is lowest
            optlow = minix == strat;
            Plowest(dist,w,strat) = sum(PD.*optlow,'all');

            % probability cost is highest
            opthigh = maxix == strat;
            Phighest(dist,w,strat) = sum(PD.*opthigh,'all');

            % worst cost (T = 1,080 days, eta = 0.1)
            % essentially no vaccine as time horizon is = T
            Worstcosts(dist,w,strat) = FWS(end,1);

        end

        % rank mean cost
        [~, ii] = sort(Meancosts(dist,w,:));
        [~, Meanranks(dist,w,:)] = sort(ii);

        % rank mode cost
        [~, ii] = sort(Modecosts(dist,w,:));
        [~, Moderanks(dist,w,:)] = sort(ii);

        % rank percentile cost
        [~, ii] = sort(Quantilecosts(dist,w,:));
        [~, Quantileranks(dist,w,:)] = sort(ii);

        % rank probability of lowest cost
        [~, ii] = sort(Plowest(dist,w,:),'descend');
        [~, Plowestranks(dist,w,:)] = sort(ii);

        % rank probability of highest cost
        [~, ii] = sort(Phighest(dist,w,:));
        [~, Phighestranks(dist,w,:)] = sort(ii);

        % rank worst cost
        [~, ii] = sort(Worstcosts(dist,w,:));
        [~, Worstranks(dist,w,:)] = sort(ii);

    end

    f = figure();
    f.Position = [200 800 800 600];
    if dist == 2
        sgtitle('(D2) Optimistic Vaccine Distribution','FontSize',20,'Interpreter','Latex')
    else
        sgtitle(strcat('Distribution',{' '},num2str(dist)),'FontSize',20,'Interpreter','Latex')
    end

    subplot(4,1,1)

    for strat = strategies
        plot(weights,Meanranks(dist,:,strat),'-o','Color',cols(strat,:),'MarkerSize',6,'MarkerFaceColor',cols(strat,:))
        hold on
    end
    set(gca, 'YDir','reverse')
    axis([min(weights) max(weights) 0.5 4.5])
    xticks(weights(1:10:end))
    xticklabels([])
    xtickangle(0)
    xlim([0.25 0.75])
    yticks([1,2,3,4])
    ylabel('Ranking');
    title('Minimise expected (mean) cost')
    grid on

    subplot(4,1,2)

    for strat = strategies
        plot(weights,Moderanks(dist,:,strat),'-o','Color',cols(strat,:),'MarkerSize',6,'MarkerFaceColor',cols(strat,:))
        hold on
    end
    set(gca, 'YDir','reverse')
    axis([min(weights) max(weights) 0.5 4.5])
    xticks(weights(1:10:end))
    xticklabels([])
    xtickangle(0)
    xlim([0.25 0.75])
    yticks([1,2,3,4])
    ylabel('Ranking')
    title('Minimise most likely cost')
    grid on

    subplot(4,1,3)

    for strat = strategies
        plot(weights,Plowestranks(dist,:,strat),'-o','Color',cols(strat,:),'MarkerSize',6,'MarkerFaceColor',cols(strat,:))
        hold on
    end
    set(gca, 'YDir','reverse')
    axis([min(weights) max(weights) 0.5 4.5])
    xticks(weights(1:10:end))
    xticklabels([])
    xtickangle(0)
    xlim([0.25 0.75])
    yticks([1,2,3,4])
    ylabel('Ranking')
    title(strcat('Maximise Probability of lowest cost'))
    grid on
    
%     subplot(6,1,4)
% 
%     for strat = strategies
%         plot(weights,Phighestranks(dist,:,strat),'-o','Color',cols(strat,:),'MarkerSize',6,'MarkerFaceColor',cols(strat,:))
%         hold on
%     end
%     set(gca, 'YDir','reverse')
%     axis([min(weights) max(weights) 0.5 4.5])
%     xticks(weights(1:10:end))
%     xticklabels([])
%     xtickangle(0)
%     yticks([1,2,3,4])
%     %xlabel('Weight $w_1$')
%     ylabel('Ranking')
%     title(strcat('Minimise Probability of highest cost'))
%     grid on

    subplot(4,1,4)

    for strat = strategies
        plot(weights,Quantileranks(dist,:,strat),'-o','Color',cols(strat,:),'MarkerSize',6,'MarkerFaceColor',cols(strat,:))
        hold on
    end
    set(gca, 'YDir','reverse')
    axis([min(weights) max(weights) 0.5 4.5])
    xticks(weights(1:10:end))
    xtickangle(0)
    xlim([0.25 0.75])
    yticks([1,2,3,4])
    xlabel('Weight $w$')
    ylabel('Ranking')
    title(strcat('Minimise',{' '},num2str(percentile*100),'th percentile cost'))
    grid on

%     subplot(6,1,6)
% 
%     for strat = strategies
%         plot(weights,Worstranks(dist,:,strat),'-o','Color',cols(strat,:),'MarkerSize',6,'MarkerFaceColor',cols(strat,:))
%         hold on
%     end
%     set(gca, 'YDir','reverse')
%     axis([min(weights) max(weights) 0.5 4.5])
%     xticks(weights(1:10:end))
%     xticklabels([])
%     xtickangle(0)
%     yticks([1,2,3,4])
%     xlabel('Weight $w$')
%     ylabel('Ranking')
%     title(strcat('Minimise worst cost (no vaccine)'))
%     grid on 


    saveas(f,strcat('./ranking_images/rankings_dist',num2str(dist),'.png'))
    close all

end

