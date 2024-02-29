% Script to evaluate control strategies with vaccination, taking the
% vaccine joint distribution outcomes into account

% NB: this is setup for running Fig5.m
clear; close all

%Plotting preferences
set(0,'defaultlinelinewidth',2)
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultTextInterpreter','latex')
set(0,'defaultaxesfontsize',18)

if not(exist('mats/Thresholds.mat','file'))
    disp('No strategy thresholds saved: Running Simulations.m')
    Simulations
end

% load default parameters
para0 = load('./mats/Parameters.mat');

% vaccination start times
vstart_times = 60:60:1080;

% Define time to run model for
t_init = 30;     % preliminary run
maxtime = max(vstart_times);  % main simulation

% define strategy numbers and switching thresholds
load('mats/Thresholds.mat')
strategies = 1:length(thresholds);

% add control thresholds defined by strategy
para = para0;
para.maxtime = maxtime;
para.Hmax = 1250;        % modify hospital capacity

% define vaccine efficacy
etas = 0.0:0.05:1.0;

% define functional weights
weights = 0.0:0.01:1.0;
w2 = 2;

% stores cost function outputs
ns = length(strategies);
nw = length(weights);
nT = length(vstart_times);
ne = length(etas);
fs = zeros(nw,ns,nT,ne);

tic
for strat = strategies
    strat
    % set switching thresholds
    para.T10 = thresholds(strat,1);
    para.T01 = thresholds(strat,2);
    para.T21 = thresholds(strat,3);
    para.T12 = thresholds(strat,4);

    % run preliminary simulation to get ICs
    [Prelim, Prelim_ICs] = Get_ICs(para0);

    % starting control state
    if sum(Prelim.IH(end,:)) < para.T12
        para.init = 1;
    else
        para.init = 2;
    end

    for T = 1:nT
        para.vstart = vstart_times(T);

        for e = 1:ne
            para.efficacy = etas(e);

            % Run main simulation
            [Classes, burden, stringency, peak_hospital] = ODEmodel(para, Prelim_ICs);
    
            for w = 1:nw
                % evaluate cost function
                fs(w,strat,T,e) = CostFunction([weights(w), w2], para, burden, stringency, peak_hospital, 0);
            end
        end
    end
end
toc

%% save tensor to compute expectations elsewhere
save('./mats/cost_tensor.mat',"fs")
