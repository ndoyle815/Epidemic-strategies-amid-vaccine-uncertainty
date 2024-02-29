% script to define model and control parameters into a mat file
clear; close all

if not(isfolder('mats'))
    mkdir('mats')
end

% model parameters
gamma = 1/5;                     % infectious period (at one point 0.4614)
epsilon = 1/5.28;                % incubation period (at one point 0.3052, 1/4)
omega = 1/800;                   % recovered period
zeta = gamma;                    % period between symptoms and hospital (at one point 1/7.1)
delta = 1/8.78;                  % hospitalisation period
tau = 0.25;                      % relative infectiousness of asymptomatic
rho = 0.1;                       % relative infectiousness of hospitalised
da = [0.0275; 0.1350; 0.5374];   % probability of symptomatic infection
ha = [0.1268; 0.1173; 0.2430];   % probability of hospitalisation given symptomatic infection
N = [233000; 580000; 187000];    % population structure
n = size(N,1);                   % number of age classes
bedsper1000 = 2.5;               % UK hospital beds per 1,000 population
Hmax = bedsper1000*sum(N)/1000;  % capacity

DIa = [0; 0; 0];                 % probability of death given (non-hospitalised) symptomatic infection
DHa = [0.0276; 0.0533; 0.1935];  % probability of death given hospitalised

% transmission matrix
beta = [6.8045, 4.0791, 0.2121; 1.9898, 8.1466, 0.55565; 0.5221, 3.4354, 2.0034].';
R0 = 3.0; % basic reproduction number

% initial exposed
E0 = round(N.*1e-4);

% control parameters
init = 0;                        % default control state (no control)
tgap = 18;                       % decision-making gap (relaxation)
tdelay = 3;                      % natural delay in implementing decision
tdiff = 7;                       % tgap-tdiff = decision-making gap (reintroduction)
ICRED = 0.4;                     % reduction in transmission when in Intermediate Control
LKRED = 0.7;                     % reduction in transmission when in Lockdown

% default (dummy) switching thresholds
T01 = 10000;                     % No Control -> Intermediate Control
T10 = 10000;                     % Intermediate Control -> No Control
T12 = 10000;                     % Intermediate Control -> Lockdown
T21 = 10000;                     % Lockdown -> Intermediate Control

% Default time to run model for (preliminary run)
t0 = 0;
maxtime = 30;

% generate fixed vaccine parameters
% parameters varying sigmoid curve (ie. rollout speed)
tc = 100;                        % time to complete 50% vaccination
kappa = 0.05;                    % logistic shaping parameter
stagger = tc/2;                  % time between vaccination commencement of age groups
efficacy = 0.9;                  % default vaccine efficacy
vstart = 2000;                   % default vaccine arrival date

save("./mats/Parameters.mat","beta","gamma","epsilon","omega","zeta","delta",...
    "tau","rho","da","ha","N","n","Hmax","E0","init","R0","DIa","DHa",...
    "tgap","tdelay","tdiff","T01","T10","T12","T21","t0","maxtime","tc","kappa",...
    "stagger","efficacy","vstart","ICRED","LKRED",'-mat')

%rescale transmission matrix to desired R0
defaultpara = load('./mats/Parameters.mat');
R0scale = Get_R0(defaultpara);

beta = (R0/R0scale).*beta;
save("./mats/Parameters.mat","beta",'-append')
