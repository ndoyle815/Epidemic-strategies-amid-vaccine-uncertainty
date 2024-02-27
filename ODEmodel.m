% ODE SEIR model code, taking in arguments
% para: model parameters
% ICs: initial conditions
% control and vaccination incorporated at discrete timesteps

function [Classes, burden, stringency, peak_hospital] = ODEmodel(para,ICs)

simlen = para.maxtime - para.t0 + 1;
tnidx = 2;

%unpack ICs
pop = zeros(simlen,20*para.n);
pop(1,:) = reshape([ICs.S ICs.E1 ICs.E2 ICs.E3 ICs.IA1 ICs.IA2 ICs.IA3 ICs.IS1 ICs.IS2 ICs.IS3 ...
                    ICs.IPH1 ICs.IPH2 ICs.IPH3 ICs.IH ICs.R1 ICs.R2 ICs.R3 ICs.Cases ICs.Hosp ICs.V], ...
                    1,20*para.n);

%indexes for classes (to take out of pop, check above line if changing classes)
Sidx = 1:para.n;
IHidx = 13*para.n+1:14*para.n;
R1idx = 14*para.n+1:15*para.n;
R2idx = 15*para.n+1:16*para.n;
R3idx = 16*para.n+1:17*para.n;
HOSPidx = 18*para.n+1:19*para.n;
Vidx = 19*para.n+1:20*para.n;

%setup
tn = para.t0;
SD = [tn, para.init];  % This class will record date and control index at every switch
opts = odeset('RelTol',1e-5);

% number to be vaccinated
Vmax = para.efficacy.*para.N';

% vaccine stagger parameter
c = para.tc + para.vstart + [2, 1, 0].*para.stagger;

% captures burden and stringency costs
burden = zeros(simlen,1);
burden(tnidx-1) = sum(ICs.Hosp(end,:));
stringency = zeros(simlen,1);
stringency(tnidx-1) = 0;

% the main iteration
while tn < para.maxtime

    % acquire end-of-day outputs
    inits = pop(tnidx-1,:);

    % begin vaccine-associated decrease in transmission
    if tn >= para.vstart

        N_eligible = inits(Sidx) + inits(R1idx) + inits(R2idx) + inits(R3idx);

        % do vaccine movements (logistic function)
        vacc = (1./(1+exp(-para.kappa.*(tn-c)))).*Vmax;

        dvacc = vacc - inits(Vidx);
        Vacc_S = dvacc.*(inits(Sidx)./N_eligible);
        Vacc_R1 = dvacc.*(inits(R1idx)./N_eligible);
        Vacc_R2 = dvacc.*(inits(R2idx)./N_eligible);
        Vacc_R3 = dvacc.*(inits(R3idx)./N_eligible);
    
        inits(Sidx) = inits(Sidx) - Vacc_S;
        inits(R1idx) = inits(R1idx) - Vacc_R1;
        inits(R2idx) = inits(R2idx) - Vacc_R2;
        inits(R3idx) = inits(R3idx) - Vacc_R3;
        inits(Vidx) = inits(Vidx) + Vacc_S + Vacc_R1 + Vacc_R2 + Vacc_R3;
    end
    

    % "social distancing" control: 70% decrease in contact rates if hospital 
    % occupancy is above a given threshold (lockdown) or a 40% decrease for softer
    % restrictions and smaller thresholds (Intermediate Control)

    % Notation for states SD: 
    % 0: No Control, 1: Intermediate Control, 2: Lockdown
    % 0.5: No Control           <-> Intermediate Control
    % 1.5: Intermediate Control <-> Lockdown
    % 2.5: No Control            -> Lockdown

    % para.RIT controls the reduction in transmission

    currently_in_hosp = sum(inits(IHidx));
    currently_in_hosp_prev = sum(pop(max(tnidx-2,1),IHidx));

    if SD(end,2) == 0
        para.RIT = 0;
        if currently_in_hosp > para.T12 && tn - SD(end,1) >= para.tgap - para.tdiff
            SD(end+1,:) = [tn, 2.5];
        elseif currently_in_hosp > para.T01 && tn - SD(end,1) >= para.tgap - para.tdiff
            SD(end+1,:) = [tn, 0.5];
        end
    elseif SD(end,2) == 0.5
        para.RIT = para.ICRED*SD(end-1,2); % = 0 if previous state 0 or 0.4 if previous state 1
        if tn - SD(end,1) == para.tdelay
            SD(end+1,:) = [tn, 1-SD(end-1,2)]; % = 1 if previous state 0 or 0 if previous state 1
        end
    elseif SD(end,2) == 1
        para.RIT = para.ICRED;
        if currently_in_hosp > para.T12 && tn - SD(end,1) >= para.tgap - para.tdiff
            SD(end+1,:) = [tn, 1.5];
        elseif currently_in_hosp < para.T10 && tn - SD(end,1) >= para.tgap && currently_in_hosp < currently_in_hosp_prev
            SD(end+1,:) = [tn, 0.5];
        end
    elseif SD(end,2) == 1.5
        para.RIT = para.ICRED + (para.LKRED - para.ICRED)*(SD(end-1,2)-1); % = 0.4 if previous state 1 or 0.7 if previous state 2
        if tn - SD(end,1) == para.tdelay
            SD(end+1,:) = [tn, 3-SD(end-1,2)]; % = 2 if previous state 1 or 1 if previous state 2
        end
    elseif SD(end,2) == 2
        para.RIT = para.LKRED;
        if currently_in_hosp < para.T21 && tn - SD(end,1) >= para.tgap && currently_in_hosp < currently_in_hosp_prev
            SD(end+1,:) = [tn, 1.5];
        end
    elseif SD(end,2) == 2.5
        para.RIT = 0;
        if tn - SD(end,1) == para.tdelay
            SD(end+1,:) = [tn, 2];
        end
    end


    % Run ODE using ODE45
    [~, newpop] = ode45(@diff_SEIR_model, tn:tn+1, inits, opts, para);
    
    % add daily outputs, update burden & stringency and time
    pop(tnidx,:) = newpop(end,:);
    burden(tnidx) = sum(pop(tnidx,HOSPidx)) - sum(pop(tnidx-1,HOSPidx));
    stringency(tnidx) = para.RIT^2;

    tn = tn + 1;
    tnidx = tnidx + 1;

end

%Convert output to struct
Classes = struct('S',pop(:,1:para.n),'E1',pop(:,para.n+1:2*para.n),'E2',pop(:,2*para.n+1:3*para.n), ...
            'E3',pop(:,3*para.n+1:4*para.n),'IA1',pop(:,4*para.n+1:5*para.n),'IA2',pop(:,5*para.n+1:6*para.n),'IA3',pop(:,6*para.n+1:7*para.n), ...
            'IS1',pop(:,7*para.n+1:8*para.n),'IS2',pop(:,8*para.n+1:9*para.n),'IS3',pop(:,9*para.n+1:10*para.n), ...
            'IPH1',pop(:,10*para.n+1:11*para.n),'IPH2',pop(:,11*para.n+1:12*para.n),'IPH3',pop(:,12*para.n+1:13*para.n), ...
            'IH',pop(:,13*para.n+1:14*para.n),'R1',pop(:,14*para.n+1:15*para.n),'R2',pop(:,15*para.n+1:16*para.n),'R3',pop(:,16*para.n+1:17*para.n), ...
            'Cases',pop(:,17*para.n+1:18*para.n),'Hosp',pop(:,18*para.n+1:19*para.n),'V',pop(:,19*para.n+1:20*para.n),'SD',SD,'t',para.t0:tn);


% final calculation is the peak active hospitalisations
peak_hospital = round(max(sum(Classes.IH,2)));

%Diff equations
function dPop = diff_SEIR_model(~,pop,para)

%dPop = zeros(size(pop));

S     = pop(1 : para.n);
E1    = pop(para.n+1 : 2*para.n);
E2    = pop(2*para.n+1 : 3*para.n);
E3    = pop(3*para.n+1 : 4*para.n);
IA1   = pop(4*para.n+1 : 5*para.n);
IA2   = pop(5*para.n+1 : 6*para.n);
IA3   = pop(6*para.n+1 : 7*para.n);
IS1   = pop(7*para.n+1 : 8*para.n);
IS2   = pop(8*para.n+1 : 9*para.n);
IS3   = pop(9*para.n+1 : 10*para.n);
IPH1  = pop(10*para.n+1 : 11*para.n);
IPH2  = pop(11*para.n+1 : 12*para.n);
IPH3  = pop(12*para.n+1 : 13*para.n);
IH    = pop(13*para.n+1 : 14*para.n);
R1    = pop(14*para.n+1 : 15*para.n);
R2    = pop(15*para.n+1 : 16*para.n);
R3    = pop(16*para.n+1 : 17*para.n);
%Cases = pop(17*para.n+1 : 18*para.n);
%Hosp  = pop(18*para.n+1 : 19*para.n);
%V     = pop(19*para.n+1 : 20*para.n);

% Force of Infection
FOI = para.beta*(para.tau.*(IA1+IA2+IA3) + (IS1+IS2+IS3) + (IPH1+IPH2+IPH3) + para.rho.*IH);

% ODE equations
dS     = -(1 - para.RIT).*S.*FOI./para.N + 3*para.omega.*R3;
dE1    = (1 - para.RIT).*S.*FOI./para.N - 3*para.epsilon.*E1;
dE2    = 3*para.epsilon.*E1 - 3*para.epsilon.*E2;
dE3    = 3*para.epsilon.*E2 - 3*para.epsilon.*E3;
dIA1   = 3*para.epsilon.*(1-para.da).*E3 - 3*para.gamma.*IA1;
dIA2   = 3*para.gamma.*IA1 - 3*para.gamma.*IA2;
dIA3   = 3*para.gamma.*IA2 - 3*para.gamma.*IA3;
dIS1   = 3*para.epsilon.*para.da.*(1-para.ha).*E3 - 3*para.gamma.*IS1;
dIS2   = 3*para.gamma.*IS1 - 3*para.gamma.*IS2;
dIS3   = 3*para.gamma.*IS2 - 3*para.gamma.*IS3;
dIPH1  = 3*para.epsilon.*para.da.*para.ha.*E3 - 3*para.zeta.*IPH1;
dIPH2  = 3*para.zeta.*IPH1 - 3*para.zeta.*IPH2;
dIPH3  = 3*para.zeta.*IPH2 - 3*para.zeta.*IPH3;
dIH    = 3*para.zeta.*IPH3 - para.delta.*IH;
dR1    = 3*para.gamma.*IA3 + 3*para.gamma.*(1-para.DIa).*IS3 + (1-para.DHa).*para.delta.*IH - 3*para.omega.*R1;
dR2    = 3*para.omega.*R1 - 3*para.omega.*R2;
dR3    = 3*para.omega.*R2 - 3*para.omega.*R3;
dCases = 3*para.epsilon.*para.da.*E3;
dHosp  = 3*para.zeta.*IPH3;
dV     = [0; 0; 0];  % S and R individuals already moved in discretised manner

dPop = [dS; dE1; dE2; dE3; dIA1; dIA2; dIA3; dIS1; dIS2; dIS3; dIPH1; dIPH2; dIPH3; dIH; dR1; dR2; dR3; dCases; dHosp; dV];

