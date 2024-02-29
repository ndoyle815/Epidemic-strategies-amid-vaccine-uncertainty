% script to plot vaccine arrival date distributions and generate Figures 
% S1-2 from the manuscript
clear; close all

%Plotting preferences
set(0,'defaultlinelinewidth',2)
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultTextInterpreter','latex')
set(0,'defaultaxesfontsize',18)

Ts = [60:60:1080]';
etas = [0.0:0.05:1.0];
Tdiff = Ts(2) - Ts(1);
etadiff = etas(2) - etas(1);

% generate marginal distributions using discretenormal.m
mue1 = 5;
mut1 = 6;

vdist_eta1 = poisspdf([0:length(etas)-1], mue1);
vdist_eta2 = flip(poisspdf([0:length(etas)-1], mue1));
vdist_T1   = poisspdf([0:length(Ts)-1], mut1)';
vdist_T2   = flip(poisspdf([0:length(Ts)-1], mut1))';

% make sure to normalise
vdist_eta1 = vdist_eta1./sum(vdist_eta1);
vdist_eta2 = vdist_eta2./sum(vdist_eta2);
vdist_T1   = vdist_T1./sum(vdist_T1);
vdist_T2   = vdist_T2./sum(vdist_T2);

% eta dists
etadists = [vdist_eta1; vdist_eta2];
% T dists
Tdists = [vdist_T1 vdist_T2];

f = figure(1);
f.Position = [400 600 1200 300];
subplot(1,2,1)
bar(Ts,Tdists,1,"grouped",'FaceAlpha',0.8)
xticks(Ts(1:2:end))
xtickangle(0)
axis([min(Ts)-Tdiff max(Ts)+Tdiff 0 0.2])
xlabel('$T$ (days)')
ylabel('Probability')
title('$T$ marginal distribution')
grid on

subplot(1,2,2)
bar(etas,etadists,1,"grouped",'FaceAlpha',0.8)
xticks(etas(1:2:end))
xtickangle(0)
axis([min(etas)-etadiff max(etas)+etadiff 0 0.2])
xlabel('$\eta$')
ylabel('Probability')
title('$\eta$ marginal distribution')
grid on

%save figure
if not(isfolder('./figs/dist_images'))
    mkdir('./figs/dist_images')
end

saveas(gcf,'./figs/dist_images/vaccdists_marg.png')


%% joint distributions

P = zeros(size(Tdists,2)*size(etadists,1),size(Tdists,1),size(etadists,2));

f = figure(3);
f.Position = [400 400 size(Tdists,2)*500 size(etadists,1)*300];
sgtitle('Vaccine $(T,\eta)$ joint distributions','Fontsize',22)

k = 1;
for i = 1:size(Tdists,2)
    
    PT = Tdists(:,i);
    for j = 1:size(etadists,1)
        
        Peta = etadists(j,:);
        P(k,:,:) = PT*Peta;

        subplot(size(Tdists,2),size(etadists,1),k,'Parent',f)
        subax = gca;
        subax.Position(1) = subax.Position(1) - 0.04;

        imagesc(etas, Ts, PT*Peta)
        %set(gca,'YDir','normal')
        yticks(Ts(1:3:end))
        xticks(etas(1:4:end))
        xtickangle(0)
        clim([0,0.03])
        
        if i == size(Tdists,2)
            xlabel('Eventual coverage $\eta$')
        end
        if j == 1
            ylabel('Arrival time $T$ (days)')
        end

        k = k + 1;
    end
end

h = axes(f,'visible','off'); 
h.Title.Visible = 'on';
h.XLabel.Visible = 'on';
h.YLabel.Visible = 'on';

c = colorbar(h,'Position',[0.88 0.168 0.02 0.7],'FontSize',18,'TickLabelInterpreter','Latex');  % attach colorbar to h
colormap(c);
clim(h,[0,0.03]);
c.Ticks = [0 0.01 0.02 0.03];
c.Label.String = "Probability";
c.Label.Rotation = 270;
c.Label.VerticalAlignment = "bottom";
c.Label.Interpreter = "latex";



save('./mats/jointdists.mat',"P")

saveas(f,'./figs/dist_images/vaccdists_joint.png')
