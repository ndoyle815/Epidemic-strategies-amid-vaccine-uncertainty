% script to reproduce Figure S3 (qualitative behaviour of exponential term
% in objective function) from the manuscript
clear

%Plotting preferences
set(0,'defaultlinelinewidth',2)
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultTextInterpreter','latex')
set(0,'defaultaxesfontsize',16)

wHs = [0.2 2 10];
Hmax = 1250;

hosp_peaks = (Hmax-15):0.01:(Hmax+5);

f = figure(1);
f.Position = [1250 400 450 300];
for wH = wHs
    plot(hosp_peaks, exp(wH*(hosp_peaks - Hmax)),'LineWidth',2.5)
    hold on
end
patch([min(hosp_peaks) min(hosp_peaks) max(hosp_peaks) max(hosp_peaks)], [0 1 1 0], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none')
xline(Hmax,'--','$H_c$','color','k','LineWidth',2,'Interpreter','latex','FontSize',16,'LabelOrientation','horizontal','LabelHorizontalAlignment','left')
axis([min(hosp_peaks) max(hosp_peaks) 0 5])
xlabel('Peak Hospitalisations')
ylabel('Cost')
legend({'$w_H = 0.2$', '$w_H = 2$', '$w_H = 10$'},'Location','west','Interpreter','latex','FontSize',16)

%save figure
if not(isfolder('sim_images'))
    mkdir('sim_images')
end

saveas(gcf,'./sim_images/expterm.png')
