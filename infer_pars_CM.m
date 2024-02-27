% script to calibrate parameters da,ha,mua,delta in addition to contact matrix
% from available data
clear

UKpop = readtable('./data/United Kingdom-2020.csv');

ageY = 1:4;
ageA = 5:13;
ageS = 13:21;
menpop = table2array(UKpop(1:end-1,'M'));
womenpop = table2array(UKpop(1:end-1,'F'));
Tpop = menpop + womenpop;

popY = Tpop(ageY);
popA = Tpop(ageA);
popS = Tpop(ageS);

pop_propY = popY./sum(popY);
pop_propA = popA./sum(popA);
pop_propS = popS./sum(popS);

% load CHESS data from Keeling et al, 2021
load('./data/CHESSdata.mat')

da = zeros(3,1);
ha = zeros(3,1);
DHa = zeros(3,1);

da(1) = Detection(ageY)*pop_propY;
ha(1) = sum(Sympt_2_hosp(ageY).*pop_propY);
DHa(1) = sum(Hosp_2_Death(ageY).*pop_propY);
da(2) = Detection(ageA)*pop_propA;
ha(2) = sum(Sympt_2_hosp(ageA).*pop_propA);
DHa(2) = sum(Hosp_2_Death(ageA).*pop_propA);
da(3) = Detection(ageS)*pop_propS;
ha(3) = sum(Sympt_2_hosp(ageS).*pop_propS);
DHa(3) = sum(Hosp_2_Death(ageS).*pop_propS);

da
ha
DHa

%% Calibrate delta from hospital time distribution

days2 = [0:59]';
f2 = fit(days2,Distribution_Hosp_Time(days2+1)','exp1','Lower',[1 -1],'Upper',[1 0])

figure(2)
bar(days2,Distribution_Hosp_Time(days2+1),"blue","FaceAlpha",0.3,"EdgeColor","none")
hold on
plot(f2,days2,Distribution_Hosp_Time(days2+1))
axis([-1 60 0 1.1])
title('Time in Hospital')

delta = 1/-f2.b

%% Calibrating contact matrix

CDATA = load('./data/prem_UKall.csv');
lastage = ageS(end);

CM = [sum(CDATA(:,ageY),2), sum(CDATA(:,ageA),2), sum(CDATA(:,ageS),2)]; %sum columns by age groupings
CM = [Tpop(ageY)'*CM(ageY,:)/sum(Tpop(ageY)); Tpop(ageA)'*CM(ageA,:)/sum(Tpop(ageA)); Tpop(ageS)'*CM(ageS,:)/sum(Tpop(ageS))] %weight rows by age groupings

