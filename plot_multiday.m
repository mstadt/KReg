% Plot simulation results from driver_multiday.m

clear all;

%% load data
f1 = './MultiDaySim/29-Sep-2023_driver_multiday_insulin-1_Kamt_meal-26_ndays-50_notes-control.mat';
f2 = './MultiDaySim/29-Sep-2023_driver_multiday_insulin-1_Kamt_meal-104_ndays-50_notes-lower_highK.mat';
f3 = './MultiDaySim/03-Oct-2023_driver_multiday_insulin-1_Kamt_meal-104_highKeff-1_ndays-50_notes-PT_GFR.mat';
f4 = './MultiDaySim/03-Oct-2023_driver_multiday_insulin-1_Kamt_meal-104_highKeff-2_ndays-50_notes-PTonly.mat';
f5 = './MultiDaySim/03-Oct-2023_driver_multiday_insulin-1_Kamt_meal-104_highKeff-3_ndays-50_notes-GFRonly.mat';

dat1 = load(f1);
dat2 = load(f2);
dat3 = load(f3);
dat4 = load(f4);
dat5 = load(f5);

lab1 = 'Control K^+';
lab2 = 'High K^+ - no PT effect';
lab3 = 'High K^+ - PT + GFR effects';
lab4 = 'High K^+ - only PT effect';
lab5 = 'High K^+ - only GFR effect';

%% All the days
T_all1 = []; Y_all1 = [];
T_all2 = []; Y_all2 = [];
T_all3 = []; Y_all3 = [];
T_all4 = []; Y_all4 = [];
T_all5 = []; Y_all5 = [];
for ii = 1:dat1.n_days
    temp = dat1.Tvals{ii}./60 + (24)*(ii - 1);
    T_all1 = [T_all1; temp];
    Y_all1 = [Y_all1; dat1.Yvals{ii}];
    temp = dat2.Tvals{ii}./60 + 24*(ii - 1);
    T_all2 = [T_all2; temp];
    Y_all2 = [Y_all2; dat2.Yvals{ii}];
    temp = dat3.Tvals{ii}./60 + 24 * (ii - 1);
    T_all3 = [T_all3; temp];
    Y_all3 = [Y_all3; dat3.Yvals{ii}];
    temp = dat4.Tvals{ii}./60 + 24 * (ii - 1);
    T_all4 = [T_all4; temp];
    Y_all4 = [Y_all4; dat4.Yvals{ii}];
    temp = dat5.Tvals{ii}./60 + 24 * (ii - 1);
    T_all5 = [T_all5; temp];
    Y_all5 = [Y_all5; dat5.Yvals{ii}];
end
% convert to days
T_all1 = T_all1./24; T_all2 = T_all2./24; 
T_all3 = T_all3./24; T_all4 = T_all4./24;
T_all5 = T_all5./24;

%% Make figures
figure(1)
clf;
nr = 2; nc = 2;
lw = 3; lwgray = 3; lsgray = ':';
f.labs = 18; f.xlab = 18; f.ylab = 18; f.gca = 18; f.leg = 16; f.title = 22;
cmap = parula(6);
c1 = cmap(1,:); c2 = cmap(2,:); c3 = cmap(3,:);c4 = cmap(4,:);c5=cmap(5,:);
cgraymap = gray(5);
cgray = cgraymap(3,:);
subplot(nr,nc,1)
hold on
plot(T_all1, Y_all1(:,1), 'linewidth', lw, 'color', c1)
plot(T_all2, Y_all2(:,1), 'linewidth', lw, 'color', c2)
plot(T_all3, Y_all3(:,1), 'linewidth', lw, 'color', c3)
plot(T_all4, Y_all4(:,1), 'linewidth', lw, 'color', c4)
plot(T_all5, Y_all5(:,1), 'linewidth', lw, 'color', c5)
xlabel('Time (days)')
ylabel('Gut amount')
title('Gut amount')
grid on

subplot(nr,nc,2)
hold on
plot(T_all1,Y_all1(:,2)/dat1.pars.V_plasma, 'linewidth',lw,'color',c1)
plot(T_all2,Y_all2(:,2)/dat2.pars.V_plasma, 'linewidth',lw,'color',c2)
plot(T_all3,Y_all3(:,2)/dat3.pars.V_plasma, 'linewidth',lw,'color',c3)
plot(T_all4,Y_all4(:,2)/dat4.pars.V_plasma, 'linewidth',lw,'color',c4)
plot(T_all5,Y_all5(:,2)/dat5.pars.V_plasma, 'linewidth',lw,'color',c5)
yline(3.5,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
yline(5.0,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
xlabel('Time (days)', 'fontsize', f.xlab)
ylabel('Plasma [K^+]', 'fontsize', f.ylab)
title('Plasma [K^+]', 'fontsize', f.title)
grid on

subplot(nr,nc,3)
hold on
plot(T_all1,Y_all1(:,3)/dat1.pars.V_interstitial,'linewidth',lw,'color',c1)
plot(T_all2,Y_all2(:,3)/dat2.pars.V_interstitial,'linewidth',lw,'color',c2)
plot(T_all3,Y_all3(:,3)/dat3.pars.V_interstitial,'linewidth',lw,'color',c3)
plot(T_all4,Y_all4(:,3)/dat4.pars.V_interstitial,'linewidth',lw,'color',c4)
plot(T_all5,Y_all5(:,3)/dat5.pars.V_interstitial,'linewidth',lw,'color',c5)
yline(3.5,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
yline(5.0,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
set(gca, 'fontsize', f.gca)
xlabel('Time (days)', 'fontsize', f.xlab)
ylabel('Interstitial [K^+]', 'fontsize', f.ylab)
title('Interstitial [K^+]', 'fontsize', f.title)
grid on

subplot(nr,nc,4)
hold on
plot(T_all1,Y_all1(:,4)/dat1.pars.V_muscle,'linewidth',lw,'color',c1)
plot(T_all2,Y_all2(:,4)/dat2.pars.V_muscle,'linewidth',lw,'color',c2)
plot(T_all3,Y_all3(:,4)/dat3.pars.V_muscle,'linewidth',lw,'color',c3)
plot(T_all4,Y_all4(:,4)/dat4.pars.V_muscle,'linewidth',lw,'color',c4)
plot(T_all5,Y_all5(:,4)/dat5.pars.V_muscle,'linewidth',lw,'color',c5)
yline(120,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
yline(140,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
set(gca, 'fontsize', f.gca)
xlabel('Time (days)', 'fontsize', f.xlab)
ylabel('Intracellular [K^+]', 'fontsize', f.ylab)
title('Intracellular [K^+]', 'fontsize', f.title)
grid on

legend({lab1, lab2,lab3,lab4,lab5})
%%
figure(2)
clf;
nr = 1; nc = 2;
lw = 3; lwgray = 3; lsgray = ':';
cmap = parula(4);
cmap2 = spring(3);
c1 = cmap2(1,:); c2 = cmap2(2,:); 
c3 = cmap(1,:);c4 = cmap(2,:);c5 = cmap(3,:);
cgraymap = gray(5);
cgray = cgraymap(2,:);
subplot(nr,nc,1)
hold on
plot(T_all1,Y_all1(:,2)/dat1.pars.V_plasma, 'linewidth',lw,'color',c1)
plot(T_all2,Y_all2(:,2)/dat2.pars.V_plasma, 'linewidth',lw,'color',c2)
plot(T_all3,Y_all3(:,2)/dat3.pars.V_plasma, 'linewidth',lw,'color',c3)
plot(T_all4,Y_all4(:,2)/dat4.pars.V_plasma, 'linewidth',lw,'color',c4)
plot(T_all5,Y_all5(:,2)/dat5.pars.V_plasma, 'linewidth',lw,'color',c5)
yline(3.5,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
yline(5.0,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
set(gca, 'fontsize', f.gca)
xlabel('Time (days)', 'fontsize', f.xlab)
ylabel('Plasma [K^+]', 'fontsize', f.ylab)
title('Plasma [K^+]', 'fontsize', f.title)
grid on
%legend({lab1, lab2,lab3,lab4,lab5}, 'fontsize', f.leg, 'location', 'northwest')

subplot(nr,nc,2)
hold on
plot(T_all1,Y_all1(:,4)/dat1.pars.V_muscle,'linewidth',lw,'color',c1)
plot(T_all2,Y_all2(:,4)/dat2.pars.V_muscle,'linewidth',lw,'color',c2)
plot(T_all3,Y_all3(:,4)/dat3.pars.V_muscle,'linewidth',lw,'color',c3)
plot(T_all4,Y_all4(:,4)/dat4.pars.V_muscle,'linewidth',lw,'color',c4)
plot(T_all5,Y_all5(:,4)/dat5.pars.V_muscle,'linewidth',lw,'color',c5)
yline(120,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
yline(140,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
set(gca, 'fontsize', f.gca)
xlabel('Time (days)', 'fontsize', f.xlab)
ylabel('Intracellular [K^+]', 'fontsize', f.ylab)
title('Intracellular [K^+]', 'fontsize', f.title)
grid on

legend({lab1, lab2,lab3,lab4,lab5}, 'fontsize', f.leg, 'location', 'northwest')

AddLetters2Plots(figure(2), {'(a)', '(b)'},...
                'HShift', -0.05, 'VShift', -0.06, ...
                'fontsize', f.labs)