% Plots 2 outputs from driverMeal.m with compute_vars output

clear all;
%% load data
f1 = './MealSim/05-Jun-2023_driverMeal_insulin-0_Kin-35_notes-KClOnly_computevars.mat';
f2 = './MealSim/05-Jun-2023_driverMeal_insulin-1_Kin-35_notes-MealKCl_computevars.mat';

dat1 = load(f1);
dat2 = load(f2);

lab1 = 'KCl Only';
lab2 = 'Meal + KCl';

%% make figures
fprintf('making figures \n')
% figure specs
lw = 3;
f.xlab = 16; f.ylab = 16; f.title = 18;
f.leg = 16;
cmap = parula(5);
c1 = cmap(1,:);
c2 = cmap(3,:);
ls1 = '-'; ls2 = '-.';
cgraymap = gray(5);
cgray = cgraymap(3,:);
lwgray = 2; lsgray = '--';
leglabs = {lab1, lab2};

%% variables
figure(10)
clf
nrows = 2; ncols  = 3;
subplot(nrows,ncols,1)
hold on
plot(dat1.t,dat1.y(:,1),'linewidth',lw,'color',c1, 'linestyle',ls1)
plot(dat2.t,dat2.y(:,1),'linewidth',lw,'color',c2, 'linestyle',ls2)
ylabel('M_{Kgut}', 'fontsize', f.ylab)
xlabel('t', 'fontsize', f.xlab)
title('Gut K', 'fontsize', f.title)
grid on

subplot(nrows,ncols,2)
hold on
plot(dat1.t,dat1.y(:,2),'linewidth',lw,'color',c1, 'linestyle',ls1)
plot(dat2.t,dat2.y(:,2),'linewidth',lw,'color',c2, 'linestyle',ls2)
ylabel('M_{Kplas}', 'fontsize', f.ylab)
xlabel('t', 'fontsize', f.xlab)
title('Plasma K', 'fontsize', f.title)
grid on

subplot(nrows,ncols,3)
hold on
plot(dat1.t,dat1.y(:,3),'linewidth',lw,'color',c1, 'linestyle',ls1)
plot(dat2.t,dat2.y(:,3),'linewidth',lw,'color',c2, 'linestyle',ls2)
ylabel('M_{Kinter}', 'fontsize', f.ylab)
xlabel('t', 'fontsize', f.xlab)
title('Interstitial K', 'fontsize', f.title)
grid on

subplot(nrows,ncols,4)
hold on
plot(dat1.t,dat1.y(:,4),'linewidth',lw,'color',c1, 'linestyle',ls1)
plot(dat2.t,dat2.y(:,4),'linewidth',lw,'color',c2, 'linestyle',ls2)
ylabel('M_{Kmuscle}', 'fontsize', f.ylab)
xlabel('t', 'fontsize', f.xlab)
title('Muscle K', 'fontsize', f.title)
grid on

subplot(nrows,ncols,5)
hold on
plot(dat1.t,dat1.y(:,5),'linewidth',lw,'color',c1, 'linestyle',ls1)
plot(dat2.t,dat2.y(:,5),'linewidth',lw,'color',c2, 'linestyle',ls2)
ylabel('N_{al}', 'fontsize', f.ylab)
xlabel('t', 'fontsize', f.xlab)
title('Normalized ALD', 'fontsize', f.title)
grid on
legend(leglabs, 'fontsize', f.leg)

%% concentrations
figure(11)
clf
nrows = 1; ncols  = 3;
subplot(nrows,ncols,1)
hold on
plot(dat1.t,dat1.y(:,2)/dat1.pars.V_plasma,'linewidth',lw,'color',c1, 'linestyle',ls1)
plot(dat2.t,dat2.y(:,2)/dat2.pars.V_plasma,'linewidth',lw,'color',c2, 'linestyle',ls2)
yline(3.5,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
yline(5.0,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
ylabel('[K^+]_{plasma}', 'fontsize', f.ylab)
xlabel('t', 'fontsize', f.xlab)
title('Plasma K concentration', 'fontsize', f.title)
grid on

subplot(nrows,ncols,2)
hold on
plot(dat1.t,dat1.y(:,3)/dat1.pars.V_interstitial,'linewidth',lw,'color',c1, 'linestyle',ls1)
plot(dat2.t,dat2.y(:,3)/dat2.pars.V_interstitial,'linewidth',lw,'color',c2, 'linestyle',ls2)
yline(3.5,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
yline(5.0,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
ylabel('[K^+]_{inter}', 'fontsize', f.ylab)
xlabel('t', 'fontsize', f.xlab)
title('Interstitial K concentration', 'fontsize', f.title)
grid on

subplot(nrows,ncols,3)
hold on
plot(dat1.t,dat1.y(:,4)/dat1.pars.V_muscle,'linewidth',lw,'color',c1, 'linestyle',ls1)
plot(dat2.t,dat2.y(:,4)/dat2.pars.V_muscle,'linewidth',lw,'color',c2, 'linestyle',ls2)
yline(120,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
yline(140,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
ylabel('[K^+]_{muscle}', 'fontsize', f.ylab)
xlabel('t', 'fontsize', f.xlab)
title('Muscle K concentration', 'fontsize', f.title)
grid on
legend(leglabs, 'fontsize', f.leg)

%% feedforward and feedback response
figure(15)
clf
nrows = 2; ncols = 2;
subplot(nrows,ncols,1)
hold on
rhoins1 = [dat1.vals1.rho_insulin; dat1.vals2.rho_insulin; dat1.vals3.rho_insulin];
rhoins2 = [dat2.vals1.rho_insulin; dat2.vals2.rho_insulin; dat2.vals3.rho_insulin];
plot(dat1.t, rhoins1, 'linewidth',lw,'color',c1,'linestyle',ls1)
plot(dat2.t, rhoins2, 'linewidth',lw,'color',c2,'linestyle',ls2)
ylabel('\rho_{ins}', 'fontsize', f.ylab)
xlabel('t', 'fontsize', f.xlab)
title('[insulin] effect on \Phi_{ECtoIC} (\rho_{insulin})', 'fontsize', f.title)
grid on

subplot(nrows,ncols,2)
hold on
rhoal1 = [dat1.vals1.rho_al; dat1.vals2.rho_al; dat1.vals3.rho_al];
rhoal2 = [dat2.vals1.rho_al; dat2.vals2.rho_al; dat2.vals3.rho_al];
plot(dat1.t, rhoal1, 'linewidth',lw,'color',c1,'linestyle',ls1)
plot(dat2.t, rhoal2,'linewidth',lw,'color',c2,'linestyle',ls2)
ylabel('\rho_{al}', 'fontsize', f.ylab)
xlabel('t','fontsize',f.xlab)
title('[ALD] effect on \Phi_{ECtoIC} (\rho_{al})', 'fontsize', f.title)
grid on

subplot(nrows,ncols,3)
hold on
gamkin1 = [dat1.vals1.gamma_Kin; dat1.vals2.gamma_Kin; dat1.vals3.gamma_Kin];
gamkin2 = [dat2.vals1.gamma_Kin; dat2.vals2.gamma_Kin; dat2.vals3.gamma_Kin];
plot(dat1.t, gamkin1, 'linewidth',lw,'color',c1,'linestyle',ls1)
plot(dat2.t, gamkin2,'linewidth',lw,'color',c2,'linestyle',ls2)
ylabel('\gamma_{Kin}', 'fontsize', f.ylab)
xlabel('t', 'fontsize', f.xlab)
title('GI feedforward effect (\gamma_{Kin})', 'fontsize', f.title)
grid on

subplot(nrows,ncols,4)
hold on
gamal1 = [dat1.vals1.gamma_al; dat1.vals2.gamma_al; dat1.vals3.gamma_al];
gamal2 = [dat2.vals1.gamma_al; dat2.vals2.gamma_al; dat2.vals3.gamma_al];
plot(dat1.t, gamal1, 'linewidth', lw,'color',c1,'linestyle',ls1)
plot(dat2.t, gamal2, 'linewidth',lw,'color',c2,'linestyle',ls2)
ylabel('gamma_{al}', 'fontsize', f.ylab)
xlabel('t','fontsize',f.xlab)
title('[ALD] effect on \Phi_{dt-Ksec} (\gamma_{al})', 'fontsize', f.title)
ylabel('\gamma_{al}', 'fontsize', f.ylab)
xlabel('t', 'fontsize', f.xlab)
grid on

legend(leglabs, 'fontsize', f.leg)