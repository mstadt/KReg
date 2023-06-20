clear all;

%--------------------------
% Begin user input
%--------------------------

% initialize parameter values
fprintf('loading params \n')
pars = set_params();
[params, parnames] = pars2vector(pars, 0);

% change parameters here!

%% set initial conditions
temp = load('./SS/SS1.mat');
IC = temp.SS;

%-----------------------------
% End user input
%-----------------------------
% ODE options
tspan = [0 4000];
options = odeset('RelTol',1.0e-6,'AbsTol',1e-9);

%% solve ODE system
fprintf('solving ODEs \n')
[t,y] = ode15s(@(t,y) kreg_eqns(t,y,params,...
                            'SS', true), ...
                            tspan, IC, options);


%% make figures
fprintf('making figures \n')
% figure specs
lw = 3;
f.xlab = 16; f.ylab = 16; f.title = 18;
cmap = parula(5);
c1 = cmap(3,:);
cgraymap = gray(5);
cgray = cgraymap(3,:);
lwgray = 2; lsgray = '--';

% variables
figure(1)
clf
nrows = 2; ncols  = 3;
subplot(nrows,ncols,1)
plot(t,y(:,1),'linewidth',lw,'color',c1)
ylabel('M_{Kgut}', 'fontsize', f.ylab)
xlabel('t', 'fontsize', f.xlab)
title('Gut K', 'fontsize', f.title)
grid on

subplot(nrows,ncols,2)
plot(t,y(:,2),'linewidth',lw,'color',c1)
ylabel('M_{Kplas}', 'fontsize', f.ylab)
xlabel('t', 'fontsize', f.xlab)
title('Plasma K', 'fontsize', f.title)
grid on

subplot(nrows,ncols,3)
plot(t,y(:,3),'linewidth',lw,'color',c1)
ylabel('M_{Kinter}', 'fontsize', f.ylab)
xlabel('t', 'fontsize', f.xlab)
title('Interstitial K', 'fontsize', f.title)
grid on

subplot(nrows,ncols,4)
plot(t,y(:,4),'linewidth',lw,'color',c1)
ylabel('M_{Kmuscle}', 'fontsize', f.ylab)
xlabel('t', 'fontsize', f.xlab)
title('Muscle K', 'fontsize', f.title)
grid on

subplot(nrows,ncols,5)
plot(t,y(:,5),'linewidth',lw,'color',c1)
ylabel('N_{al}', 'fontsize', f.ylab)
xlabel('t', 'fontsize', f.xlab)
title('Normalized ALD', 'fontsize', f.title)
grid on

% concentrations
figure(2)
clf
nrows = 1; ncols  = 3;
subplot(nrows,ncols,1)
plot(t,y(:,2)/pars.V_plasma,'linewidth',lw,'color',c1)
hold on
yline(3.5,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
yline(5.0,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
ylabel('[K^+]_{plasma}', 'fontsize', f.ylab)
xlabel('t', 'fontsize', f.xlab)
title('Plasma K concentration', 'fontsize', f.title)
grid on

subplot(nrows,ncols,2)
plot(t,y(:,3)/pars.V_interstitial,'linewidth',lw,'color',c1)
hold on
yline(3.5,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
yline(5.0,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
ylabel('[K^+]_{inter}', 'fontsize', f.ylab)
xlabel('t', 'fontsize', f.xlab)
title('Interstitial K concentration', 'fontsize', f.title)
grid on

subplot(nrows,ncols,3)
plot(t,y(:,4)/pars.V_muscle,'linewidth',lw,'color',c1)
hold on
yline(120,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
yline(140,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
ylabel('[K^+]_{muscle}', 'fontsize', f.ylab)
xlabel('t', 'fontsize', f.xlab)
title('Muscle K concentration', 'fontsize', f.title)
grid on

