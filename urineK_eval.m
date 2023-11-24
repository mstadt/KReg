% Used to analyze urine output

clear all;

% get control K data
%fname = './MultiDaySim/22-Oct-2023_driver_multiday_insulin-1_Kamt_meal-26_TGFeff-3_alphaTGF-0.11694_etaPTKreab-0.67_ndays-50_notes-control.mat';

% TGF only data
fname = './MultiDaySim/22-Oct-2023_driver_multiday_insulin-1_Kamt_meal-104_TGFeff-2_alphaTGF-0.11694_etaPTKreab-0.36_ndays-50_notes-GFRonly.mat';
dat = load(fname);

% Get first 5 days of data separately

% IDEAS:
%   may need to use trapz to get the area under curve
%   what is the actual urine K rate needed (if stayed stable) to 
%     match the control K?

% Convert data to days
T_all = []; Y_all = []; Urine = []; 


%% Get Day 1 values
day_num = 40;
fprintf('Day %i \n', day_num)
Y = dat.Yvals{day_num};
v = compute_kidney_vars(Y, dat.params, ...
                'do_MKX', [dat.doMKX, dat.MKXslope], ...
                'do_FF', dat.doFF, ...
                'TGF_eff', [dat.TGF_eff, ...
                            dat.alpha_TGF, ...
                            dat.eta_ptKreab]);
T = dat.Tvals{day_num}./60; % time in hours

%% plot results
f.labs = 18; f.xlab = 18; f.ylab = 18; f.gca = 18; f.leg = 16; f.title = 22;
lw = 3; lwgray = 3; lsgray = ':';
cgraymap = gray(5);
cgray = cgraymap(3,:);

% Urine figure
figure(1)
clf;
plot(T, v.UrineK, 'linewidth', lw)
set(gca, 'fontsize', f.gca)
xlabel('Time (hours)')
ylabel('Urine K (mmol/min)')
grid on

% Urine sum
x = T * 60;
y = v.UrineK;
A = trapz(x,y);
disp(A)
fprintf('area of urine K: %.5f \n', A)

%fprintf('Sum of urine K: %.5f \n', sum(v.UrineK))

% K plasma
figure(2)
clf;
hold on
plot(T, Y(:,2)/dat.pars.V_plasma, 'linewidth', lw)
yline(3.5, 'color', cgray, 'linestyle', lsgray, 'linewidth',lwgray)
yline(5.0, 'color', cgray, 'linestyle', lsgray, 'linewidth',lwgray)
set(gca, 'fontsize', f.gca)
xlabel('Time (hours)')
ylabel('Plasma [K^+] (mmol/L)')
grid on

% K intracellular
figure(3)
clf;
hold on
plot(T, Y(:,4)/dat.pars.V_muscle, 'linewidth', lw)
yline(120, 'color', cgray, 'linestyle', lsgray, 'linewidth',lwgray)
yline(140, 'color', cgray, 'linestyle', lsgray, 'linewidth',lwgray)
set(gca, 'fontsize', f.gca)
xlabel('Time (hours)')
ylabel('Intracellular [K^+] (mmol/L)')
grid on

%% Total body K
figure(5)
clf;
total_K = Y(:,1) + Y(:,2) + Y(:,3) + Y(:,4);
plot(T, total_K, 'linewidth', lw)
yline(total_K(1), 'color', cgray, 'linestyle', lsgray, 'linewidth', lwgray)
set(gca, 'fontsize', f.gca)
xlabel('Time (hours)')
ylabel('Total body K^+ (mmol)')
grid on
fprintf('Total body K+ diff: %0.3f \n', total_K(1) - total_K(end))

%% Gut K
figure(4)
plot(T, Y(:,1), 'linewidth',lw)
set(gca, 'fontsize', f.gca)
xlabel('Time (hours)')
ylabel('Gut K^+ amount (mmol)')
grid on
fprintf('Gut amount diff: %.3f \n', Y(1,1) - Y(end,1))

%% Plasma K amount
figure(6)
clf;
plot(T, Y(:,2), 'linewidth',lw)
yline(Y(1,2),'linewidth',lwgray, 'color',cgray, 'linestyle', lsgray)
set(gca, 'fontsize', f.gca)
xlabel('Time (hours)')
ylabel('Plasma K^+ amount (mmol)')
grid on
fprintf('Plasma amount diff: %0.3f \n', Y(1,2) - Y(end,2))

%% Intracellular K amount
figure(7)
clf;
plot(T, Y(:,4), 'linewidth',lw)
yline(Y(1,4),'linewidth',lwgray, 'color',cgray, 'linestyle', lsgray)
set(gca, 'fontsize', f.gca)
xlabel('Time (hours)')
ylabel('Intracellular K^+ amount (mmol)')
grid on
fprintf('Intracellular amount diff: %0.3f \n', Y(1,4) - Y(end,4))

%% Intracellular K amount
figure(8)
clf;
plot(T, Y(:,3), 'linewidth',lw)
yline(Y(1,3),'linewidth',lwgray, 'color',cgray, 'linestyle', lsgray)
set(gca, 'fontsize', f.gca)
xlabel('Time (hours)')
ylabel('Interstitial K^+ amount (mmol)')
grid on
fprintf('Interstitial amount diff: %0.3f \n', Y(1,3) - Y(end,3))

%% GFR
figure(9)
clf;
plot(T, v.GFR, 'linewidth',lw)
set(gca, 'fontsize', f.gca)
xlabel('Time (hrs)')
ylabel('GFR')
grid on
