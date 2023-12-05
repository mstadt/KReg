% slope tries (0.005, 0.01, 0.025, 0.05, 0.075, 0.1)
% notes = 'MKX', 'MKXPTeff'
clear all;
slope_vals = [0.005, 0.01, 0.025, 0.05, 0.075, 0.1];

%% endpoint for MKX, MKXPTeff K plas, K int
ends_MKX = zeros(2, length(slope_vals) + 1); % row 1: intracellular, row 2: Kplas
ends_MKXPTeff = zeros(2, length(slope_vals) + 1); % row 1: Kint, row 2: Kplas
%% Figure 1a
% MKX without PT effects
% intracellular only
% K plasma could be bar plot as end K plasma value
figure(1)
clf;
nc = 2; nr = 1;
subplot(nr,nc,1)
f.labs = 18; f.xlab = 18; f.ylab = 18; f.gca = 18; f.leg = 16; f.title = 22;
lw = 3; lwgray = 4.5; lsgray = ':';
ls = '-';
cmap = parula(length(slope_vals));
cmap2 = spring(4);
cgraymap = gray(5);
cgray = cgraymap(1,:);
hold on
for ii = 1:length(slope_vals)
    MKXslope = slope_vals(ii);
    data_date = '03-Oct-2023';
    fname = strcat('./MultiDaySim/', data_date, '_driver_multiday',...
                        '_insulin-', num2str(1),...
                        '_Kamt_meal-', num2str(104),...
                        '_MKX-', num2str(1),...
                        '_MKXSlope-', num2str(MKXslope),...
                        '_highKeff-', num2str(0),...
                        '_ndays-', num2str(50),...
                        '_notes-', 'MKX',...
                        '.mat');
    dat = load(fname);
    lab = sprintf('m_{Kic} = %0.3f', MKXslope);

    [T_all, Y_all] = convert_dat(dat);
    ends_MKX(1, ii+1) = Y_all(end, 4)./dat.pars.V_muscle; % store end intracellular value
    % plot result
    plot(T_all, Y_all(:,4)/dat.pars.V_muscle, ...
        'linewidth', lw, 'linestyle', ls, 'color', cmap(ii,:), ...
        'DisplayName',lab);
end
% plot MKX result without MKX
fname = './MultiDaySim/03-Oct-2023_driver_multiday_insulin-1_Kamt_meal-104_highKeff-0_ndays-50_notes-noMKX.mat';
dat = load(fname);
lab = 'no MKX';
[T_all, Y_all] = convert_dat(dat);
ends_MKX(1, 1) = Y_all(end,4)/dat.pars.V_muscle;
plot(T_all, Y_all(:,4)/dat.pars.V_muscle, ...
    'linewidth', lw, 'linestyle', ls, 'color', cmap2(1,:), ...
    'DisplayName', 'no MKX');
legend('fontsize', f.leg, 'location', 'northwest')
yline(120,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray, 'HandleVisibility', 'off')
yline(140,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray, 'HandleVisibility', 'off')
set(gca, 'fontsize', f.gca)
xlabel('Time (days)', 'fontsize', f.xlab)
ylabel('Intracellular [K^+] (mmol/L)', 'fontsize', f.ylab)
title('MKX only (no PT + TGF effects)', 'fontsize', f.title)
grid on
ylim([110,280])

%% Figure 1b
% MKX with PT effects
% intracellular only
% K plasma could be bar plot as end K plasma value
subplot(nr,nc,2)
hold on
for ii = 1:length(slope_vals)
    MKXslope = slope_vals(ii);
    fname = strcat('./MultiDaySim/', data_date, '_driver_multiday',...
                        '_insulin-', num2str(1),...
                        '_Kamt_meal-', num2str(104),...
                        '_MKX-', num2str(1),...
                        '_MKXSlope-', num2str(MKXslope),...
                        '_highKeff-', num2str(1),...
                        '_ndays-', num2str(50),...
                        '_notes-', 'MKXPTeff',...
                        '.mat');
    dat = load(fname);
    lab = sprintf('m_{Kic} = %0.3f', MKXslope);

    [T_all, Y_all] = convert_dat(dat);
    ends_MKXPTeff(1, ii+1) = Y_all(end, 4)/dat.pars.V_muscle;

    % plot result
    plot(T_all, Y_all(:,4)/dat.pars.V_muscle, ...
        'linewidth', lw, 'linestyle', ls, 'color', cmap(ii,:), ...
        'DisplayName',lab);
end
% plot PT GFR result without MKX
fname = './MultiDaySim/03-Oct-2023_driver_multiday_insulin-1_Kamt_meal-104_highKeff-1_ndays-50_notes-PTGFRonly.mat';
dat = load(fname);
lab = 'no MKX';
[T_all, Y_all] = convert_dat(dat);
ends_MKXPTeff(1, 1) = Y_all(end,4)/dat.pars.V_muscle;
plot(T_all, Y_all(:,4)/dat.pars.V_muscle, ...
    'linewidth', lw, 'linestyle', ls, 'color', cmap2(2,:), ...
    'DisplayName', 'PT + TGF effect (no MKX)');
legend('fontsize', f.leg, 'location', 'northwest')
yline(120,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray, 'HandleVisibility', 'off')
yline(140,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray, 'HandleVisibility', 'off')
set(gca, 'fontsize', f.gca)
xlabel('Time (days)', 'fontsize', f.xlab)
ylabel('Intracellular [K^+] (mmol/L)', 'fontsize', f.ylab)
title('MKX with PT + TGF effects', 'fontsize', f.title)
grid on
ylim([110,280])

AddLetters2Plots(figure(1), {'(A)', '(B)'},...
                'HShift', -0.05, 'VShift', -0.06, ...
                'fontsize', f.labs)

%% Figure 2 (K plasma)
figure(2)
clf;
nc = 2; nr = 1;
subplot(nr,nc,1)
f.labs = 18; f.xlab = 18; f.ylab = 18; f.gca = 18; f.leg = 16; f.title = 22;
lw = 3; lwgray = 4.5; lsgray = ':';
ls = '-';
cmap = parula(length(slope_vals));
cmap2 = spring(4);
cgraymap = gray(5);
cgray = cgraymap(1,:);
hold on
for ii = 1:length(slope_vals)
    MKXslope = slope_vals(ii);
    fname = strcat('./MultiDaySim/', data_date, '_driver_multiday',...
                        '_insulin-', num2str(1),...
                        '_Kamt_meal-', num2str(104),...
                        '_MKX-', num2str(1),...
                        '_MKXSlope-', num2str(MKXslope),...
                        '_highKeff-', num2str(0),...
                        '_ndays-', num2str(50),...
                        '_notes-', 'MKX',...
                        '.mat');
    dat = load(fname);
    lab = sprintf('m_{Kic} = %0.3f', MKXslope);

    [T_all, Y_all] = convert_dat(dat);
    ends_MKX(2,ii+1) = Y_all(end,2)/dat.pars.V_plasma;

    % plot result
    plot(T_all, Y_all(:,2)/dat.pars.V_plasma, ...
        'linewidth', lw, 'linestyle', ls, 'color', cmap(ii,:), ...
        'DisplayName',lab);
end
% plot MKX result without MKX
fname = './MultiDaySim/03-Oct-2023_driver_multiday_insulin-1_Kamt_meal-104_highKeff-0_ndays-50_notes-noMKX.mat';
dat = load(fname);
lab = 'no MKX';
[T_all, Y_all] = convert_dat(dat);
ends_MKX(2, 1) = Y_all(end,2)/dat.pars.V_plasma;
plot(T_all, Y_all(:,2)/dat.pars.V_plasma, ...
    'linewidth', lw, 'linestyle', ls, 'color', cmap2(1,:), ...
    'DisplayName', 'no MKX');
legend('fontsize', f.leg, 'location', 'northwest')
yline(3.5,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray, 'HandleVisibility', 'off')
yline(5.0,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray, 'HandleVisibility', 'off')
set(gca, 'fontsize', f.gca)
xlabel('Time (days)', 'fontsize', f.xlab)
ylabel('Plasma [K^+]', 'fontsize', f.ylab)
title('MKX with no PT + TGF effect', 'fontsize', f.title)
grid on
ylim([3.5,8.0])

% Figure 2b
subplot(nr,nc,2)
hold on
for ii = 1:length(slope_vals)
    MKXslope = slope_vals(ii);
    fname = strcat('./MultiDaySim/', data_date, '_driver_multiday',...
                        '_insulin-', num2str(1),...
                        '_Kamt_meal-', num2str(104),...
                        '_MKX-', num2str(1),...
                        '_MKXSlope-', num2str(MKXslope),...
                        '_highKeff-', num2str(1),...
                        '_ndays-', num2str(50),...
                        '_notes-', 'MKXPTeff',...
                        '.mat');
    dat = load(fname);
    lab = sprintf('m_{Kic} = %0.3f', MKXslope);

    [T_all, Y_all] = convert_dat(dat);
    ends_MKXPTeff(2,ii+1) = Y_all(end,2)/dat.pars.V_plasma;

    % plot result
    plot(T_all, Y_all(:,2)/dat.pars.V_plasma, ...
        'linewidth', lw, 'linestyle', ls, 'color', cmap(ii,:), ...
        'DisplayName',lab);
end
% plot PT GFR result without MKX
fname = './MultiDaySim/03-Oct-2023_driver_multiday_insulin-1_Kamt_meal-104_highKeff-1_ndays-50_notes-PTGFRonly.mat';
dat = load(fname);
lab = 'no MKX';
[T_all, Y_all] = convert_dat(dat);
plot(T_all, Y_all(:,2)/dat.pars.V_plasma, ...
    'linewidth', lw, 'linestyle', ls, 'color', cmap2(2,:), ...
    'DisplayName', 'PT + TGF effect (no MKX)');
ends_MKXPTeff(2, 1) = Y_all(end,2)/dat.pars.V_plasma;
legend('fontsize', f.leg, 'location', 'northwest')
yline(3.5,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray, 'HandleVisibility', 'off')
yline(5.0,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray, 'HandleVisibility', 'off')
set(gca, 'fontsize', f.gca)
xlabel('Time (days)', 'fontsize', f.xlab)
ylabel('Plasma [K^+] (mmol/L)', 'fontsize', f.ylab)
title('PT + TGF effect with MKX', 'fontsize', f.title)
grid on
ylim([3.5,8.0])

AddLetters2Plots(figure(2), {'(A)', '(B)'},...
                'HShift', -0.05, 'VShift', -0.06, ...
                'fontsize', f.labs)


%% Figure with end points
figure(3)
cmap = parula(6);
b1_c = cmap(1,:); b2_c = cmap(4,:);
clf;
subplot(1,2,1)
labs = {};
labs{1} = 'no MKX';
for ii = 2:(length(slope_vals)+1)
    labs{ii} = sprintf('m_{Kic} = %0.3f', slope_vals(ii-1));
end
x = 1:length(labs);
b = bar(x, [ends_MKX(1,:); ends_MKXPTeff(1,:)]);
b(1).FaceColor = b1_c;
b(2).FaceColor = b2_c;
xticklabels(labs)
yline(120,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray, 'HandleVisibility', 'off')
yline(140,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray, 'HandleVisibility', 'off')
set(gca, 'fontsize', f.gca)
legend({'No PT + TGF effects', 'with PT + TGF effects'})
grid on
ylabel('Final intracellular [K^+] (mmol/L)')
xlabel('m_{Kic}')

% K plas bar
subplot(1,2,2)
b = bar(x,[ends_MKX(2,:); ends_MKXPTeff(2,:)]);
b(1).FaceColor = b1_c;
b(2).FaceColor = b2_c;
xticklabels(labs)
yline(3.5,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray, 'HandleVisibility', 'off')
yline(5.0,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray, 'HandleVisibility', 'off')
set(gca, 'fontsize', f.gca)
legend({'No PT + TGF effects', 'with PT + TGF effects'})
grid on
ylabel('Final plasma [K^+] (mmol/L)')
xlabel('m_{Kic}')

AddLetters2Plots(figure(3), {'(A)', '(B)'},...
                'HShift', -0.05, 'VShift', -0.06, ...
                'fontsize', f.labs)

%---------
% functions
%---------
function [T_all, Y_all] = convert_dat(dat)
    % Put all data together
    T_all = []; Y_all = [];
    for i = 1:dat.n_days
        temp = dat.Tvals{i}./60 + 24 * (i - 1);
        T_all = [T_all; temp];
        Y_all = [Y_all; dat.Yvals{i}];
    end
    % convert to days
    T_all = T_all./24;
end