% Plot the results from alpha_TGF changing simulations
clear all;

% Fold changes of alpha_TGF parameter
alpha_fold = [0.5000, 0.6250, 0.7500, 0.8750, 1.0000, 1.1250, 1.2500, 1.3750, 1.5000]; 

% Parameters needed
pars = set_params();
V_plasma = pars.V_plasma;
V_muscle = pars.V_muscle;


%% Get data
Tvals = cell(length(alpha_fold), 1);
Yvals = cell(length(alpha_fold), 1);
labs = cell(length(alpha_fold), 1);
% values at end of simulation
Kend = zeros(2, length(alpha_fold)); % row 1: Kplas, row 2: KIC
for ii = 1:length(alpha_fold)
    fold_val = alpha_fold(ii);

    alpha_TGF = fold_val * pars.alpha_TGF;

    notes = 'alpha_tgf';
    n_days = 50;
    date_fixed = '25-Oct-2023'; % change date if updating data

    fname = strcat('./MultiDaySim/', date_fixed, '_driver_multiday',...
            '_insulin-', num2str(1),...
            '_Kamt_meal-', num2str(104),... % High K intake
            '_TGFeff-', num2str(1),... % Both PT + GFR
            '_alphaTGF-', num2str(alpha_TGF),... % change alpha_TGF
            '_etaPTKreab-', num2str(0.36),...
            '_ndays-', num2str(n_days),...
            '_notes-', notes,...
            '.mat');
    dat = load(fname);

    % Get T and Y values
    T = [];
    Y = [];
    for jj = 1:n_days
        temp = dat.Tvals{jj}./60 + (24)*(jj - 1);
        T = [T; temp];
        Y = [Y; dat.Yvals{jj}];
    end

    % save endpoint
    Kend(1, ii) = Y(end, 2) / V_plasma; % K_plas
    Kend(2, ii) = Y(end, 4) / V_muscle; % K_IC
 
    % save T values and Y values
    Tvals{ii} = T./24; % convert to days
    Yvals{ii} = Y;

    if fold_val == 1.0
        labs{ii} = strcat('baseline \alpha_{TGF}');
    else
        labs{ii} = strcat(num2str(fold_val), ' x \alpha_{TGF}^{base}');
    end
end

%% Make figure
% Full simulation plot
figure(5)
clf;
nr = 2; nc =2;
cmap = turbo(length(alpha_fold) + 2);
lw = 3;  ls = '-';
cgraymap = gray(5); cgray = cgraymap(1,:);
lwgray = 4.5; lsgray = ':';
f.labs = 18; f.xlab = 18; f.ylab = 16; f.gca = 16; f.leg = 14; f.title = 22;

ylims_plas = [3.4,6.5];
ylims_muscle = [110,200];

% K_plasma
subplot(nr,nc,1)
hold on
for ii = 1:length(alpha_fold)
    T = Tvals{ii};
    Y = Yvals{ii};
    plot(T, Y(:,2)/V_plasma, 'linewidth', lw, 'linestyle', ls, 'color', cmap(ii,:))
end
yline(3.5, 'color', cgray, 'linestyle', lsgray, 'linewidth', lwgray, 'HandleVisibility', 'off')
yline(5.0, 'color', cgray, 'linestyle', lsgray, 'linewidth', lwgray, 'HandleVisibility', 'off')
set(gca, 'fontsize', f.gca)
%legend(labs, 'location', 'best', 'fontsize', f.leg)
xlabel('Time (days)', 'fontsize', f.xlab)
ylabel('Plasma [K^+] (mmol/L)', 'fontsize', f.ylab)
ylim(ylims_plas)
%title('Plasma [K^+] for changing \eta_{pt-Kreab}', 'fontsize', f.title)
grid on
legend(labs, 'location', 'best', 'fontsize', f.leg)

% K_intracellular
subplot(nr,nc,2)
hold on
for ii = 1:length(alpha_fold)
    T = Tvals{ii};
    Y = Yvals{ii};
    plot(T, Y(:,4)/V_muscle, 'linewidth', lw, 'linestyle', ls, 'color', cmap(ii,:))
end
yline(120, 'color', cgray, 'linestyle', lsgray, 'linewidth', lwgray)
yline(140, 'color', cgray, 'linestyle', lsgray, 'linewidth', lwgray)
set(gca, 'fontsize', f.gca)
legend(labs, 'location', 'best', 'fontsize', f.leg)
xlabel('Time (days)', 'fontsize', f.xlab)
ylabel('Intracellular [K^+] (mmol/L)', 'fontsize', f.ylab)
ylim(ylims_muscle)
%title('Intracellular [K^+] for changing \eta_{pt-Kreab}', 'fontsize', f.title)
grid on
legend(labs, 'location', 'best', 'fontsize', f.leg)


% Only the points
ms = '^';
lw = 2.0;
marksize = 15;
% K_plasma
subplot(nr,nc,3)
hold on
for ii = 1:length(alpha_fold)
    plot(alpha_fold(ii), Kend(1,ii), 'marker', ms, ...
                    'markerfacecolor', cmap(ii,:), 'color', cmap(ii,:), 'linewidth', lw, ...
                    'markersize', marksize)
end
xticks(alpha_fold)
yline(3.5,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray, 'HandleVisibility', 'off')
yline(5.0,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray, 'HandleVisibility', 'off')
set(gca, 'fontsize', f.gca)
xlabel("\alpha_{TGF} / \alpha_{TGF}^{base}", 'fontsize', f.xlab)
ylabel('Final Plasma [K^+] (mmol/L)')
ylim(ylims_plas)
grid on

% K_intracellular
subplot(nr,nc,4)
hold on
for ii = 1:length(alpha_fold)
    plot(alpha_fold(ii), Kend(2,ii), 'marker', ms, ...
                    'markerfacecolor', cmap(ii,:), 'color', cmap(ii,:), 'linewidth', lw, ...
                    'markersize', marksize)
end
%xticklabels(labs)
xticks(alpha_fold)
yline(120,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray, 'HandleVisibility', 'off')
yline(140,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray, 'HandleVisibility', 'off')
set(gca, 'fontsize', f.gca)
xlabel("\alpha_{TGF} / \alpha_{TGF}^{base}", 'fontsize', f.xlab)
ylabel('Final Intracellular [K^+] (mmol/L)')
ylim(ylims_muscle)
grid on

AddLetters2Plots(figure(5), {'(A1)', '(B1)', '(A2)', '(B2)'},...
                'HShift', -0.05, ...
                'VShift', -0.06, ...
                'fontsize', f.labs)