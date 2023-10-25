% Plot the results from the eta_PTKreab changing simulations
clear all;
eta_PTKreab_vals = [0.18,0.24,0.30,0.36,0.43,0.49,0.55,0.61,0.67]; %round(linspace(0.5 * 0.36, 0.67, 9),2);
% Parameters needed
V_plasma = 4.5;
V_muscle = 24;

% Get data
Tvals = cell(length(eta_PTKreab_vals), 1);
Yvals = cell(length(eta_PTKreab_vals), 1);
labs = cell(length(eta_PTKreab_vals), 1);
% values at end of simulation
Kend = zeros(2, length(eta_PTKreab_vals)); % row 1: K plasma end point; row 2: K_IC end point
for ii = 1:length(eta_PTKreab_vals)
    eta_ptKreab = eta_PTKreab_vals(ii);
    notes = 'etaPTKreab';
    n_days = 50;
    date_fixed = '25-Oct-2023'; % CHANGE DATE IF UPDATING DATA
    fname = strcat('./MultiDaySim/', date_fixed, '_driver_multiday',...
            '_insulin-', num2str(1),...
            '_Kamt_meal-', num2str(104),... % High K intake
            '_TGFeff-', num2str(1),... % Both PT + GFR
            '_alphaTGF-', num2str(0.11694),... % Baseline alpha_TGF
            '_etaPTKreab-', num2str(eta_ptKreab),...
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

    % set label
%     if eta_ptKreab == 0.36
%         labs{ii} = strcat('\eta_{pt-Kreab} =  ', ' ', num2str(eta_ptKreab,2));
%     elseif eta_ptKreab == 0.67
%         labs{ii} = strcat('\eta_{pt-Kreab} =  ', ' ', num2str(eta_ptKreab,2));
%     else
%         labs{ii} = strcat('\eta_{pt-Kreab} =  ', ' ', num2str(eta_ptKreab,2));
%     end
    
    %labs{ii} = strcat('\eta_{pt-Kreab} =  ', ' ', num2str(eta_ptKreab,2));
    labs{ii} = strcat('\eta_{pt} = ', num2str(eta_ptKreab,2));
end





%% Make figure
% Full simulation plot
figure(5)
clf;
nr = 2; nc =2;
cmap = turbo(length(eta_PTKreab_vals) + 2);
lw = 3;  ls = '-';
cgraymap = gray(5); cgray = cgraymap(2,:);
lwgray = 3; lsgray = ':';
f.labs = 18; f.xlab = 18; f.ylab = 16; f.gca = 16; f.leg = 14; f.title = 22;

% K_plasma
subplot(nr,nc,1)
hold on
for ii = 1:length(eta_PTKreab_vals)
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
%title('Plasma [K^+] for changing \eta_{pt-Kreab}', 'fontsize', f.title)
grid on
legend(labs, 'location', 'best', 'fontsize', f.leg)

% K_intracellular
subplot(nr,nc,2)
hold on
for ii = 1:length(eta_PTKreab_vals)
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
%title('Intracellular [K^+] for changing \eta_{pt-Kreab}', 'fontsize', f.title)
grid on
legend(labs, 'location', 'best', 'fontsize', f.leg)


% Only the points
ms = '^';
lw = 2.0;
marksize = 15;
xlims = [0.15, 0.7];
% K_plasma
subplot(nr,nc,3)
hold on
for ii = 1:length(eta_PTKreab_vals)
    plot(round(eta_PTKreab_vals(ii),2), Kend(1,ii), 'marker', ms, ...
                    'markerfacecolor', cmap(ii,:), 'color', cmap(ii,:), 'linewidth', lw, ...
                    'markersize', marksize)
end
xticks(round(eta_PTKreab_vals,2))
yline(3.5,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray, 'HandleVisibility', 'off')
yline(5.0,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray, 'HandleVisibility', 'off')
set(gca, 'fontsize', f.gca)
xlim(xlims)
ylim([3.0,8.0])
xlabel("\eta_{pt-Kreab}", 'fontsize', f.xlab)
ylabel('Final Plasma [K^+] (mmol/L)')
grid on

% K_intracellular
subplot(nr,nc,4)
hold on
for ii = 1:length(eta_PTKreab_vals)
    plot(round(eta_PTKreab_vals(ii),2), Kend(2,ii), 'marker', ms, ...
                    'markerfacecolor', cmap(ii,:), 'color', cmap(ii,:), 'linewidth', lw, ...
                    'markersize', marksize)
end
%xticklabels(labs)
xticks(round(eta_PTKreab_vals,2))
yline(120,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray, 'HandleVisibility', 'off')
yline(140,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray, 'HandleVisibility', 'off')
set(gca, 'fontsize', f.gca)
xlim(xlims)
xlabel("\eta_{pt-Kreab}", 'fontsize', f.xlab)
ylabel('Final Intracellular [K^+] (mmol/L)')
grid on

AddLetters2Plots(figure(5), {'(A1)', '(B1)', '(A2)', '(B2)'},...
                'HShift', -0.05, ...
                'VShift', -0.06, ...
                'fontsize', f.labs)









