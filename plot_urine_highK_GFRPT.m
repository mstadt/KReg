% Plot urine and kidney results from driver_multiday.m
clear all;

%% load data
f1 = './MultiDaySim/22-Nov-2023_driver_multiday_insulin-1_Kamt_meal-104_TGFeff-0_alphaTGF-0.11694_etaPTKreab-0.67_ndays-50_notes-noeffect_new.mat';
f2 = './MultiDaySim/22-Oct-2023_driver_multiday_insulin-1_Kamt_meal-104_TGFeff-1_alphaTGF-0.11694_etaPTKreab-0.36_ndays-50_notes-TongHighK.mat';
f3 = './MultiDaySim/22-Oct-2023_driver_multiday_insulin-1_Kamt_meal-104_TGFeff-3_alphaTGF-0.11694_etaPTKreab-0.36_ndays-50_notes-PTonly.mat';
f4 = './MultiDaySim/22-Oct-2023_driver_multiday_insulin-1_Kamt_meal-26_TGFeff-3_alphaTGF-0.11694_etaPTKreab-0.67_ndays-50_notes-control.mat';

dat1 = load(f1);
dat2 = load(f2);
dat3 = load(f3);
dat4 = load(f4);


lab1 = 'High K^+ - no PT/TGF effect';
lab2 = 'High K^+ - PT + TGF effects';
lab3 = 'High K^+ - only PT effect'; 
lab4 = 'Control K^+';

%% Convert data to days
T_all1 = []; Y_all1 = []; Urine1 = []; Urine1_frac = []; FilK1 = [];
T_all2 = []; Y_all2 = []; Urine2 = []; Urine2_frac = []; FilK2 = [];
T_all3 = []; Y_all3 = []; Urine3 = []; Urine3_frac = []; FilK3 = [];
T_all4 = []; Y_all4 = []; Urine4 = []; Urine4_frac = []; FilK4 = [];

% Mean Urine per day
Urine1_stats = zeros(dat1.n_days, 6); % col1: mean, col2: median, col3: max, col4: min, col5: cumulative
Urine2_stats = zeros(dat1.n_days, 6);
Urine3_stats = zeros(dat1.n_days, 6);
Urine4_stats = zeros(dat1.n_days, 6);

FilK1_stats = zeros(dat1.n_days, 6); % col1: mean, col2: median, col3: max, col4: min, col5: cumulative
FilK2_stats = zeros(dat1.n_days, 6);
FilK3_stats = zeros(dat1.n_days, 6);
FilK4_stats = zeros(dat1.n_days, 6);

UrineFrac1_mean = zeros(dat1.n_days,1);
UrineFrac2_mean = zeros(dat1.n_days,1);
UrineFrac3_mean = zeros(dat1.n_days,1);
UrineFrac4_mean = zeros(dat1.n_days,1);
for ii = 1:dat1.n_days
    % Data set 1
    Y = dat1.Yvals{ii};
    v = compute_kidney_vars(Y, dat1.params, ...
                                            'do_MKX', [dat1.doMKX, dat1.MKXslope], ...
                                            'do_FF', dat1.doFF, ...
                                            'TGF_eff', [dat1.TGF_eff, ...
                                                        dat1.alpha_TGF, ...
                                                        dat1.eta_ptKreab]);
    Urine1 = [Urine1; v.UrineK];
    
    Urine_frac = v.UrineK ./ v.filK;
    UrineFrac1_mean(ii) = mean(Urine_frac);
    Urine1_frac = [Urine1_frac; Urine_frac];

    % Compute cumulative urine
    x = dat1.Tvals{ii};
    y = v.UrineK;
    cum_UK = trapz(x,y); % cumulative urine

    FilK1 = [FilK1; v.filK];
    y = v.filK;
    cum_filK = trapz(x,y); % cumulative filK


    FilK1_stats(ii,:) = [mean(v.filK), ...
                            median(v.filK), ...
                            max(v.filK), ...
                            min(v.filK), ...
                            std(v.filK),...
                            cum_filK];

    Urine1_stats(ii,:) = [mean(v.UrineK), ...
                            median(v.UrineK), ...
                            max(v.UrineK), ...
                            min(v.UrineK), ...
                            std(v.UrineK),...
                            cum_UK];

    temp = dat1.Tvals{ii}./60 + (24)*(ii - 1);
    T_all1 = [T_all1; temp];
    Y_all1 = [Y_all1; dat1.Yvals{ii}];

    % Data set 2
    Y = dat2.Yvals{ii};
    v = compute_kidney_vars(Y, dat2.params, ...
                                        'do_MKX', [dat2.doMKX, dat2.MKXslope], ...
                                        'do_FF', dat2.doFF, ...
                                        'TGF_eff', [dat2.TGF_eff, ...
                                                    dat2.alpha_TGF, ...
                                                    dat2.eta_ptKreab]);
    Urine2 = [Urine2; v.UrineK];
    Urine_frac = v.UrineK ./ v.filK;
    Urine2_frac = [Urine2_frac; Urine_frac];
    UrineFrac2_mean(ii) = mean(Urine_frac);
    % Compute cumulative urine
    x = dat2.Tvals{ii};
    y = v.UrineK;
    cum_UK = trapz(x,y); % cumulative urine
    
    FilK2 = [FilK2; v.filK];
    y = v.filK;
    cum_filK = trapz(x,y); % cumulative filK


    FilK2_stats(ii,:) = [mean(v.filK), ...
                            median(v.filK), ...
                            max(v.filK), ...
                            min(v.filK), ...
                            std(v.filK),...
                            cum_filK];


    Urine2_stats(ii,:) = [mean(v.UrineK), ...
                            median(v.UrineK), ...
                            max(v.UrineK), ...
                            min(v.UrineK),...
                            std(v.UrineK),...
                            cum_UK];

    temp = dat2.Tvals{ii}./60 + 24*(ii - 1);
    T_all2 = [T_all2; temp];
    Y_all2 = [Y_all2; dat2.Yvals{ii}];

    % Data set 3
    Y = dat3.Yvals{ii};
    v = compute_kidney_vars(Y, dat3.params, ...
                                        'do_MKX', [dat3.doMKX, dat3.MKXslope], ...
                                        'do_FF', dat3.doFF, ...
                                        'TGF_eff', [dat3.TGF_eff, ...
                                                    dat3.alpha_TGF, ...
                                                    dat3.eta_ptKreab]);
    Urine3 = [Urine3; v.UrineK];
    Urine_frac = v.UrineK ./ v.filK;
    Urine3_frac = [Urine3_frac; Urine_frac];
    UrineFrac3_mean(ii) = mean(Urine_frac);
    % Compute cumulative urine
    x = dat3.Tvals{ii};
    y = v.UrineK;
    cum_UK = trapz(x,y); % cumulative urine
    FilK3 = [FilK3; v.filK];
    y = v.filK;
    cum_filK = trapz(x,y); % cumulative filK


    FilK3_stats(ii,:) = [mean(v.filK), ...
                            median(v.filK), ...
                            max(v.filK), ...
                            min(v.filK), ...
                            std(v.filK),...
                            cum_filK];
    Urine3_stats(ii,:) = [mean(v.UrineK), ...
                            median(v.UrineK), ...
                            max(v.UrineK), ...
                            min(v.UrineK)...
                            std(v.UrineK),...
                            cum_UK];
    temp = dat3.Tvals{ii}./60 + 24 * (ii - 1);
    T_all3 = [T_all3; temp];
    Y_all3 = [Y_all3; dat3.Yvals{ii}];

    % Data set 4
    Y = dat4.Yvals{ii};
    v = compute_kidney_vars(Y, dat4.params, ...
                                        'do_MKX', [dat4.doMKX, dat4.MKXslope], ...
                                        'do_FF', dat4.doFF, ...
                                        'TGF_eff', [dat4.TGF_eff, ...
                                                    dat4.alpha_TGF, ...
                                                    dat4.eta_ptKreab]);
    Urine4 = [Urine4; v.UrineK];
    Urine_frac = v.UrineK ./ v.filK;
    Urine4_frac = [Urine4_frac; Urine_frac];
    UrineFrac4_mean(ii) = mean(Urine_frac);
    % Compute cumulative urine
    x = dat4.Tvals{ii};
    y = v.UrineK;
    cum_UK = trapz(x,y); % cumulative urine
    FilK4 = [FilK4; v.filK];
    y = v.filK;
    cum_filK = trapz(x,y); % cumulative filK


    FilK4_stats(ii,:) = [mean(v.filK), ...
                            median(v.filK), ...
                            max(v.filK), ...
                            min(v.filK), ...
                            std(v.filK),...
                            cum_filK];
    Urine4_stats(ii,:) = [mean(v.UrineK), ...
                            median(v.UrineK), ...
                            max(v.UrineK), ...
                            min(v.UrineK)...
                            std(v.UrineK),...
                            cum_UK];

    temp = dat4.Tvals{ii}./60 + 24 * (ii - 1);
    T_all4 = [T_all4; temp];
    Y_all4 = [Y_all4; dat4.Yvals{ii}];

end
% convert to days
T_all1 = T_all1./24; T_all2 = T_all2./24; 
T_all3 = T_all3./24; T_all4 = T_all4./24;

%% Make figure
figure(10)
clf;
clf;
nr = 1; nc = 2;
f.labs = 18; f.xlab = 18; f.ylab = 18; f.gca = 18; f.leg = 16; f.title = 22;
lw = 3; lwgray = 3; lsgray = ':';
ls1 = '-'; ls2 = '-'; ls3 = '-'; ls4 = '-'; 
cmap = parula(3);
cmap2 = spring(3);
c1 = cmap2(1,:); c2 = cmap(1,:); 
c3 = cmap(2,:);c4 = cmap2(2,:);
cgraymap = gray(5);
cgray = cgraymap(2,:);

hold on
plot(T_all1, Urine1, 'linewidth',lw,'linestyle',ls1,'color',c1)
plot(T_all2, Urine2, 'linewidth',lw,'linestyle',ls2,'color',c2)
plot(T_all3, Urine3, 'linewidth',lw,'linestyle',ls3,'color',c3)
plot(T_all4, Urine4, 'linewidth',lw,'linestyle',ls4,'color',c4)
set(gca, 'fontsize', f.gca)
xlabel('Time (days)', 'fontsize', f.xlab)
ylabel('Urine K^+ excretion (mmol/min)', 'fontsize', f.ylab)
title('Urine K^+ excretion', 'fontsize', f.title)
grid on

legend({lab1, lab2,lab3,lab4}, 'fontsize', f.leg, 'location', 'northwest')


%% Bar plot final 
cg = gray(6);
cgray = cg(3,:);
UK_frac_means = [mean(Urine1_frac*100), ...
                        mean(Urine2_frac*100), ...
                        mean(Urine3_frac*100),...
                        mean(Urine4_frac*100)];
UK_frac_stds = [std(Urine1_frac*100),...
                        std(Urine2_frac*100), ...
                        std(Urine3_frac*100),...
                        std(Urine4_frac*100)];
xvals=  [1,2,3,4];
figure(25)
clf;
subplot(1,2,1)
hold on
% Simulation data

% High K --  PT + TGF  effects
b1 = bar(xvals(1), ...
        UK_frac_means(2), ...
    'FaceColor', c2);
% High K - only PT effect
b2 = bar(xvals(2),...
    UK_frac_means(3),...
    'FaceColor', c3);
% High K+ - no PT/TGF effect
b3 = bar(xvals(3),...
            UK_frac_means(1),...
            'FaceColor', c1);
b4 = bar(xvals(4),...
            UK_frac_means(4),...
            'FaceColor', c4);
errorbar(xvals, UK_frac_means([2,3,1,4]), UK_frac_stds([2,3,1,4]), ...
                'linestyle','none', ...
                'linewidth', 3, ...
                'color', cgray)



set(gca, 'fontsize',f.gca)
ylabel('Fractional urine K^+ excretion (%)', 'fontsize',f.ylab)
ylim([0.0,100])
xlim([0.5,4.5])

xticks([1,2,3,4])
labs = {lab2, lab3,...
    lab1,lab4};
xticklabels(labs)
grid on

% Bar plot -- Cumulative urine K
% Cumulative urine K average
UK_means = [mean(Urine1_stats(:,6)), ...
                        mean(Urine2_stats(:,6)), ...
                        mean(Urine3_stats(:,6)),...
                        mean(Urine4_stats(:,6))];
UK_stds = [std(Urine1_stats(:,6)),...
                        std(Urine2_stats(:,6)), ...
                        std(Urine3_stats(:,6)),...
                        std(Urine4_stats(:,6))];
xvals=  [1,2,3,4,5];
subplot(1,2,2)
hold on
% Simulation data

% High K --  PT + TGF  effects
b1 = bar(xvals(1), ...
        UK_means(2), ...
    'FaceColor', c2);
% High K - only PT effect
b2 = bar(xvals(2),...
    UK_means(3),...
    'FaceColor', c3);
% High K+ - no PT/TGF effect
b3 = bar(xvals(3),...
            UK_means(1),...
            'FaceColor', c1);
b4 = bar(xvals(4),...
            UK_means(4),...
            'FaceColor', c4);
% errorbar(xvals, UK_means([2,3,1,4,5]), UK_stds([2,3,1,4,5]), ...
%                 'linestyle','none', ...
%                 'linewidth', 3, ...
%                 'color', cgray)



set(gca, 'fontsize',f.gca)
ylabel('Cumulative urine K^+ (mmol/24 hr)', 'fontsize',f.ylab)
% ylim([0.0,100])
xlim([0.5,4.5])

xticks([1,2,3,4])
labs = {lab2, lab3,...
    lab1,lab4};
xticklabels(labs)
grid on

AddLetters2Plots(figure(25), {'(A)', '(B)'},...
                'HShift', -0.05, 'VShift', -0.06, ...
                'fontsize', f.labs)


%% Filtered K
% Bar plot -- Cumulative filtered K
% Cumulative filterd K average
filK_means = [mean(FilK1_stats(:,6)), ...
            mean(FilK2_stats(:,6)), ...
            mean(FilK3_stats(:,6)),...
            mean(FilK4_stats(:,6))];
filK_stds = [std(FilK1_stats(:,6)),...
                        std(FilK2_stats(:,6)), ...
                        std(FilK3_stats(:,6)),...
                        std(FilK4_stats(:,6))];
xvals=  [1,2,3,4];
figure(50)
clf;
subplot(1,2,2)
hold on
% Simulation data

% High K --  PT + TGF  effects
b1 = bar(xvals(1), ...
        filK_means(2), ...
    'FaceColor', c2);
% High K - only PT effect
b2 = bar(xvals(2),...
    filK_means(3),...
    'FaceColor', c3);
% High K+ - no PT/TGF effect
b3 = bar(xvals(3),...
            filK_means(1),...
            'FaceColor', c1);
% control
b4 = bar(xvals(4),...
            filK_means(4),...
            'FaceColor', c4);
% errorbar(xvals, filK_means([2,3,1,4]), filK_stds([2,3,1,4]), ...
%                 'linestyle','none', ...
%                 'linewidth', 3, ...
%                 'color', cgray)


set(gca, 'fontsize',f.gca)
ylabel('Cumulative filtered K^+ (mmol/24 hr)', 'fontsize',f.ylab)
xlim([0.5,4.5])

xticks([1,2,3,4])
labs = {lab2, lab3,...
    lab1,lab4};
xticklabels(labs)
grid on

% Filtered K all
% figure(51)
% clf;
subplot(1,2,1)
nr = 1; nc = 2;
f.labs = 18; f.xlab = 18; f.ylab = 18; f.gca = 18; f.leg = 16; f.title = 22;
lw = 3; lwgray = 3; lsgray = ':';
ls1 = '-'; ls2 = '-'; ls3 = '-'; ls4 = '-'; 
cmap = parula(3);
cmap2 = spring(3);
c1 = cmap2(1,:); c2 = cmap(1,:); 
c3 = cmap(2,:);c4 = cmap2(2,:);
cgraymap = gray(5);
cgray = cgraymap(2,:);

hold on
plot(T_all1, FilK1, 'linewidth',lw,'linestyle',ls1,'color',c1)
plot(T_all2, FilK2, 'linewidth',lw,'linestyle',ls2,'color',c2)
plot(T_all3, FilK3, 'linewidth',lw,'linestyle',ls3,'color',c3)
plot(T_all4, FilK4, 'linewidth',lw,'linestyle',ls4,'color',c4)
set(gca, 'fontsize', f.gca)
xlabel('Time (days)', 'fontsize', f.xlab)
ylabel('Filtered K^+ (mmol/min)', 'fontsize', f.ylab)
ylim([0.3,1.1])
%title('Urine K^+ excretion', 'fontsize', f.title)
grid on

legend({lab1, lab2,lab3,lab4}, 'fontsize', f.leg, 'location', 'northwest')

AddLetters2Plots(figure(50), {'(A)', '(B)'},...
                'HShift', -0.05, 'VShift', -0.06, ...
                'fontsize', f.labs)