% Plot urine and kidney results from driver_multiday.m
clear all;

%% load data
f1 = './MultiDaySim/22-Oct-2023_driver_multiday_insulin-1_Kamt_meal-104_TGFeff-1_alphaTGF-0.11694_etaPTKreab-0.67_ndays-50_notes-baseline.mat';
f2 = './MultiDaySim/22-Oct-2023_driver_multiday_insulin-1_Kamt_meal-104_TGFeff-1_alphaTGF-0.11694_etaPTKreab-0.36_ndays-50_notes-TongHighK.mat';
f3 = './MultiDaySim/22-Oct-2023_driver_multiday_insulin-1_Kamt_meal-104_TGFeff-3_alphaTGF-0.11694_etaPTKreab-0.36_ndays-50_notes-PTonly.mat';
f4 = './MultiDaySim/22-Oct-2023_driver_multiday_insulin-1_Kamt_meal-104_TGFeff-2_alphaTGF-0.11694_etaPTKreab-0.36_ndays-50_notes-GFRonly.mat';
f5 = './MultiDaySim/22-Oct-2023_driver_multiday_insulin-1_Kamt_meal-26_TGFeff-3_alphaTGF-0.11694_etaPTKreab-0.67_ndays-50_notes-control.mat';

dat1 = load(f1);
dat2 = load(f2);
dat3 = load(f3);
dat4 = load(f4);
dat5 = load(f5);


lab1 = 'High K^+ - no PT/TGF effect';
lab2 = 'High K^+ - PT + TGF effects';
lab3 = 'High K^+ - only PT effect'; 
lab4 = 'High K^+ - only TGF effect';
lab5 = 'Control K^+';

%% Convert data to days
T_all1 = []; Y_all1 = []; Urine1 = []; Urine1_frac = [];
T_all2 = []; Y_all2 = []; Urine2 = []; Urine2_frac = [];
T_all3 = []; Y_all3 = []; Urine3 = []; Urine3_frac = [];
T_all4 = []; Y_all4 = []; Urine4 = []; Urine4_frac = [];
T_all5 = []; Y_all5 = []; Urine5 = []; Urine5_frac = [];

% Mean Urine per day
Urine1_stats = zeros(dat1.n_days, 5); % col1: mean, col2: median, col3: max, col4: min
Urine2_stats = zeros(dat1.n_days, 5);
Urine3_stats = zeros(dat1.n_days, 5);
Urine4_stats = zeros(dat1.n_days, 5);
Urine5_stats = zeros(dat1.n_days, 5);

UrineFrac1_mean = zeros(dat1.n_days,1);
UrineFrac2_mean = zeros(dat1.n_days,1);
UrineFrac3_mean = zeros(dat1.n_days,1);
UrineFrac4_mean = zeros(dat1.n_days,1);
UrineFrac5_mean = zeros(dat1.n_days,1);
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
    Urine1_stats(ii,:) = [mean(v.UrineK), ...
                            median(v.UrineK), ...
                            max(v.UrineK), ...
                            min(v.UrineK), ...
                            std(v.UrineK)];

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
    Urine2_stats(ii,:) = [mean(v.UrineK), ...
                            median(v.UrineK), ...
                            max(v.UrineK), ...
                            min(v.UrineK),...
                            std(v.UrineK)];

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
    Urine3_stats(ii,:) = [mean(v.UrineK), ...
                            median(v.UrineK), ...
                            max(v.UrineK), ...
                            min(v.UrineK)...
                            std(v.UrineK)];
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
    Urine4_stats(ii,:) = [mean(v.UrineK), ...
                            median(v.UrineK), ...
                            max(v.UrineK), ...
                            min(v.UrineK)...
                            std(v.UrineK)];

    temp = dat4.Tvals{ii}./60 + 24 * (ii - 1);
    T_all4 = [T_all4; temp];
    Y_all4 = [Y_all4; dat4.Yvals{ii}];

    % Data set 5
    Y = dat5.Yvals{ii};
    v = compute_kidney_vars(Y, dat5.params, ...
                                        'do_MKX', [dat5.doMKX, dat5.MKXslope], ...
                                        'do_FF', dat5.doFF, ...
                                        'TGF_eff', [dat5.TGF_eff, ...
                                                    dat5.alpha_TGF, ...
                                                    dat5.eta_ptKreab]);
    Urine5 = [Urine5; v.UrineK];
    Urine_frac = v.UrineK ./ v.filK;
    UrineFrac5_mean(ii) = mean(Urine_frac);
    Urine5_frac = [Urine5_frac; Urine_frac];
    Urine5_stats(ii,:) = [mean(v.UrineK), ...
                            median(v.UrineK), ...
                            max(v.UrineK), ...
                            min(v.UrineK)...
                            std(v.UrineK)];

    temp = dat5.Tvals{ii}./60 + 24 * (ii - 1);
    T_all5 = [T_all5; temp];
    Y_all5 = [Y_all5; dat5.Yvals{ii}];
end
% convert to days
T_all1 = T_all1./24; T_all2 = T_all2./24; 
T_all3 = T_all3./24; T_all4 = T_all4./24;
T_all5 = T_all5./24;

%% Make figure
figure(10)
clf;
f.labs = 18; f.xlab = 18; f.ylab = 18; f.gca = 18; f.leg = 16; f.title = 22;
lw = 3; lwgray = 3; lsgray = ':';
ls1 = '-'; ls2 = '-'; ls3 = '-'; ls4 = '-'; ls5 = '-';
cmap = parula(4);
cmap2 = spring(3);
c1 = cmap2(1,:); c2 = cmap(1,:); 
c3 = cmap(2,:);c4 = cmap(3,:);c5 = cmap2(2,:);
cgraymap = gray(5);
cgray = cgraymap(2,:);

hold on
plot(T_all1, Urine1, 'linewidth',lw,'linestyle',ls1,'color',c1)
plot(T_all2, Urine2, 'linewidth',lw,'linestyle',ls2,'color',c2)
plot(T_all3, Urine3, 'linewidth',lw,'linestyle',ls3,'color',c3)
plot(T_all4, Urine4, 'linewidth',lw,'linestyle',ls4,'color',c4)
plot(T_all5, Urine5, 'linewidth',lw,'linestyle',ls5,'color',c5)
set(gca, 'fontsize', f.gca)
xlabel('Time (days)', 'fontsize', f.xlab)
ylabel('Urine K^+ excretion (mmol/min)', 'fontsize', f.ylab)
title('Urine K^+ excretion', 'fontsize', f.title)
grid on

legend({lab1, lab2,lab3,lab4,lab5}, 'fontsize', f.leg, 'location', 'northwest')

%%
% Median Urine K
figure(11)
clf;
marksize = 15;
day_vals = 1:50;
lw = 2; cs = 10;
hold on
% plot(day_vals,Urine1_stats(:,2), 'marker', 'o', 'markersize',marksize,...
%                             'linestyle', 'none', 'color', c1,...
%                             'markerfacecolor', c1)
errorbar(day_vals, Urine1_stats(:,2), Urine1_stats(:,5), 'color',c1, ...
                                'linestyle','none', ...
                                'marker', 'o', 'markerfacecolor', c1,...
                                'markersize', marksize, ...
                                'linewidth', lw, 'capsize', cs)
errorbar(day_vals, Urine2_stats(:,2), Urine2_stats(:,5), 'marker', 'o', ...
                            'markersize',marksize,...
                            'linestyle', 'none', 'color', c2,...
                            'markerfacecolor', c2, ...
                            'linewidth', lw, 'capsize', cs)
errorbar(day_vals, Urine3_stats(:,2), Urine3_stats(:,5), ...
    'marker', 'o', 'markersize',marksize,...
                            'linestyle', 'none', 'color', c3,...
                            'markerfacecolor', c3, ...
                            'linewidth', lw, 'capsize', cs)
errorbar(day_vals, Urine4_stats(:,2), Urine4_stats(:,5),...
    'marker', 'o', 'markersize',marksize,...
                            'linestyle', 'none', 'color', c4,...
                            'markerfacecolor', c4, ...
                            'linewidth', lw, 'capsize', cs)
errorbar(day_vals, Urine5_stats(:,2), Urine5_stats(:,5),...
    'marker', 'o', 'markersize',marksize,...
                            'linestyle', 'none', 'color', c5,...
                            'markerfacecolor', c5, ...
                            'linewidth', lw, 'capsize', cs)
set(gca, 'fontsize', f.gca)
xlabel('Time (days)', 'fontsize', f.xlab)
ylabel('Mean Urine K^+ excretion (mmol/min)', 'fontsize', f.ylab)
title('Urine K^+ excretion', 'fontsize', f.title)
grid on

legend({lab1, lab2,lab3,lab4,lab5}, 'fontsize', f.leg, 'location', 'northwest')

%% Urine fractional excretion
figure(12)
clf;
f.labs = 18; f.xlab = 18; f.ylab = 18; f.gca = 18; f.leg = 16; f.title = 22;
lw = 3; lwgray = 3; lsgray = ':';
ls1 = '-'; ls2 = '-'; ls3 = '-'; ls4 = '-'; ls5 = '-';
cmap = parula(4);
cmap2 = spring(3);
c1 = cmap2(1,:); c2 = cmap(1,:); 
c3 = cmap(2,:);c4 = cmap(3,:);c5 = cmap2(2,:);
cgraymap = gray(5);
cgray = cgraymap(2,:);

hold on
plot(T_all1, Urine1_frac, 'linewidth',lw,'linestyle',ls1,'color',c1)
plot(T_all2, Urine2_frac, 'linewidth',lw,'linestyle',ls2,'color',c2)
plot(T_all3, Urine3_frac, 'linewidth',lw,'linestyle',ls3,'color',c3)
plot(T_all4, Urine4_frac, 'linewidth',lw,'linestyle',ls4,'color',c4)
plot(T_all5, Urine5_frac, 'linewidth',lw,'linestyle',ls5,'color',c5)
set(gca, 'fontsize', f.gca)
xlabel('Time (days)', 'fontsize', f.xlab)
ylabel('Fractional urine K^+ excretion (mmol/min)', 'fontsize', f.ylab)
title('Fractional urine K^+ excretion', 'fontsize', f.title)
grid on

legend({lab1, lab2,lab3,lab4,lab5}, 'fontsize', f.leg, 'location', 'northwest')


%% Box plots of the data
figure(20)
clf;
Urine_frac_all = [Urine1_frac;Urine2_frac;Urine3_frac;Urine4_frac;Urine5_frac];
x = Urine_frac_all;
labs = {lab1, lab2,lab3,lab4,lab5};
g1 = repmat({lab1},length(Urine1_frac),1);
g2 = repmat({lab2},length(Urine2_frac),1);
g3 = repmat({lab3},length(Urine3_frac),1);
g4 = repmat({lab4},length(Urine4_frac),1);
g5 = repmat({lab5},length(Urine5_frac),1);

g = [g1; g2; g3;g4;g5];
boxplot(x,g);
hold on
% b1 = boxplot(Urine1_frac, 'positions', 1); %boxplot(Urine1_frac, 'Labels', lab1);
% b2 = boxplot(Urine2_frac, 'positions', 2); %boxplot(Urine2_frac, 'Labels', lab2);
% b3 = boxplot(Urine3_frac, 'positions', 3);
% b4 = boxplot(Urine4_frac, 'positions', 4);
% b5 = boxplot(Urine5_frac, 'positions', 5);

set(gca, 'fontsize', f.gca)
ylabel('Fractional urine K^+ excretion')
grid on

%% Bar plot
UK_frac_means = [mean(Urine1_frac), mean(Urine2_frac), ...
                        mean(Urine3_frac),...
                        mean(Urine4_frac),...
                        mean(Urine5_frac)];
UK_frac_stds = [std(Urine1_frac), std(Urine2_frac), ...
                        std(Urine3_frac),...
                        std(Urine4_frac),...
                        std(Urine5_frac)];
xvals= 1:length(UK_frac_means);
figure(21)
clf;
hold on
b = bar(UK_frac_means);

errorbar(xvals, UK_frac_means, UK_frac_stds, ...
                'linestyle','none', ...
                'linewidth', 3, ...
                'color', 'black')

set(gca, 'fontsize',f.gca)
ylabel('Fractional urine K^+ excretion', 'fontsize',f.ylab)
ylim([0.0,1.0])
xlim([0.5,5.5])
labs = {lab1, lab2,lab3,lab4,lab5};
xticks(xvals)
xticklabels(labs)
grid on
%
%labs = {lab1, lab2,lab3,lab4,lab5};
%(labs)

%% Bar plot 2
cmap = parula(4);
cmap2 = spring(3);
c1 = cmap2(1,:); c2 = cmap(1,:); 
c3 = cmap(2,:);c4 = cmap(3,:);c5 = cmap2(2,:);
cg = gray(6);
cgray = cg(3,:);
UK_frac_means = [mean(Urine1_frac*100), ...
                        mean(Urine2_frac*100), ...
                        mean(Urine3_frac*100),...
                        mean(Urine4_frac*100),...
                        mean(Urine5_frac*100)];
UK_frac_stds = [std(Urine1_frac*100),...
                        std(Urine2_frac*100), ...
                        std(Urine3_frac*100),...
                        std(Urine4_frac*100),...
                        std(Urine5_frac*100)];
xvals=  [2,3,4,5,7];
figure(22)
clf;
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
b5 = bar(xvals(5), UK_frac_means(5), ...
    'FaceColor', c5);
errorbar(xvals, UK_frac_means([2,3,1,4,5]), UK_frac_stds([2,3,1,4,5]), ...
                'linestyle','none', ...
                'linewidth', 3, ...
                'color', cgray)

% Wang Data
xvals2 = [1,6];
Wang_mean = [80, 21]; % High K, Control FEk
Wang_std = [50, 24];
errorbar(xvals2, Wang_mean, Wang_std,...
                'linestyle', 'none', ...
                'linewidth', 3,...
                'marker', 'o',...
                'markerfacecolor', 'black',...
                'markersize', 10,...
                'color', 'black')

set(gca, 'fontsize',f.gca)
ylabel('Fractional urine K^+ excretion (%)', 'fontsize',f.ylab)
ylim([0.0,140])
xlim([0.5,7.5])

xticks([1,2,3,4,5,6,7])
labs = {'High K^+ - Wang et al., 2023',lab2, lab3,...
    lab1,lab4,'Control K^+ - Wang et al., 2023', lab5};
xticklabels(labs)
grid on

%% Bar plot 3
cmap = parula(4);
cmap2 = spring(3);
c1 = cmap2(1,:); c2 = cmap(1,:); 
c3 = cmap(2,:);c4 = cmap(3,:);c5 = cmap2(2,:);
cg = gray(6);
cgray = cg(3,:);
UK_frac_means = [mean(Urine1_frac*100), ...
                        mean(Urine2_frac*100), ...
                        mean(Urine3_frac*100),...
                        mean(Urine4_frac*100),...
                        mean(Urine5_frac*100)];
UK_frac_stds = [std(Urine1_frac*100),...
                        std(Urine2_frac*100), ...
                        std(Urine3_frac*100),...
                        std(Urine4_frac*100),...
                        std(Urine5_frac*100)];
xvals=  [1,2,3,4,5];
figure(23)
clf;
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
b5 = bar(xvals(5), UK_frac_means(5), ...
    'FaceColor', c5);
errorbar(xvals, UK_frac_means([2,3,1,4,5]), UK_frac_stds([2,3,1,4,5]), ...
                'linestyle','none', ...
                'linewidth', 3, ...
                'color', cgray)



set(gca, 'fontsize',f.gca)
ylabel('Fractional urine K^+ excretion (%)', 'fontsize',f.ylab)
ylim([0.0,100])
xlim([0.5,5.5])

xticks([1,2,3,4,5])
labs = {lab2, lab3,...
    lab1,lab4, lab5};
xticklabels(labs)
grid on