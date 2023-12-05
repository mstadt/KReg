clear all

% This script can be used to postprocess the local sensitivity analysis
% results

fname = "2023-11-30_LocalAnalysis_50days.csv";
dat = readtable(fname);
fname = "2023-11-30_LocalAnalysis_50days_base_vals.csv";
base_vals = readtable(fname);
base_Kplas = base_vals.K_plasma;
base_Kmusc = base_vals.K_musc;

%% Get nice parameter names
ParNames = [];
for ii = 1:length(dat.ParameterName)
    pname = dat.ParameterName{ii};
    %disp(pname)
    temp = change_parname(pname);
    ParNames = [ParNames; temp];
end

%% Get percent change
up_per = zeros(length(dat.ParameterName), 2);
down_per = zeros(size(up_per));
up_down_per = zeros(size(up_per));

% increase by 10%
up_per(:,1) = (dat.Kplas_Increase - base_Kplas)./base_Kplas * 100.0;
up_per(:,2) = (dat.Kmusc_Increase - base_Kmusc)./base_Kmusc * 100.0;
% decrease by 10%
down_per(:,1) = (dat.Kplas_Decrease - base_Kplas)./base_Kplas * 100.0;
down_per(:,2) = (dat.Kmusc_Decrease - base_Kmusc)./base_Kmusc * 100.0;
% increase by 10% - decrease by 10%
up_down_per(:,1) = (dat.Kplas_Increase - dat.Kplas_Decrease)./base_Kplas * 100.0;
up_down_per(:,2) = (dat.Kmusc_Increase - dat.Kmusc_Decrease)./base_Kmusc * 100.0;

%% Make figures
xlabels = ParNames; %dat.ParameterName;
ylabels = {'K_{plasma}', 'K_{intracellular}'};

% Up percent
% Find values less than 2%
up_per_h = round(up_per, 2,"significant"); % up_per for heatmap
[r,c] = find(abs(up_per_h) < 2.0);
for ii = 1:length(r)
    up_per_h(r(ii), c(ii)) = NaN;
end
figure(1)
clf;
h = heatmap(xlabels, ylabels, up_per_h',...
                'colormap', parula, ...
                'MissingDataColor', 'w', 'MissingDataLabel', '<2%');
h.FontSize = 14;

% Down percent
down_per_h = round(down_per, 2, "significant");
[r,c] = find(abs(down_per_h) < 2.0);
for ii = 1:length(r)
    down_per_h(r(ii),c(ii)) = NaN;
end
figure(2)
clf;
h = heatmap(xlabels, ylabels, down_per_h',...
                'colormap', parula, ...
                'MissingDataColor', 'w', 'MissingDataLabel', '<2%');
h.FontSize = 14;

% Up/Down percent
up_down_per_h = round(up_down_per, 2, "significant");
[r,c] = find(abs(up_down_per_h) < 2.01);
for ii = 1:length(r)
    up_down_per_h(r(ii), c(ii)) = NaN;
end
figure(3)
clf;
h = heatmap(xlabels, ylabels, up_down_per_h',...
                'colormap', parula, ...
                'MissingDataColor', 'w', 'MissingDataLabel', '<2%');
h.FontSize = 14;

fprintf('done \n')

%%
% Remove where both NaN
rows_notAllNan = find(~all(isnan(up_down_per_h), 2));

xlabels = ParNames(rows_notAllNan);
up_down_per_h2 = up_down_per_h(rows_notAllNan, :);
up_per_h2 = up_per_h(rows_notAllNan, :);
down_per_h2 = down_per_h(rows_notAllNan, :);

% Switch the order by largest impact on K_muscle
[sortvals, sortIDs] = sort(up_down_per_h2(:,2)); 

xlabels = xlabels(sortIDs);
up_down_per_h2 = up_down_per_h2(sortIDs,:);
up_per_h2 = up_per_h2(sortIDs,:);
down_per_h2 = down_per_h2(sortIDs,:);

figure(4)
clf
h = heatmap(xlabels, ylabels, up_per_h2',...
                'colormap', parula, ...
                'MissingDataColor', 'w', 'MissingDataLabel', '<2%');
h.FontSize = 14;

figure(5)
clf
h = heatmap(xlabels, ylabels, down_per_h2',...
                'colormap', parula, ...
                'MissingDataColor', 'w', 'MissingDataLabel', '<2%');
h.FontSize = 14;

figure(6)
clf;
h = heatmap(xlabels, ylabels, up_down_per_h2',...
                'colormap', parula, ...
                'MissingDataColor', 'w', 'MissingDataLabel', '<2%');
h.FontSize = 14;

fprintf('done \n')

% Up_per with rows_notAllNan
