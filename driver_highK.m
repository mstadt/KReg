% High K simulation experiment
% Options: PT impact of high K intake
clear all;

%------------- 
% Settings
%-------------
MealInsulin = 1; % insulin + KCl meals
len_meal = 30; % length of meal in minutes
doFF = 1; % do FF effect on DT
Kamt_high = 5 * 70 / 3; % high K intake, per meal
Kamt_control = 70 / 3; % control K intake, per meal
Kamt_meal = Kamt_control; % control sim

PT_highK = 0; % PT effect of high K intake
GFR_highK = 0; % GFR effect of high K intake
%------------------
%------------------

%% initialize parameter values
fprintf('loading params \n')
pars = set_params();
[params, parnames] = pars2vector(pars, 0);

%% set initial conditions
temp = load('./SS/SS_4vars.mat');
SS = temp.SS;
[IC, ~, ~] = getSS(SS, params, 'do_figs', 0); % start at SS

%% Intake for one day 
% TO DO: make in loop for multiple days
% sim settings
opts.do_insulin = MealInsulin;
opts.do_FF = doFF; 



MealTimes = [6, 12, 18]*60;
tspan = [0, MealTimes(1)];
fprintf('start sim \n')
[t1, y1, vals1] = fast_sim(IC, tspan, -6*60, params, opts);


% Meal 1
fprintf('meal 1 \n')
[t2, y2, vals2] = meal_sim(y1(end, :), t1(end), len_meal, Kamt_meal, ...
                            params, opts);

% Fasting post meal 1
[t3,y3,vals3] = fast_sim(y2(end,:), [t2(end), MealTimes(2)], t2(1), ...
                            params, opts);

fprintf('meal 2 \n')
% Meal 2
[t4, y4, vals4] = meal_sim(IC, t3(end), len_meal, Kamt_meal, ...
                            params, opts);

% Fasting post meal 2
[t5, y5, vals5] = fast_sim(y4(end,:), [t4(end), MealTimes(3)], t4(1), ...
                            params, opts);

% Meal 3
[t6, y6, vals6] = meal_sim(IC,t5(end),len_meal,Kamt_meal,...
                            params,opts);

fprintf('meal 3 \n')
% Fasting post meal 3
[t7, y7, vals7] = fast_sim(y6(end,:), [t6(end), 24*60+MealTimes(1)], t6(1),...
                            params, opts);

T = [t1;t2;t3;t4;t5;t6;t7]./60; Y = [y1;y2;y3;y4;y5;y6;y7];
%---------------
% plot results
%---------------
figure(1);
clf;
nr = 2; nc = 2;
lw = 3; lwgray = 2; lsgray = '--';
cmap = parula(6);
c1 = cmap(1,:); c2 = cmap(2,:); c3 = cmap(3,:); c4 = cmap(4,:);
cgraymap = gray(5);
cgray = cgraymap(3,:);
subplot(nr,nc,1)
plot(T, Y(:,1), 'linewidth', lw, color = c1)
xlabel('Time (hrs)')
ylabel('Gut amount')
title('Gut amount')
grid on

subplot(nr,nc,2)
hold on
plot(T,Y(:,2)/pars.V_plasma, 'linewidth',lw,'color',c2)
yline(3.5,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
yline(5.0,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
xlabel('Time (hrs)')
ylabel('Plasma [K^+]')
title('Plasma [K^+]')
grid on

subplot(nr,nc,3)
hold on
plot(T,Y(:,3)/pars.V_interstitial,'linewidth',lw,'color',c3)
yline(3.5,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
yline(5.0,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
xlabel('Time (hrs)')
ylabel('Interstitial [K^+]')
title('Interstitial [K^+]')
grid on

subplot(nr,nc,4)
hold on
plot(T,Y(:,4)/pars.V_muscle,'linewidth',lw,'color',c4)
yline(120,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
yline(140,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
xlabel('Time (hrs)')
ylabel('Intracellular [K^+]')
title('Intracellular [K^+]')
grid on

%------------
% functions
%------------
% Meal simulation
function [t, y, vals] = meal_sim(IC, t0, len_meal, Kamt, params, opts)
    tf = t0 + len_meal; % meal length
    tspan = [t0, tf];
    Kintake = Kamt / (tf - t0);
    options = odeset('RelTol', 1.0e-6, 'AbsTol', 1e-9); % ode solver settings
    [t, y] = ode15s(@(t,y) kreg_eqns(t,y,params,...
                                'do_insulin', opts.do_insulin,...
                                'do_FF', opts.do_FF,...
                                'meal_time', t0,... % t0 is start of meal
                                'Kintake', Kintake),...
                                tspan, IC, options);
    vals = compute_vars(t,y,params,...
                            'do_insulin', opts.do_insulin,...
                            'do_FF', opts.do_FF, ...
                            'meal_time', t0);
end

% Fasting simulation
function [t, y, vals] = fast_sim(IC, tspan, last_meal, params, opts)
    options = odeset('RelTol', 1.0e-6, 'AbsTol', 1e-9); % ode solver settings
    [t, y] = ode15s(@(t,y) kreg_eqns(t,y,params,...
                                'do_insulin', opts.do_insulin,...
                                'do_FF', opts.do_FF,...
                                'meal_time', last_meal,...
                                'Kintake', 0),... % fasting state
                                tspan, IC, options);
    vals = compute_vars(t,y,params,...
                            'do_insulin', opts.do_insulin,...
                            'do_FF', opts.do_FF,...
                            'meal_time', last_meal,...
                            'Kintake', 0); 
end