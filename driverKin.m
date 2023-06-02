clear all;
%% initialize parameter values
fprintf('loading params \n')
pars = set_params();
[params, parnames] = pars2vector(pars,0);

%% set initial conditions
temp = load('./SS/SS1.mat');
[IC, ~, ~] = getSS(temp.SS, params, 'do_figs', 0); % start at SS

%% Fasting state
% ODE options
t0 = 0;
tf = 6*60; % 6 hours of fasting
tspan = [t0, tf];
options = odeset('RelTol',1.0e-6,'AbsTol',1e-9); % ode solver settings

% simulation settings
MKX = 0; MKXslope = 0;
do_insulin = 0; % set to 1 when doing meal
do_FF      = 1; % set to 1 unless no FF effect
alt_sim    = 0; % only if have other options

fprintf('solving ODEs \n')
[t1,y1] = ode15s(@(t,y) kreg_eqns(t,y,params,...
                            'SS', false, ...
                            'alt_sim', alt_sim,...
                            'do_MKX', [MKX, MKXslope],...
                            'do_insulin', do_insulin,...
                            'do_FF', do_FF), ...
                            tspan, IC, options);