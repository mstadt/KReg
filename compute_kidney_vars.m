function v = compute_kidney_vars(yvals, params, varargin)
% computes algebraic equations given output of ODE solver
% Input: yvals -- output from ODE solver, params -- parameter vector

% set variable names
MKgut_vals    = yvals(:,1); % amount of K in gut
MKplas_vals   = yvals(:,2); % amount of K in plasma
MKinter_vals  = yvals(:,3); % amount of K in interstitial space
MKmuscle_vals = yvals(:,4); % amount of K in muscle

% set parameter names
Phi_Kin_ss = params(1);
t_insulin_ss = params(2);
fecal_excretion = params(3);
kgut = params(4);
MKgutSS = params(5);
V_plasma = params(6);
V_interstitial = params(7);
V_muscle = params(8);
Kecf_total = params(9);
P_ECF = params(10);
Kmuscle_baseline = params(11);
Vmax = params(12);
Km = params(13);
P_muscle = params(14);
GFR_base = params(15);
eta_ptKreab_base = params(16);
eta_LoHKreab = params(17);
dtKsec_eq = params(18);
A_dtKsec = params(19);
B_dtKsec = params(20);
cdKsec_eq = params(21);
A_cdKsec = params(22);
B_cdKsec = params(23);
alpha_TGF = params(24);
A_cdKreab = params(25);
ALD_eq = params(26);
m_K_ALDO = params(27);
FF = params(28);
A_insulin = params(29);
B_insulin = params(30);

%% Get variable inputs
% default settings, varargin is used to change settings
SS = false; % compute SS solution
alt_sim = false; % use alternate equations
do_insulin = true;
do_FF = true;
MKX = 0;
Kintake = 0;
meal_start = 0;
highK_eff = 0;
TGF_eff = 0;
for i = 1:2:length(varargin)
    temp = varargin{i+1};
    if strcmp(varargin{i}, 'SS')
        SS = temp;
    elseif strcmp(varargin{i}, 'alt_sim')
        alt_sim = temp;
    elseif strcmp(varargin{i}, 'do_MKX')
        MKX = temp(1);
        MKslope = temp(2);
    elseif strcmp(varargin{i}, 'do_insulin')
        do_insulin = temp(1);
    elseif strcmp(varargin{i}, 'do_FF')
        do_FF = temp(1);
    elseif strcmp(varargin{i}, 'Kintake')
        Kintake = temp(1);
    elseif strcmp(varargin{i}, 'meal_time')
        meal_start = temp(1);
%     elseif strcmp(varargin{i}, 'highK_eff')
%         highK_eff = temp(1);
    elseif strcmp(varargin{i}, 'TGF_eff')
        TGF_eff = temp(1);
        alpha_TGF = temp(2);
        eta_ptKreab = temp(3);
%         if TGF_eff
%             fprintf('doing TGF_eff \n')
%         end
    else
        disp('WRONG VARARGIN INPUT')
        fprintf('What is this varargin input? %s \n', varargin{i})
        error('wrong varargin input')
    end % if
end %for



% concentrations
v.K_plas = MKplas_vals./V_plasma;
v.K_inter = MKinter_vals./V_interstitial;
v.K_muscle = MKmuscle_vals./V_muscle;
v.K_ECFtot = (MKplas_vals + MKinter_vals)./(V_plasma + V_interstitial);


% ALD impact
v.Nal_vals = exp(m_K_ALDO .* (v.K_ECFtot - Kecf_total));
v.C_al = v.Nal_vals .* ALD_eq;
v.gamma_al = A_dtKsec .* v.C_al .^B_dtKsec;
v.lambda_al = A_cdKsec .* v.C_al.^B_cdKsec;

if do_FF
    v.gamma_Kin = max(1, FF.*(MKgut_vals - MKgutSS));
else
    v.gamma_Kin = ones(size(MKgut_vals));
end

if MKX > 0
    v.omegaKic = max(0, (MKslope.*(v.K_muscle - Kmuscle_baseline) + 1));
else
    v.omegaKic = ones(size(v.K_muscle));
end

% renal K handling

eta_psKreab_base = eta_ptKreab_base + eta_LoHKreab;
if TGF_eff == 1
    eta_psKreab = eta_ptKreab + eta_LoHKreab;
    v.GFR = GFR_base + alpha_TGF * (eta_psKreab - eta_psKreab_base);
elseif TGF_eff == 2 % GFR only
    eta_ptKreab = eta_ptKreab_base; % PT K reab is baseline value
    eta_psKreab = eta_ptKreab + eta_LoHKreab;
    v.GFR = GFR_base + alpha_TGF * (eta_psKreab - eta_psKreab_base);
elseif TGF_eff == 3 % PT only
    eta_psKreab = eta_ptKreab + eta_LoHKreab;
    v.GFR = GFR_base;
else
    eta_ptKreab = eta_ptKreab_base; % PT K reab is baseline value
    eta_psKreab = eta_ptKreab + eta_LoHKreab;
    v.GFR = GFR_base; 
end

v.filK = v.GFR .* v.K_plas;

v.psKreab = eta_psKreab * v.filK;


% distal tubule
if MKX == 1
    v.eta_dtKsec = v.gamma_al .* v.gamma_Kin .* v.omegaKic;
else
    v.eta_dtKsec = v.gamma_al .* v.gamma_Kin;
end
v.dtKsec = dtKsec_eq .* v.eta_dtKsec;

% collecting duct
if MKX == 2
    v.eta_cdKsec = v.lambda_al .* v.omegaKic;
else
    v.eta_cdKsec = v.lambda_al;
end
v.cdKsec = cdKsec_eq .* v.eta_cdKsec;

if MKX == 3
    v.eta_cdKreab = v.omegaKic;
else
    v.eta_cdKreab = ones(size(v.omegaKic));
end

dtK = v.filK - v.psKreab + v.dtKsec;
v.cdKreab = dtK .* A_cdKreab.*v.eta_cdKreab;

v.UrineK = dtK + v.cdKsec - v.cdKreab;


end % compute_vars