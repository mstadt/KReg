clear all
% script to write my own kidney equations
% work from a desired steady state
% IDEAS:
    % what should the DCT, CNT, CCD cellular concentrations be?
        % This can be found in the nephron model
        % NOTE: intracellular K shouldn't really change much... is this
        % necessary to have as an ODE?

Kplas = 4.2/1000; % mEq/mL, plasma K concentration
Naplas = 140/1000; % mEq/mL, plasma Na concentration

GFR    = 105; % ml/min/1.73 m^2 (normal range 90-120 ml/min)
Nnephs = 2e6; % each kidney has about 1e6 nephrons
SNGFR  = GFR/Nnephs; % single nephron GFR

fil_K  = Kplas  * SNGFR; % mEq/min/nephron, filtration in a single nephron
fil_Na = Naplas * SNGFR; % mEq/min/nephron, filtration in a single nephron

%% proximal segment (reabsorbs ~92% of filtered load)
frac_psNareab = 0.92; % Maddah has this dependent on Na intake, see if want to add later...
Phi_psNareab  = frac_psNareab * fil_Na;
Phi_psKreab   = frac_psNareab * fil_K; % ps K reabsorption depends on Na

LoH_Na_out    = fil_Na - Phi_psNareab; % Na flow at start of DT
LoH_K_out     = fil_K  - Phi_psKreab; % K flow at start of DT
LoH_vol_out   = SNGFR * (1-frac_psNareab); % water flow at start of DT

%% distal convoluted tubule 
% Maddah parameters
GFR = 105; % ml/min, code
N_nephs = 2e6; % count of nephrons, code/Table 1
C_sod = 140/1000; % mEq/mL, code/Table 1
intake_Na = 100/24/60; %0.07; % mEq/min, code
F = 97; % Faraday constant, C/mmol, code
R = 8.3145; % universal gas constant, J/mol.K, code
T = 310.6; % normal body temp, K, code
V_L = -18.4; % mV, luminal electrical potential
V_B = -78.2; % mV, basolateral electrical potential


v_L = (F*V_L)/(R*T); % normalized luminal membrane potential
v_B = (F*V_B)/(R*T); % normalized basolateral membrane potential
h_L = 2.4935e-5; % cm/min; luminal cell membrane permeability to K through ROMK channels
h_B = 3.43e-5; % cm/min; basolateral cell membrane permeability to K through Kir channels
J_NKAmax = 14.66e-5; % mmol/min/cm^2

% DCT
SV_DCT = 7.5e-3; % m^3/m^2, surface to volume ratio of DCT cell
DCT_diam = 0.0015; % cm
DCT_len = 0.5; % cm % NOTE: in Layton & Layton is 0.2 ---- NOT 0.5....
Fr_DCT = 1.0; % fraction of DCT luminal surfact area covered by secreting cells
DCT_volume = pi*((DCT_diam/2)^2)*DCT_len; % 
DCT_area = pi*DCT_diam*DCT_len*Fr_DCT; %
eta_vol_DCT = 0; % fractional volume reabsorption along DCT (almost none)
DCT_lumen_K_conc0 = 0.0055; % mEq/ml, calcNormParams.R
cell_Naconc = 8/1000; % mEq/mL


DCT_cell_Kcon0   = 0.15; DCT_cell_Kcon = DCT_cell_Kcon0;
DCT_lumen_Kamt0  = 4.8597e-09; DCT_lumen_Kamt = DCT_lumen_Kamt0;

DCT_lumen_Kcon = max(0, DCT_lumen_Kamt) / DCT_volume;
J_DCT_lumen = (h_L * v_L * (1/(1-exp(-v_L)))*...
                        (DCT_cell_Kcon * exp(-v_L) - DCT_lumen_Kcon)); % eq 7, DCT

J_DCT_pass_baso = (h_B*v_B*(1/(1-exp(-v_B)))*...
                                (Kplas- DCT_cell_Kcon * exp(-v_B))); % eq 10, DCT

DCT_K_Na = 0.2*(1+DCT_cell_Kcon/8.33);
K_K = 0.1*(1+(C_sod/18.5));

J_DCT_act_baso = ((-2/3)*J_NKAmax*...
                        (cell_Naconc/(cell_Naconc - DCT_K_Na))^3*...
                         (Kplas/(Kplas + K_K))^2); % eq 12, DCT

J_DCT_baso = J_DCT_pass_baso + J_DCT_act_baso; % Eq 9, DCT


ddt_DCT_cell_Kcon = ((1/SV_DCT)*(J_DCT_baso - J_DCT_lumen)) 

ddt_DCT_lumen_Kamt = (LoH_K_out + J_DCT_lumen*DCT_area ...
                               - DCT_lumen_Kcon * LoH_vol_out * (1-eta_vol_DCT))

% d/dt (DCT_cell_Kcon) = ((1/SV_DCT)*(J_DCT_baso - J_DCT_lumen)); % eq 15, DC

% This is zero when J_DCT_baso = J_DCT_lumen




% d/dt (DCT_lumen_Kamt) = (LoH_K_out + J_DCT_lumen*DCT_area ...
%                                - DCT_lumen_Kcon * LoH_vol_out * (1-eta_vol_DCT));

%% connecting tubule 

% DCT + CNT: 
%     normal or high K intake: secretes 10-50% of filtered load
%     low K: reabsorbs 3% of filtered load

%% cortical collecting duct
% CCD: 
%    normal/high K intake: secretes 5-30% of filtered K
%    low K: reabsorbs ~10% of filted load

%% medullary collecting duct



% NOTE: this doesn't seem like should be a time derivative, but should
% rather just be the general fluxes..... like J_DCT_lumen is just the flux
% from inside the cells into the DCT luminal area rather than this being
% dynamic.....
% --- let's just try to get the kidney working and then see what I think I
% should do....