% Test equations for kidney at baseline
% [K+]_p before putting in the function
% Notes:
%  - likely need to convert some parameters to get the urinary excretion
%      matching what I need it to
%  - I have this partially coded in K-QSP repo on my github


% Maddah parameters
GFR = 105; % ml/min, code
N_nephs = 2e6; % count of nephrons, code/Table 1
C_sod = 140/1000; % mEq/mL, code/Table 1
intake_Na = 100/24/60; %0.07; % mEq/min, code
Kplas = 4.0/1000; % mEq/mL -- baseline plasma [K+] in mmol/mL
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
SV_DCT = 0.75e-3; % m^3/m^2, surface to volume ratio of DCT cell
DCT_diam = 0.0015; % cm
DCT_len = 0.5; % cm
Fr_DCT = 1.0; % fraction of DCT luminal surfact area covered by secreting cells
DCT_volume = pi*((DCT_diam/2)^2)*DCT_len; % 
DCT_area = pi*DCT_diam*DCT_len*Fr_DCT; %
eta_vol_DCT = 0; % fractional volume reabsorption along DCT (almost none)
DCT_lumen_K_conc0 = 0.0055; % mEq/ml, calcNormParams.R
cell_Naconc = 8/1000; % mEq/mL

% CNT
CNT_diam      = 2*12*10^-4; % cm
CNT_len       = 0.4; % cm
Fr_CNT        = 0.6; % fraction of CNT luminal surface area covered by secreting cells
CNT_volume  = pi*((CNT_diam/2)^2)*CNT_len;
CNT_lumen_K_conc0 = 0.0072; % mEq/ml, calcNormParams.R

% CCD
CCD_diam      = 0.0025; % cm
CCD_len       = 0.2; % cm
Fr_CCD        = 0.6; % fraction of CCD luminal surface area covered by secreting cells
CCD_volume  = pi*((CCD_diam/2)^2)*CCD_len;
CCD_lumen_K_conc0 = 0.032; % mEq/ml, calcNormParams.R


%% renal handling
% filtration
filK = (GFR/N_nephs)*Kplas; % K filtration, Eq 1
filNa = (GFR/N_nephs)*C_sod; % Na filtration, Eq 2

%% proximal segments
LoH_Na_out = max(intake_Na/((1-0.95)*N_nephs), 0.96*filNa); % Eq 3, Note: in code maximized reabsorption to 0.96
gamma_psNareab = (filNa - LoH_Na_out)/(filNa); % fractional Na reabsorption in PT/LoH, Eq 4

LoH_K_out = filK*(1-gamma_psNareab); % K flow into DCT, Eq 5
LoH_vol_out = (GFR/N_nephs)*(1-gamma_psNareab); % volume flow into DCT, Eq. 6

%% DCT
% NOTE: okay....so this is actually done through ODES... so may need to
% pick my steady state values to be able to figure this out!
% variables: Lumen K amount, Intracellular [K+]
% % fluxes

% steady state variable values
DCT_cell_Kcon_ss = 0.15; % model code, (runToEquilibrium.R)
DCT_lumen_Kamt_ss = DCT_lumen_K_conc0 * DCT_volume;
% ODE equations
% d/dt (DCT_cell_Kcon) = ((1/SV_DCT)*(J_DCT_baso - J_DCT_lumen)); % eq 15, DCT
% d/dt (DCT_lumen_Kamt) = (LoH_K_out + J_DCT_lumen*DCT_area ...
%                                - DCT_lumen_Kcon * LoH_vol_out * (1-eta_vol_DCT)); % eq 16, DCT


DCT_lumen_Kcon = max(0, DCT_lumen_Kamt_ss) / DCT_volume;
J_DCT_lumen = (h_L * v_L * (1/(1-exp(-v_L)))*...
                        (DCT_cell_Kcon_ss * exp(-v_L) - DCT_lumen_Kcon)); % eq 7, DCT
J_DCT_pass_baso = (h_B*v_B*(1/(1-exp(-v_B)))*...
                                (Kplas- DCT_cell_Kcon_ss * exp(-v_B))); % eq 10, DCT

DCT_K_Na = 0.2*(1+DCT_cell_Kcon_ss/8.33);
K_K = 0.1*(1+(C_sod/18.5));

J_DCT_act_baso = ((-2/3)*J_NKAmax*...
                        (cell_Naconc/(cell_Naconc - DCT_K_Na))^3*...
                         (Kplas/(Kplas + K_K))^2); % eq 12, DCT

J_DCT_baso = J_DCT_pass_baso + J_DCT_act_baso; % Eq 9, DCT


%% CNT
CNT_cell_Kcon_ss = 0.15; % model code (runToEqulibrium.R)
CNT_lumen_Kamt_ss = CNT_lumen_K_conc0 * CNT_volume;

%% CCD
CCD_cell_Kcon_ss = 0.15; % model code (runToEqulibrium.R)
CCD_lumen_Kamt_ss = CCD_lumen_K_conc0 * CCD_volume;


