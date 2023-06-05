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
eta_vol_CNT = 0.7; % fractional vol reabsorption

% CCD
CCD_diam      = 0.0025; % cm
CCD_len       = 0.2; % cm
Fr_CCD        = 0.6; % fraction of CCD luminal surface area covered by secreting cells
CCD_volume  = pi*((CCD_diam/2)^2)*CCD_len;
CCD_lumen_K_conc0 = 0.032; % mEq/ml, calcNormParams.R
eta_vol_CCD = 0.75; % fractional volume reabsorption in CCD

% MCD
m_plasK_MCD = 8.83e-7; % unitless, NOTE: this was in Maddah code, not table
K_plas0     = 4.2/1000;   % mEq/mL, normal plasma K concentration
Kreab_MCD0  = 7.293333e-09; % mEq/min
Nareab_MCD0 = 0.625/N_nephs; %3.125e-7; % mEq/min
m_NaK_MCD   = 0.69775; %min/mEq
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

DCT_vol_out = LoH_vol_out*(1-eta_vol_DCT);
DCT_K_out   = DCT_lumen_Kcon * DCT_vol_out;


%% CNT
CNT_cell_Kcon_ss = 0.15; % model code (runToEqulibrium.R)
CNT_lumen_Kamt_ss = CNT_lumen_K_conc0 * CNT_volume;

CNT_lumen_Kcon = CNT_lumen_Kamt_ss/CNT_volume;

J_CNT_lumen = (h_L * v_L * (1/(1-exp(-v_L)))*...
                           (CNT_cell_Kcon_ss*exp(-v_L) - CNT_lumen_Kcon)); % eq 7, CNT
J_CNT_pass_baso = (h_B*v_B*(1/(1-exp(-v_B)))*...
                                (Kplas - CNT_cell_Kcon_ss * exp(-v_B))); % eq 10, CNT

CNT_K_Na = (0.2*(1+CNT_cell_Kcon_ss/8.33)); % eq 13, CNT
 
J_CNT_act_baso = ((-2/3)*J_NKAmax*...
                        (cell_Naconc/(cell_Naconc - CNT_K_Na))^3*...
                         (Kplas/(Kplas + K_K))^2); % eq 12, CNT

J_CNT_baso = (J_CNT_pass_baso + J_CNT_act_baso); % eq 9, CNT

CNT_vol_out = DCT_vol_out*(1-eta_vol_CNT);
CNT_K_out = CNT_lumen_Kcon * CNT_vol_out;

%% CCD
CCD_cell_Kcon_ss = 0.15; % model code (runToEqulibrium.R)
CCD_lumen_Kamt_ss = CCD_lumen_K_conc0 * CCD_volume;

CCD_lumen_Kcon = CCD_lumen_Kamt_ss/CCD_volume;

J_CCD_lumen = (h_L * v_L * (1/(1-exp(-v_L)))*...
                           (CCD_cell_Kcon_ss*exp(-v_L) - CCD_lumen_Kcon)); % eq 7, CCD

CCD_K_Na = (0.2*(1+CCD_cell_Kcon_ss/8.33)); % eq 13, CCD

J_CCD_pass_baso = (h_B*v_B*(1/(1-exp(-v_B)))*...
                                (Kplas - CCD_cell_Kcon_ss * exp(-v_B))); % eq 10, CCD

J_CCD_act_baso = ((-2/3)*J_NKAmax*...
                        (cell_Naconc/(cell_Naconc - CCD_K_Na))^3*...
                         (Kplas/(Kplas + K_K))^2); % eq 12, CCD

CCD_vol_out = (CNT_vol_out * (1-eta_vol_CCD));

CCD_K_out = (CCD_lumen_Kcon * CCD_vol_out);

Phi_Na_CCD_out = (0.5 * LoH_Na_out); % line 173 Maddah model code

%% MCD (medullary collecting duct)
K_MCD_effect = (max(0, exp(m_plasK_MCD * (K_plas0 - Kplas)/K_plas0) - 1));
Phi_MCD_Nareab = max(0, Phi_Na_CCD_out - intake_Na/N_nephs);% eq 20, line 173
Phi_MCD_Ksec = max(0, m_NaK_MCD * (Nareab_MCD0 - Phi_MCD_Nareab)); % eq 19, NOTE: line 178 Maddah code has reversed difference so using here
Phi_MCD_Kreab_rate = (Kreab_MCD0 + K_MCD_effect + Phi_MCD_Ksec); % line 181 Maddah model code
K_transport_MCD = min(Phi_MCD_Kreab_rate, CCD_K_out);
Phi_K_out = (N_nephs*(CCD_K_out - K_transport_MCD)); % eq 21