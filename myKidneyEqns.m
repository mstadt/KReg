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

GFR = 105; % ml/min/1.73 m^2 (normal range 90-120 ml/min)
Nnephs = 2e6; % each kidney has about 1e6 nephrons
SNGFR = GFR/Nnephs; % single nephron GFR

fil_K = Kplas * SNGFR; % mEq/min/nephron, filtration in a single nephron
fil_Na = Naplas * SNGFR; % mEq/min/nephron, filtration in a single nephron

%% proximal segment
frac_psNareab = 0.92; % Maddah has this dependent on Na intake, see if want to add later...
Phi_psNareab = frac_psNareab * fil_Na;
Phi_psKreab  = frac_psNareab * fil_K; % ps K reabsorption depends on Na

LoH_Na_out = fil_Na - Phi_psNareab; % Na flow at start of DT
LoH_K_out = fil_K - Phi_psKreab; % K flow at start of DT
LoH_vol_out = SNGFR * (1-frac_psNareab); % water flow at start of DT

% d/dt (DCT_cell_Kcon) = ((1/SV_DCT)*(J_DCT_baso - J_DCT_lumen)); % eq 15, DC
% d/dt (DCT_lumen_Kamt) = (LoH_K_out + J_DCT_lumen*DCT_area ...
%                                - DCT_lumen_Kcon * LoH_vol_out * (1-eta_vol_DCT));

% NOTE: this doesn't seem like should be a time derivative, but should
% rather just be the general fluxes..... like J_DCT_lumen is just the flux
% from inside the cells into the DCT luminal area rather than this being
% dynamic.....
% --- let's just try to get the kidney working and then see what I think I
% should do....