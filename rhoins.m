clear all;

%% get data
tvals = 0:10:600;
pars = set_params();
Cins_vals = zeros(size(tvals));
rhoins_vals = zeros(size(tvals));

rhoins_lin_vals = zeros(size(tvals));

max_rho = 1.11

for ii = 1:length(tvals)
    t = tvals(ii);
    [C, rho] = get_insvals_orig(t,pars);
    
    Cins_vals(ii) = C;
    rhoins_vals(ii) = rho;

    [~, rho] = get_insvals_lin(t,pars, max_rho);
    rhoins_lin_vals(ii) = rho;
end

%% make figures
lw = 2;
figure(11);
clf;
subplot(1,2,1)
plot(tvals, Cins_vals, 'linewidth', lw)
xlabel('t')
ylabel('C_{insulin}')

subplot(1,2,2)
hold on
plot(tvals, rhoins_vals, 'linewidth',lw)
plot(tvals, rhoins_lin_vals, 'linewidth', lw)
legend('orig', 'linear')
xlabel('t')
ylabel('\rho_{insulin}')

figure(12)
clf;
hold on
plot(Cins_vals, rhoins_vals, 'linewidth',lw)
plot(Cins_vals, rhoins_lin_vals, 'linewidth',lw)
legend('orig', 'linear')
xlabel('C_{insulin}')
ylabel('\rho_{insulin}')
%% functions
%------------------------------------------------------
function [C_insulin, rho_insulin] = get_insvals_orig(t, pars)
    meal_start = 0;
    t_insulin = t- meal_start;
    C_insulin = get_Cinsulin(t_insulin);
    L = 100; x0 = 0.5381; k = 1.069;
    ins_A = pars.A_insulin; ins_B = 100*pars.B_insulin;
    temp = (ins_A.*(L./(1+exp(-k.*(log10(C_insulin)-log10(x0)))))+ ins_B)./100;
    rho_insulin = max(1.0,temp);
end

function [C_insulin, rho_insulin] = get_insvals_lin(t,pars,max_rho)
    meal_start = 0;
    t_insulin = t- meal_start;
    C_insulin = get_Cinsulin(t_insulin);
    
    %alpha = 1.15;
    %m = (1.15 - 1.0)/(0.325 - get_Cinsulin(pars.t_insulin_ss));
    %b = 1.106 - 0.325 * m;
    m = (max_rho - 1.0)/(0.325 - get_Cinsulin(pars.t_insulin_ss));
    b = max_rho - 0.325 * m;

    rho_insulin = max(1.0, m*C_insulin + b);
end
