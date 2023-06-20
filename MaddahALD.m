% function for Maddah ALD equation
norm_Aldo = 0.49; % nmol/L
m_K_ALDO = 951.2;
m_Na_ALDO = 15.569;
norm_plasma_K = 0.0042;

pars = [norm_Aldo, m_K_ALDO, m_Na_ALDO, norm_plasma_K];
pars2 = pars;
pars2(2) = 100;

Kplas_vals = 3.5:0.01:5.0;
ALDO_vals = zeros(size(Kplas_vals));
Cal_vals = zeros(size(Kplas_vals));
ALDO_vals2 = zeros(size(Kplas_vals));
for ii = 1:length(Kplas_vals)
    Kplas = Kplas_vals(ii);
    ALDO_vals(ii) = get_ALDO(Kplas/1000, pars);
    Cal_vals(ii) = get_Cal(Kplas);
    ALDO_vals2(ii) = get_ALDO(Kplas/1000, pars2);
end

figure(25)
clf;
plot(Kplas_vals, ALDO_vals/norm_Aldo, 'linewidth', 2)
hold on
plot(Kplas_vals, Cal_vals/85, 'linewidth', 2)
plot(Kplas_vals, ALDO_vals2/norm_Aldo, 'linewidth', 2)
xlabel('[K^+]_p (mmol/L)')
ylabel('N_{al}')
legend('Maddah', 'Stadt', 'Maddah2')
grid on


function ALDO = get_ALDO(Kplas, pars)
    ALDO_0  = pars(1);
    mKALDO  = pars(2);
    %mNaALDO = pars(3); % left out b/c will be zero
    Kplas0  = pars(4);

    ALDO = ALDO_0 * exp(mKALDO * (Kplas - Kplas0));
end

function Cal = get_Cal(Kecf)
    N_al0 = 1;
    tspan = [0,5000];

    options = odeset('RelTol',1.0e-6,'AbsTol',1e-9);
    [t,y] = ode15s(@(t,y) ALD_eqns(t, y, Kecf),...
                    tspan, N_al0, options);

    Cal = 85*y(end);
end

function dydt = ALD_eqns(t,y,Kecf)
    
    K_ECFtot = Kecf;

    Kecf_total = 4.2;
    Csod =144;
    xi_par = 2;
    T_al = 60;

    N_al = y;

    xi_ksod = max(0,((K_ECFtot/Csod)/(Kecf_total/144/(xi_par+1))-xi_par));
    N_als = xi_ksod;
    dydt = (1/T_al)*(N_als - N_al);
end

