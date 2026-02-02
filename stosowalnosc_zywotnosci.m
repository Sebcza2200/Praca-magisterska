% stosowalnosc.m - Analiza porównawcza metod obliczania żywotności
% Wykresy w jednym oknie (subplot 1x2):
% 1. Lewy: W funkcji wysokości perygeum (h_p: 100-2000 km, e=0.001)
% 2. Prawy: W funkcji mimośrodu (e: 10^-8 - 0.5, h_p=500 km)
clc; clear; close all;
fprintf('--- ROZPOCZYNANIE ANALIZY STOSOWALNOŚCI METOD ---\n');

% =========================================================================
% 1. INICJALIZACJA DANYCH I STAŁYCH
% =========================================================================
try
    data = readmatrix("MSISE.txt");
    Const.h0 = data(:,1); Const.H = data(:,2); Const.rho0 = data(:,3);
    load('coeffs_X.mat');
    compute_X_handle = @(e_query, h_query) compute_X(e_query, h_query, ...
        coeff1, deg1, mx1, sx1, my1, sy1, coeff2, deg2, mx2, sx2, my2, sy2);
catch
    error('Brakuje plików MSISE.txt lub coeffs_X.mat');
end

% Stałe fizyczne
Const.R_E = 6.378e6; Const.mi = 3.986e14; Const.J_2 = 1082.6e-6; Const.c = 2.9979e8;
Const.epsilon_deg = 23.44; Const.a_S = 1.496e11; Const.r_S = 1.496e11; Const.G = 6.674e-11;
Const.n_sun = 1.99e-7; Const.n_moon = 2.66e-6; Const.i_sundeg = 23.49; Const.i_moondeg = 28.6;
Const.hpp = 78000; Const.r_pp = Const.hpp;

% =========================================================================
% 2. PARAMETRY SATELITY
% =========================================================================
m = 575; A_D = 8.0; A_SRP = 10.0; C_D = 2.2; C_P = 1.2;
C_B = m / (C_D * A_D);
e_const = 0.001;
i_const = deg2rad(43); 
OMEGA_const = deg2rad(56);
omega_const = deg2rad(275);
t_date = datetime('now', 'TimeZone', 'Europe/Warsaw');

% #########################################################################
% CZĘŚĆ A: OBLICZENIA W FUNKCJI WYSOKOŚCI (100 - 2000 km)
% #########################################################################
fprintf('\n>>> CZĘŚĆ A: Analiza w funkcji wysokości h_p...\n');
h_vec = 100:50:2000; 
n_steps = length(h_vec);
Z1_res = nan(1, n_steps); Z2_res = nan(1, n_steps); 
Z3_res = nan(1, n_steps); Z4_res = nan(1, n_steps);

for k = 1:n_steps
    h_p_km = h_vec(k);
    h_p = h_p_km * 1000; 
    r_p = h_p + Const.R_E;
    a = r_p / (1 - e_const);
    n = sqrt(Const.mi / a^3);
    T_S = 2 * pi * sqrt(a^3 / Const.mi);
    
    [~, res_vec] = ObliczPerturbacje(a, e_const, i_const, OMEGA_const, omega_const, ...
                                     m, A_D, A_SRP, C_D, C_P, t_date, Const);
    
    [z1, z2, z3, z4] = ObliczZywotnosci(a, e_const, n, r_p, h_p, C_B, T_S, ...
                                        res_vec.dadt, res_vec.dedt, ...
                                        Const.r_pp + Const.R_E, ...
                                        Const.h0, Const.H, Const.rho0, compute_X_handle);
    Z1_res(k) = z1; Z2_res(k) = z2; Z3_res(k) = z3; Z4_res(k) = z4;
end

% #########################################################################
% CZĘŚĆ B: OBLICZENIA W FUNKCJI MIMOŚRODU (h_p = 500 km)
% #########################################################################
fprintf('\n>>> CZĘŚĆ B: Analiza w funkcji mimośrodu e...\n');
h_p_fixed_km = 500;
h_p_fixed = h_p_fixed_km * 1000;
r_p_fixed = h_p_fixed + Const.R_E;
e_vec = logspace(-8, log10(0.5), 100); 
n_steps_e = length(e_vec);
Z1_e = nan(1, n_steps_e); Z2_e = nan(1, n_steps_e); 
Z3_e = nan(1, n_steps_e); Z4_e = nan(1, n_steps_e);

for k = 1:n_steps_e
    e_curr = e_vec(k);
    a_curr = r_p_fixed / (1 - e_curr); 
    n_curr = sqrt(Const.mi / a_curr^3);
    T_S_curr = 2 * pi * sqrt(a_curr^3 / Const.mi);
    
    [~, res_vec] = ObliczPerturbacje(a_curr, e_curr, i_const, OMEGA_const, omega_const, ...
                                     m, A_D, A_SRP, C_D, C_P, t_date, Const);
    
    [z1, z2, z3, z4] = ObliczZywotnosci(a_curr, e_curr, n_curr, r_p_fixed, h_p_fixed, C_B, T_S_curr, ...
                                        res_vec.dadt, res_vec.dedt, ...
                                        Const.r_pp + Const.R_E, ...
                                        Const.h0, Const.H, Const.rho0, compute_X_handle);
    Z1_e(k) = z1; Z2_e(k) = z2; Z3_e(k) = z3; Z4_e(k) = z4;
end

% Sanityzacja dla loglog
Z1_e(Z1_e <= 0) = NaN; Z2_e(Z2_e <= 0) = NaN; 
Z3_e(Z3_e <= 0) = NaN; Z4_e(Z4_e <= 0) = NaN;


% #########################################################################
% CZĘŚĆ C: RYSOWANIE WYKRESÓW
% #########################################################################
figure('Name', 'Analiza Stosowalności Metod', 'Color', 'w', 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.6]);

% --- LEWY WYKRES: Wysokość ---
subplot(1, 2, 1);
semilogy(h_vec, Z1_res, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Metoda 1 (King-Hele)'); hold on;
semilogy(h_vec, Z2_res, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Metoda 2 (Approx)');
semilogy(h_vec, Z3_res, 'g-.', 'LineWidth', 1.5, 'DisplayName', 'Metoda 3 (Empiryczna)');
semilogy(h_vec, Z4_res, 'k:', 'LineWidth', 1.8, 'DisplayName', 'Metoda 4 (Własna)');

xlabel('Wysokość perygeum h_p [km]');
ylabel('Żywotność [lata]');
title({'Porównanie metod wyznaczania żywotności dla zmiennej wysokości perigeum'; ...
       ['(e = ' num2str(e_const) ', m = ' num2str(m) ' kg)']});
grid on;
xlim([100 2000]);

% Linie graniczne (Dodane do legendy poprzez DisplayName, usunięte napisy tekstowe)
lim_1 = 150; lim_2 = 550;
xline(lim_1, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5, 'DisplayName', 'Dolna granica M3');
xline(lim_2, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5, 'DisplayName', 'Górna granica M3');

legend('Location', 'best');


% --- PRAWY WYKRES: Mimośród ---
subplot(1, 2, 2);
loglog(e_vec, Z1_e, 'r-', 'LineWidth', 2.0, 'DisplayName', 'Metoda 1 (King-Hele)'); hold on;
loglog(e_vec, Z2_e, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Metoda 2 (Approx)');
loglog(e_vec, Z3_e, 'g-.', 'LineWidth', 1.5, 'DisplayName', 'Metoda 3 (Empiryczna)');
loglog(e_vec, Z4_e, 'k:', 'LineWidth', 1.5, 'DisplayName', 'Metoda 4 (Własna)');

xlabel('Mimośród e [-]');
ylabel('Żywotność [lata]');
title({'Porównanie metod wyznaczania żywotności dla zmiennej mimośrodowości'; ...
       ['(h_p = ' num2str(h_p_fixed_km) ' km, m = ' num2str(m) ' kg)']});
grid on;
xlim([1e-8 0.05]);

% Granica Metody 2 (Dodana do legendy poprzez DisplayName, usunięty napis tekstowy)
xline(0.005, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'DisplayName', 'Granica stosowalności M2');

legend('Location', 'best');

disp('Gotowe.');


% =========================================================================
% FUNKCJE LOKALNE
% =========================================================================
function [rho_MSIS, H_ref] = MSIS(r_query_meters, h0, H, rho0)
    R_E = 6378000;
    h_query = (r_query_meters - R_E) / 1000; 
    
    idx = find(h0 <= h_query, 1, 'last');
    if isempty(idx), idx = 1; end
    
    if h_query == h0(idx)
        rho_MSIS = rho0(idx); H_ref = H(idx); return;
    end
    
    H_i = H(idx); h0_i = h0(idx); rho0_i = rho0(idx);
    rho_MSIS = rho0_i * exp(-(h_query - h0_i)/H_i);
    H_ref = H_i;
end

function [Zywotnosc_1, Zywotnosc_2, Zywotnosc_3, Zywotnosc_4] = ObliczZywotnosci(a, e, n, r_p, h_p, C_B, T_S, dadt, dedt, rpp, h0, H, rho0, compute_X_handle)
    mi = 3.986004418e14; 
    R_E = 6378137;       
    
    % METODA 1
    [rhop0, Hp0_km] = MSIS(r_p, h0, H, rho0); 
    Hp0 = Hp0_km * 1000;
    B = n / C_B * rhop0 * a * e * besselj(1, a * e / Hp0) * exp(-e * (1 + a / Hp0));
    if B ~= 0
        Zywotnosc_1s = e ^ 2 * (1 - 11 / 6 * e + 29 / 16 * e ^ 2 + 7 / 8 * Hp0 / a) / 2 / B; 
    else
        Zywotnosc_1s = Inf;
    end
    Zywotnosc_1 = Zywotnosc_1s / 3600 / 24 / 365;
    
    % METODA 2
    if e > 0.5
        Zywotnosc_2 = NaN;
    else
        Zywotnosc_2s = Hp0 * C_B / rhop0 / sqrt(mi * a) * (1 - exp(-(a - R_E) / Hp0) * (1 + (a - R_E) / 2 / a));
        Zywotnosc_2 = Zywotnosc_2s / 3600 / 24 / 365;
    end
    
    % METODA 3
    try
        if ~isempty(compute_X_handle)
            NC = compute_X_handle(e, h_p / 1000);
            Zywotnosc_3s = NC * C_B * T_S; 
            Zywotnosc_3 = Zywotnosc_3s / 3600 / 24 / 365;
        else
            Zywotnosc_3 = NaN;
        end
    catch
        Zywotnosc_3 = NaN;
    end
    
    % METODA 4
    dadt_srednia = mean(dadt);
    dedt_srednia = mean(dedt);
    
    if dadt_srednia >= 0
         Zywotnosc_4 = NaN; % Nie spada
    else
        a_delta = dadt_srednia * dedt_srednia;
        b_delta = a * dedt_srednia + dadt_srednia * (1 - e);
        c_delta = r_p - rpp;
        Delta = b_delta ^ 2 - 4 * a_delta * c_delta;
        
        if Delta >= 0
            Zywotnosc_4x = (-b_delta - sqrt(Delta)) / 2 / a_delta;
            Zywotnosc_4s = (-b_delta + sqrt(Delta)) / 2 / a_delta;
            if Zywotnosc_4x > 0 && Zywotnosc_4s > 0
                 Zywotnosc_4_raw = min(Zywotnosc_4x, Zywotnosc_4s); 
            else
                 Zywotnosc_4_raw = max(Zywotnosc_4x, Zywotnosc_4s);
            end
            Zywotnosc_4 = Zywotnosc_4_raw / 365; 
        else
            Zywotnosc_4 = NaN;
        end
    end
end

function [avg, vec] = ObliczPerturbacje(a, e, i, OMEGA, omega, m, A_D, A_SRP, C_D, C_P, t_date, Const)
    mi = Const.mi; R_E = Const.R_E; J_2 = Const.J_2; c = Const.c;
    epsilon_deg = Const.epsilon_deg; a_S = Const.a_S; r_S = Const.r_S;
    n_sun = Const.n_sun; n_moon = Const.n_moon;
    h0 = Const.h0; H = Const.H; rho0 = Const.rho0;
    epsilon = deg2rad(epsilon_deg);
    n = sqrt(mi / a ^ 3);
    n_obr = n * 3600 * 24 / 2 / pi;
    
    if isempty(t_date.TimeZone), t_ref = datetime(year(t_date), 7, 4);
    else, t_ref = datetime(year(t_date), 7, 4, 'TimeZone', t_date.TimeZone); end
    
    D = 2 * pi / 365 * days(t_date - t_ref);
    day_of_year = day(t_date, 'dayofyear');
    lambda_deg = mod((day_of_year - 80) * 360 / 365, 360);
    lambda = deg2rad(lambda_deg);
    f_SRP = 1358 / (1.0004 + 0.0334 * cos(D)) / c * A_SRP / m * C_P * (a_S / r_S) ^ 2;
    dOMEGAniejed = (-3 / 2 * J_2 * n * cos(i) / (1 - e ^ 2) ^ 2 * (R_E / a) ^ 2) * 24 * 3600 * 180 / pi;
    domeganiejed = (-3 / 4 * n * J_2 * (1 - 5 * cos(i) ^ 2) / (1 - e ^ 2) ^ 2 * (R_E / a) ^ 2) * 24 * 3600 * 180 / pi;
    dMdtniejed = (n + 3 / 4 * n * J_2 * (3 * cos(i) ^ 2 - 1) / (1 - e ^ 2) ^ (3 / 2) * (R_E / a) ^ 2) * 24 * 3600 * 180 / pi;
    dOMEGAmoon = -3.4 * 10 ^ (-3) * cos(i) / n_obr;
    domegamoon = 1.7 * 10 ^ (-3) * (5 * cos(i) ^ 2 - 1) / n_obr;
    dOMEGAsun = -1.5 * 10 ^ (-3) * cos(i) / n_obr;
    domegasun = 0.8 * 10 ^ (-3) * (5 * cos(i) ^ 2 - 1) / n_obr;
    
    theta_step = 2.0; 
    theta_deg = 0:theta_step:359; 
    len_th = length(theta_deg);
    
    dadt = zeros(1, len_th); dadtatm = zeros(1, len_th); dadtSRP = zeros(1, len_th);
    dedt = zeros(1, len_th); dedtatm = zeros(1, len_th); dedtSRP = zeros(1, len_th);
    dMdt = zeros(1, len_th); dMdtSRP = zeros(1, len_th);
    
    for k = 1:len_th
        theta_val = deg2rad(theta_deg(k));
        u_val = omega + theta_val;
        r_val = a * (1 - e^2) / (1 + e * cos(theta_val));
        v_k = sqrt(mi * (2/r_val - 1/a));
        E = 2 * atan(sqrt((1 - e)/(1 + e)) * tan(theta_val / 2));
        [rho_val, ~] = MSIS(r_val, h0, H, rho0);
        F_D_val = 0.5 * rho_val * A_D * v_k ^ 2 * C_D;
        
        dadtatm(k) = (- 2 * F_D_val / m / n / sqrt(1 - e ^ 2) * sqrt(e ^ 2 + 1 + 2 * e * cos(theta_val))) * 3600 * 24;
        dedtatm(k) = (- F_D_val / m / n / a / sqrt(e ^ 2 + 1 + 2 * e * cos(theta_val)) * sqrt(1 - e ^ 2) * (2 * e + 2 * cos(theta_val))) * 3600 * 24;
        
        R_val = f_SRP * ( - cos(i/2)^2 * cos(epsilon/2)^2 * cos(lambda - (u_val) - OMEGA) - sin(i/2)^2 * sin(epsilon/2)^2 * cos(lambda - (u_val) + OMEGA) - 0.5 * sin(i) * sin(epsilon) * (cos(lambda - (u_val)) - cos(-lambda - (u_val))) - sin(i/2)^2 * cos(epsilon/2)^2 * cos(-lambda - (u_val) + OMEGA) - cos(i/2)^2 * sin(epsilon/2)^2 * cos(-lambda - (u_val) - OMEGA));
        T_val = f_SRP * ( - cos(i/2)^2 * cos(epsilon/2)^2 * sin(lambda - (u_val) - OMEGA) - sin(i/2)^2 * sin(epsilon/2)^2 * sin(lambda - (u_val) + OMEGA) - 0.5 * sin(i) * sin(epsilon) * (sin(lambda - (u_val)) - sin(-lambda - (u_val))) - sin(i/2)^2 * cos(epsilon/2)^2 * sin(-lambda - (u_val) + OMEGA) - cos(i/2)^2 * sin(epsilon/2)^2 * sin(-lambda - (u_val) - OMEGA));
        
        dadtSRP(k) = (2 / n / sqrt(1 - e ^ 2) * (e * sin(theta_val) * R_val + (1 + e * cos(theta_val))* T_val)) * 24 * 3600;
        dedtSRP(k) = (sqrt(1 - e ^ 2) / n / a * (sin(theta_val) * R_val + (cos(E) + cos(theta_val)) * T_val)) * 24 * 3600;
        dMdtSRP(k) = (n + (1 - e ^ 2) / n / a / e * ((-2 * e / (1 + e * cos(theta_val)) + cos(theta_val)) * R_val - (1 + 1 / (1 + e * cos(theta_val))) * T_val * sin(theta_val))) * 24 * 3600 * 180 / pi;
        
        dadt(k) = dadtatm(k) + dadtSRP(k);
        dedt(k) = dedtatm(k) + dedtSRP(k);
        dMdt(k) = dMdtSRP(k) + dMdtniejed;
    end
    
    avg.da_dt = mean(dadt);
    avg.de_dt = mean(dedt);
    vec.dadt = dadt;
    vec.dedt = dedt;
end
function X = compute_X(e_query, h_query, coeff1, deg1, mx1, sx1, my1, sy1, coeff2, deg2, mx2, sx2, my2, sy2)
    X = zeros(size(e_query));
    for k = 1:numel(e_query)
        e = e_query(k); h = h_query(k);
        if e <= 0.1
            xn = (h - mx1)/sx1; yn = (e - my1)/sy1;
            A = buildPolyXY([xn],[yn], deg1);
            log10N = A * coeff1(:);
        else
            xn = (h - mx2)/sx2; yn = (e - my2)/sy2;
            A = buildPolyXY([xn],[yn], deg2);
            log10N = A * coeff2(:);
        end
        X(k) = 10.^log10N;
    end
end
function A = buildPolyXY(x, y, deg)
    count = 0;
    for d = 0:deg
        for i = 0:d
            j = d - i; count = count + 1;
            A(:,count) = (x(:).^i) .* (y(:).^j);
        end
    end
end