clc; clear; close all;
fprintf('=== ANALIZA ŻYWOTNOŚCI KONSTELACJI ===\n');
fprintf('Wybierz konstelację do analizy:\n');
fprintf('1. STARLINK\n'); fprintf('2. ONEWEB\n'); fprintf('3. KUIPER\n');
fprintf('4. GUOWANG\n'); fprintf('5. QIANFAN\n');
choice_num = input('Twój wybór (1-5): ');
switch choice_num
    case 1, target_constellation = 'STARLINK';
    case 2, target_constellation = 'ONEWEB';
    case 3, target_constellation = 'KUIPER';
    case 4, target_constellation = 'GUOWANG';
    case 5, target_constellation = 'QIANFAN';
    otherwise, error('Nieprawidłowy wybór (wpisz liczbę 1-5).');
end
fprintf('Inicjalizacja modelu atmosfery...\n');
data = readmatrix("MSISE.txt");
Const.h0 = data(:,1); Const.H = data(:,2); Const.rho0 = data(:,3); Const.rho01 = data(:,4); 
Const.rhomin = data(:,5); Const.rhomid = data(:,6); Const.rhomax = data(:,7);
Const.R_E = 6.378e6; Const.mi = 3.986e14; Const.J_2 = 1082.6e-6; Const.c = 2.9979e8; 
Const.epsilon_deg = 23.44; Const.a_S = 1.496e11; Const.r_S = 1.496e11; Const.hpp = 78000; 
Const.r_pp = Const.hpp + Const.R_E; C_D = 2.2; C_P = 1.2;

fprintf('Otwieranie pliku Excel i szukanie: %s...\n', target_constellation);
satellites = WczytajDane(target_constellation, Const);
num_sats = length(satellites);
if num_sats == 0
    error('Nie znaleziono satelitów dla tej konstelacji.');
end
fprintf('Załadowano %d satelitów. Rozpoczynam obliczenia.\n', num_sats);
res_hp = nan(num_sats, 1); res_Z1 = nan(num_sats, 1);
res_Z2 = nan(num_sats, 1); res_Z5 = nan(num_sats, 1);
fprintf('Rozpoczynanie symulacji... (może to zająć kilka minut)\n');

for i = 1:num_sats
    sat = satellites(i);
    m = sat.m; A_D = sat.A_D; A_SRP = sat.A_SRP;
    
    if isnan(m) || m <= 0 || isnan(A_D), continue; end

    a = sat.a; ecc = sat.e; inc = sat.i; raan = sat.OMEGA; omega = sat.omega;
    t_sim_start = sat.t_date;
    
    n = sqrt(Const.mi / a^3); r_p = a * (1 - ecc); hp_km = (r_p - Const.R_E) / 1000;
    
    if hp_km < 120 || hp_km > 2000, continue; end
    
    C_B = m / A_D / C_D; 
    
    [z1_val, z2_val] = ObliczZywotnosci(a, ecc, n, r_p, C_B, Const);
    z5_val = ObliczZ5(a, ecc, inc, raan, omega, m, ...
        A_D, A_SRP, C_D, C_P, t_sim_start, Const.rho0, Const);

    res_hp(i) = hp_km; res_Z1(i) = z1_val; res_Z2(i) = z2_val; res_Z5(i) = z5_val; e_val(i) = ecc;
    
    if mod(i, ceil(num_sats/10)) == 0 || i == num_sats
        fprintf('Postęp: %d / %d (%.1f%%) | ID: %d | h_p: %.1f km | Z5: %.2f lat\n', ...
            i, num_sats, (i/num_sats)*100, sat.ID, hp_km, z5_val);
    end
end

nazwa_pliku_txt = [target_constellation, '.txt']; %Zapis do pliku tekstowego
fid = fopen(nazwa_pliku_txt, 'w');
fprintf(fid, 'ID_Satelity Nazwa Satelity hp Z1 Z2 Z5 e\n');
for i = 1:num_sats
    fprintf(fid, '%d %s %.4f %.4f %.4f %.4f %8f\n', satellites(i).ID, satellites(i).Nazwa, ...
        res_hp(i), res_Z1(i), res_Z2(i), res_Z5(i), e_val(i));
end
fclose(fid);
fprintf('Wyniki zapisano do pliku: %s\n', nazwa_pliku_txt);

fprintf('Generowanie wykresów...\n'); %Generowanie wykresów
figure('Name', ['Analiza: ' target_constellation], 'Units', 'normalized', ...
    'Position', [0.05 0.1 0.9 0.6]);
sgtitle(sprintf('Analiza Żywotności: %s (Populacja: %d)', target_constellation, ...
    sum(~isnan(res_Z5))), 'FontSize', 16, 'FontWeight', 'bold');

valid_mask = ~isnan(res_Z5) & res_hp > 0;
hp_plot = res_hp(valid_mask); z1_plot = res_Z1(valid_mask); 
z2_plot = res_Z2(valid_mask); z5_plot = res_Z5(valid_mask);

max_y_val = max([max(z1_plot), max(z2_plot), max(z5_plot)]); % wspólny zakres osi y
if isempty(max_y_val) || isnan(max_y_val), max_y_val = 100; end
shared_ylim = [0, max_y_val * 1.1];

subplot(1, 3, 1);
Wykresy(hp_plot, z1_plot, 'b', 'Metoda 1 (King-Hele)', shared_ylim);
subplot(1, 3, 2);
Wykresy(hp_plot, z2_plot, 'm', 'Metoda 2 (Uproszczona)', shared_ylim);
subplot(1, 3, 3);
Wykresy(hp_plot, z5_plot, 'r', 'Metoda 3 (Numeryczna)', shared_ylim);

fprintf('Zakończono.\n');

function satelity = WczytajDane(nazwa_konstelacji, Const)
    
    plik_excel = 'Baza_Satelitow_LEO_2025.xlsx';
    nazwa_upper = upper(nazwa_konstelacji);
    szukana_nazwa_w_pliku = '';
    arkusz_docelowy = '';
    
    switch nazwa_upper
        case 'STARLINK'
            arkusz_docelowy = 'Starlinki';
        case 'ONEWEB'
            arkusz_docelowy = 'OneWeb';
        case 'KUIPER'
            arkusz_docelowy = 'Inne_Konstelacje'; szukana_nazwa_w_pliku = 'KUIPER';
        case 'QIANFAN'
            arkusz_docelowy = 'Inne_Konstelacje'; szukana_nazwa_w_pliku = 'QIANFAN';
        case {'GUOWANG', 'HULIANWANG'}
            arkusz_docelowy = 'Inne_Konstelacje'; szukana_nazwa_w_pliku = 'HULIANWANG';
        otherwise
            error('Nieznana konstelacja: %s', nazwa_konstelacji);
    end
    
    fprintf('  -> Wczytywanie arkusza "%s" z pliku %s...\n', arkusz_docelowy, plik_excel);
    
    opts = detectImportOptions(plik_excel, 'Sheet', arkusz_docelowy);
    opts.VariableNamingRule = 'preserve';
    
    T = readtable(plik_excel, opts);
    
    if ~isempty(szukana_nazwa_w_pliku)
        if ismember('Nazwa_TLE', T.Properties.VariableNames)
            idx = contains(string(T.Nazwa_TLE), szukana_nazwa_w_pliku, 'IgnoreCase', true);
            T = T(idx, :);
        else
            warning('Brak kolumny Nazwa_TLE w arkuszu. Nie można przefiltrować.');
        end
    end
    satelity = struct(); R_E = Const.R_E; mu = Const.mi;
    
    for k = 1:height(T)
        row = T(k, :);
        
        try
            if isdatetime(row.Czas_PL)
                satelity(k).t_date = row.Czas_PL;
            else
                satelity(k).t_date = datetime(row.Czas_PL, 'InputFormat', ...
                    'yyyy-MM-dd HH:mm:ss');
            end
        catch
            satelity(k).t_date = datetime('now');
        end
        
        satelity(k).e = row.Ekscentrycznosc; satelity(k).i = deg2rad(row.Inklinacja_deg);
        satelity(k).OMEGA = deg2rad(row.RAAN_deg); satelity(k).m = row.Mass;
        satelity(k).omega = deg2rad(row.Arg_Perygeum_deg);
        satelity(k).meanMotion = row.Obiegi_na_dzien;
        
        n_rad_s(k) = satelity(k).meanMotion * (2 * pi) / 86400;
        satelity(k).a = (mu / n_rad_s(k)^2)^(1/3);
        
        L = row.Length; D = row.Diameter; S = row.Span;
        if isnan(L), L = 1.0; end
        if isnan(D), D = 1.0; end
        if isnan(S), S = L; end

        satelity(k).A_SRP = L * S; satelity(k).A_D = satelity(k).A_SRP * 0.2;      
        satelity(k).Nazwa = row.Nazwa_TLE; satelity(k).ID = row.ID;

        valName = row.Nazwa_TLE;
        if iscell(valName)
            satelity(k).Nazwa = char(valName{1});
        elseif isstring(valName)
            satelity(k).Nazwa = char(valName);
        else
            satelity(k).Nazwa = char(valName);
        end
    end
end

function [Z1, Z2] = ObliczZywotnosci(a, e, n, r_p, C_B, Const)
    [rhop0, Hp0_km] = Gestosc(r_p, Const, Const.rhomax);
    
    Hp0 = Hp0_km * 1000;
    
    B = n / C_B * rhop0 * a * e * besseli(1, a * e / Hp0) * exp(-e * (1 + a / Hp0));
    if B ~= 0
        Z1s = e^2 * (1 - 11/6*e + 29/16*e^2 + 7/8*Hp0/a) / 2 / B;
        Z1 = Z1s / 3600 / 24 / 365;
    else
        Z1 = NaN;
    end
    Z2s = Hp0 * C_B / rhop0 / sqrt(Const.mi * a) * (1 - exp(-(a - ...
        Const.R_E)/Hp0) * (1 + (a - Const.R_E)/2/a));
    Z2 = Z2s / 3600 / 24 / 365;
end

function Z_num = ObliczZ5(a0, e0, i0, O0, w0, m, AD, ASRP, CD, CP, t0, rho_vec, C)
    
    ca = a0; ce = e0; ci = i0; cO = O0; cw = w0; ct = t0; elapsed_days = 0;
    
    max_iter = 100000; min_h = C.hpp; max_lifetime_years = 50000; rec_idx = 0;
    
    for iter = 1:max_iter
        hp = ca * (1 - ce) - C.R_E; rec_idx = rec_idx + 1;

        if hp < min_h, break; end
        if isnan(hp) || isinf(hp)
            fprintf('BŁĄD: Wykryto NaN w iteracji %d.\n', iter); break;
        end
        if elapsed_days > max_lifetime_years * 365.25
             Z_num = max_lifetime_years; break;
        end
        
        res = ObliczPerturbacje(ca, ce, ci, cO, cw, m, AD, ASRP, CD, CP, ct, C);
        
        da_dt = res.da_dtS; de_dt = res.de_dtS; 
        
        if abs(da_dt) < 1e-12, dt_a = 3650; 
        else, dt_a = abs(5000 / da_dt); end
        
        target_de = 0.002; 
        if abs(de_dt) < 1e-12, dt_e = 3650; 
        else, dt_e = abs(target_de / de_dt); end
        
        dt = min([dt_a, dt_e]);
        if dt > 60, dt = 60; end 
        if dt < 0.001, dt = 0.001; end
        
        next_ca = ca + da_dt * dt; next_ce = ce + de_dt * dt;
        
        if next_ce < 1e-12, next_ce = 1e-12; 
        elseif next_ce > 0.999, next_ce = 0.999; dt = dt * 0.5; end
        if next_ca < C.R_E, hp = 0; break; end
        
        ca = next_ca; ce = next_ce; ci = ci + deg2rad(res.di_dtS * dt);
        cO = cO + deg2rad(res.dOMEGA_dtS * dt); cw = cw + deg2rad(res.domega_dtS * dt);
        cO = mod(cO, 2*pi); cw = mod(cw, 2*pi);
        
        elapsed_days = elapsed_days + dt; ct = ct + days(dt);
    end
    Z_num = elapsed_days / 365.25;
end

function Wykresy(x, y, color, titleText, fixed_ylim)
    yyaxis left
    h_pts = scatter(x, y, 20, color, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.6);
    hold on; grid on;
    xlabel('Wysokość perigeum [km]');
    ylabel('Żywotność [Lata]');
    ylim(fixed_ylim);
    
    ax = gca; ax.YColor = 'k';
    h_trend = []; 
    
    if length(x) > 5
        [x_sort, idx] = sort(x);
        y_sort = y(idx);
        valid = y_sort > 0 & isfinite(y_sort);
        
        if sum(valid) > 5
            try
                f = fit(x_sort(valid), y_sort(valid), 'exp1');
                h_trend = plot(x_sort(valid), f(x_sort(valid)), 'k--', 'LineWidth', 1);
            catch
                p = polyfit(x_sort(valid), y_sort(valid), 2);
                h_trend = plot(x_sort, polyval(p, x_sort), 'k--', 'LineWidth', 1);
            end
        end
    end
    
    yyaxis right
    
    max_h = max(x);
    if isempty(max_h), max_h = 0; end
    
    if max_h <= 700
        bin_step = 50;
    else
        bin_step = 100;
    end
    
    start_bin = floor(min(x) / bin_step) * bin_step;
    end_bin = ceil(max(x) / bin_step) * bin_step;
    if start_bin == end_bin, end_bin = start_bin + bin_step; end
    
    edges = start_bin : bin_step : end_bin;
    
    [counts, ~] = histcounts(x, edges); total_count = length(x);
    percentages = (counts / total_count) * 100;
    
    centers = edges(1:end-1) + bin_step/2;
    
    h_bar = bar(centers, percentages, 0.8, 'FaceColor', [0.4 0.4 0.4], ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none');
    
    ylabel('Rozkład satelitów [%]');
    ax.YColor = [0.4 0.4 0.4];
    
    title(titleText);
    label_hist = sprintf('Rozkład (co %dkm)', bin_step);
    
    if ~isempty(h_trend)
        legend([h_pts, h_trend, h_bar], {'Żywotności', 'Linia trendu', ...
            label_hist}, 'Location', 'northwest');
    else
        legend([h_pts, h_bar], {'Satelity', label_hist}, 'Location', 'northwest');
    end
end

function [rho_val, H_i] = Gestosc(r_query_meters, Const, rho_input)
    
    R_E = Const.R_E; h0 = Const.h0; H = Const.H;
    h_query = (r_query_meters - R_E) / 1000; % [km]
    
    idx = find(h0 <= h_query, 1, 'last');
    
    if isempty(idx)
        idx = 1; 
    end
    
    if h_query == h0(idx)
        rho_val = rho_input(idx);
        return;
    end
    
    H_i    = H(idx); h0_i   = h0(idx); rho_base = rho_input(idx);
    
    rho_val = rho_base * exp(-(h_query - h0_i)/H_i);
end

function avg = ObliczPerturbacje(a, e, i, OMEGA, omega, m, A_D, A_SRP, C_D, C_P, t_date, Const)
    mi = Const.mi; R_E = Const.R_E; J_2 = Const.J_2; c = Const.c;
    epsilon_deg = Const.epsilon_deg; a_S = Const.a_S; r_S = Const.r_S;
    h0 = Const.h0; H = Const.H; rho0 = Const.rho0; rho01 = Const.rho01;
    
    epsilon = deg2rad(epsilon_deg);
    n = sqrt(mi / a ^ 3); n_deg_day = n * (180/pi) * 86400; n_obr = n * 3600 * 24 / 2 / pi; 
    
    if isempty(t_date.TimeZone), t_ref = datetime(year(t_date), 7, 4);
    else, t_ref = datetime(year(t_date), 7, 4, 'TimeZone', t_date.TimeZone); end
    D = 2 * pi / 365 * days(t_date - t_ref);
    day_of_year = day(t_date, 'dayofyear');
    lambda = deg2rad(mod((day_of_year - 80) * 360 / 365, 360));
    
    f_SRP = 1358 / (1.0004 + 0.0334 * cos(D)) / c * A_SRP / m * C_P * (a_S / r_S) ^ 2;
    
    dOMEGA_S = (-3 / 2 * J_2 * n * cos(i) / (1 - e ^ 2) ^ 2 * (R_E / a) ...
        ^ 2) * 24 * 3600 * 180 / pi;
    domega_S = (-3 / 4 * n * J_2 * (1 - 5 * cos(i) ^ 2) / (1 - e ^ ...
        2) ^ 2 * (R_E / a) ^ 2) * 24 * 3600 * 180 / pi;
    dM_perturb_S = (3 / 4 * n * J_2 * (3 * cos(i) ^ 2 - 1) / (1 - e ^ ...
        2) ^ (3 / 2) * (R_E / a) ^ 2) * 24 * 3600 * 180 / pi;
    dM_total_S   = n_deg_day + dM_perturb_S;
    
    dOMEGAmoon = -3.4e-3 * cos(i) / n_obr; domegamoon = 1.7e-3 * (5 * ...
        cos(i)^2 - 1) / n_obr;
    dOMEGAsun = -1.5e-3 * cos(i) / n_obr; domegasun = 0.8e-3 * (5 * ...
        cos(i)^2 - 1) / n_obr;
    dOMEGA_3body = dOMEGAmoon + dOMEGAsun; domega_3body = domegamoon + domegasun;
    
    theta_deg = 0:10:359; len_th = length(theta_deg);
    dadt_S = zeros(1, len_th); dedt_S = zeros(1, len_th); 
    didt_S = zeros(1, len_th); dOMEGAdt_S = zeros(1, len_th); 
    domegadt_S = zeros(1, len_th); dMdt_S = zeros(1, len_th);
    
    for k = 1:len_th
        theta_val = deg2rad(theta_deg(k)); u_val = omega + theta_val;
        r_val = a * (1 - e^2) / (1 + e * cos(theta_val)); 
        v_k = sqrt(mi * (2/r_val - 1/a));
        E_anom = 2 * atan(sqrt((1 - e)/(1 + e)) * tan(theta_val / 2));
        [rho_val, ~] = Gestosc(r_val, Const, Const.rhomax);
        F_D_val = 0.5 * rho_val * A_D * v_k ^ 2 * C_D;
        
        dadtatm = (- 2 * F_D_val / m / n / sqrt(1 - e ^ 2) * sqrt(e ^ 2 + ...
            1 + 2 * e * cos(theta_val))) * 86400;
        dedtatm = (- F_D_val / m / n / a / sqrt(e ^ 2 + 1 + 2 * e * ...
            cos(theta_val)) * sqrt(1 - e ^ 2) * (2 * e + 2 * ...
            cos(theta_val))) * 86400;
        domegadtatm = (- 2 * F_D_val * sin(theta_val) / m / n / a / e / ...
            sqrt(e ^ 2 + 1 + 2 * e * cos(theta_val)) * sqrt(1 - e ^ ...
            2)) * 86400 * 180 / pi;
        
        R_val = f_SRP * ( - cos(i/2)^2 * cos(epsilon/2)^2 * cos(lambda - ...
            (u_val) - OMEGA) - sin(i/2)^2 * sin(epsilon/2)^2 * cos(lambda - ...
            (u_val) + OMEGA) - 0.5 * sin(i) * sin(epsilon) * (cos(lambda - ...
            (u_val)) - cos(-lambda - (u_val))) - sin(i/2)^2 * cos(epsilon/ ...
            2)^2 * cos(-lambda - (u_val) + OMEGA) - cos(i/2)^2 * sin(epsilon/ ...
            2)^2 * cos(-lambda - (u_val) - OMEGA));
        T_val = f_SRP * ( - cos(i/2)^2 * cos(epsilon/2)^2 * sin(lambda - ...
            (u_val) - OMEGA) - sin(i/2)^2 * sin(epsilon/2)^2 * sin(lambda - ...
            (u_val) + OMEGA) - 0.5 * sin(i) * sin(epsilon) * (sin(lambda - ...
            (u_val)) - sin(-lambda - (u_val))) - sin(i/2)^2 * cos(epsilon/ ...
            2)^2 * sin(-lambda - (u_val) + OMEGA) - cos(i/2)^2 * sin(epsilon/ ...
            2)^2 * sin(-lambda - (u_val) - OMEGA));
        W_val = f_SRP * (sin(i) * cos(epsilon/2)^2 * sin(lambda - OMEGA) - ...
            sin(i) * sin(epsilon/2)^2 * sin(lambda + OMEGA) - cos(i) * ...
            sin(epsilon) * sin(lambda));
        
        dadtSRP = (2 / n / sqrt(1 - e ^ 2) * (e * sin(theta_val) * R_val ...
            + (1 + e * cos(theta_val))* T_val)) * 86400;
        dedtSRP = (sqrt(1 - e ^ 2) / n / a * (sin(theta_val) * R_val + ...
            (cos(E_anom) + cos(theta_val)) * T_val)) * 86400;
        didtSRP = (r_val / a ^ 2 / n / sqrt(1 - e ^ 2) * cos(theta_val + ...
            omega) * W_val) * 86400 * 180 / pi;
        dOMEGAdtSRP = (r_val / a ^ 2 / n / sqrt(1 - e ^ 2) * sin(theta_val ...
            + omega) / sin(i) * W_val) * 86400 * 180 / pi;
        domegadtSRP = (sqrt(1 - e ^ 2) / n / a / e * (- cos(theta_val) * ...
            R_val + (1 + 1 / (1 + e * cos(theta_val))) * sin(theta_val) * ...
            T_val - cos(i) * (dOMEGAdtSRP / 86400 / 180 * pi))) * (86400 * ...
            180 / pi);
        dMdtSRP = ((1 - e ^ 2) / n / a / e * ((-2 * e / (1 + e * ...
            cos(theta_val)) + cos(theta_val)) * R_val - (1 + 1 / ...
            (1 + e * cos(theta_val))) * T_val * sin(theta_val))) * (86400 * ...
            180 / pi);
        
        dadt_S(k) = dadtatm + dadtSRP; dedt_S(k) = dedtatm + dedtSRP; 
        didt_S(k) = didtSRP; dOMEGAdt_S(k) = dOMEGAdtSRP + dOMEGA_S + dOMEGA_3body;
        domegadt_S(k) = domegadtatm + domegadtSRP + domega_S + domega_3body;
        dMdt_S(k) = dM_total_S + dMdtSRP; 
    end
    avg.da_dtS = mean(dadt_S); avg.de_dtS = mean(dedt_S); 
    avg.di_dtS = mean(didt_S); avg.dOMEGA_dtS = mean(dOMEGAdt_S); 
    avg.domega_dtS = mean(domegadt_S); avg.dM_dtS = mean(dMdt_S);
end
