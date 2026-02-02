clc; clear; close all;
fprintf('=== INICJALIZACJA ===\n');

% 1. DANE I STAŁE
data = readmatrix("MSISE.txt");
Const.h0 = data(:,1); Const.H = data(:,2); Const.rho0 = data(:,3); 
Const.rho01 = data(:,4); Const.rhomin = data(:,5); Const.rhomid = data(:,6); 
Const.rhomax = data(:,7);

load('coeffs_X.mat');
ObliczZ3_handle = @(e_query, h_query) ObliczZ3(e_query, h_query, coeff1, ...
   deg1, mx1, sx1, my1, sy1, coeff2, deg2, mx2, sx2, my2, sy2);

Const.R_E = 6.378e6; Const.mi = 3.986e14; Const.J_2 = 1082.6e-6; Const.c = 2.9979e8;
Const.epsilon_deg = 23.44; Const.a_S = 1.496e11; Const.r_S = 1.496e11; 
Const.G = 6.674e-11; Const.n_sun = 1.99e-7; Const.n_moon = 2.66e-6; 
Const.i_sundeg = 23.49; Const.i_moondeg = 28.6 + 23.49; Const.hpp = 78000; 
Const.r_pp = Const.hpp + Const.R_E; Const.J_3 = -2.5327e-6; Const.J_4 = -1.6196e-6;

% 2. POBRANIE DANYCH
fprintf('Numery przypadków testowych: 38083, 43814, 27391, 43216, 37820 \n');
ID_Satelity = input('Podaj ID satelity (SATCAT): ');
test_IDs = [38083, 43814, 27391, 43216, 37820];

if ismember(ID_Satelity, test_IDs)
    sat = Dane_Excel(ID_Satelity, 'DeorbitowaneSatelity.xlsx', Const);
else
    sat = Dane_Full(ID_Satelity, Const);
end

m = sat.Masa;
len = sat.Wymiary_Raw(1); span = sat.Wymiary_Raw(3);
A_SRP = len*span; A_D = A_SRP*0.2;

h_p_km = sat.Orbita.h_per; h_a_km = sat.Orbita.h_apo;
i_deg = sat.Orbita.inc; OMEGA_deg = sat.Orbita.raan; omega_deg = sat.Orbita.argPer;
C_D = 2.5; C_P = 1.2; t_date = sat.Czas_PL;

r_p = h_p_km*1000 + Const.R_E; r_a = h_a_km*1000 + Const.R_E;
a = 0.5*(r_p + r_a); e = (r_a - r_p)/(r_a + r_p);
i = deg2rad(i_deg); OMEGA = deg2rad(OMEGA_deg); omega = deg2rad(omega_deg);
n = sqrt(Const.mi/a^3); C_B = m/A_D/C_D; T_S = 2*pi/sqrt(Const.mi)*a^(3/2);

% --- Metody 1, 2, 3, 4 ---
fprintf('\n=== ANALIZA WSTĘPNA ===\n');
[res_avg, res_vec] = ObliczPerturbacje(a, e, i, OMEGA, omega, m, A_D, A_SRP, ...
    C_D, C_P, t_date, Const.rho0, Const);

%WykresyPerturbacji(res_vec, sat, ID_Satelity);
%PorownajModele(res_vec);

[Z1, Z2, Z3, Z4] = ObliczZywotnosci(a, e, n, r_p, h_p_km*1000, C_B, T_S, ...
    res_avg, Const.rho0, Const, ObliczZ3_handle, 1);
[Z1min, Z2min, ~, ~] = ObliczZywotnosci(a, e, n, r_p, h_p_km*1000, C_B, T_S, ...
    res_avg, Const.rhomin, Const, ObliczZ3_handle, 0);
[Z1max, Z2max, ~, ~] = ObliczZywotnosci(a, e, n, r_p, h_p_km*1000, C_B, T_S, ...
    res_avg, Const.rhomax, Const, ObliczZ3_handle, 0);

% --- METODA 5 ---
[Z5, t_nom, h_nom, a_nom, e_nom] = ObliczZ5(a, e, i, OMEGA, ...
    omega, m, A_D, A_SRP, C_D, C_P, t_date, Const.rho0, Const, 1);
[Z5min, t_nommin, h_nommin, a_nommin, e_nommin] = ObliczZ5(a, ...
    e, i, OMEGA, omega, m, A_D, A_SRP, C_D, C_P, t_date, Const.rhomin, Const, 0);
[Z5max, t_nommax, h_nommax, a_nommax, e_nommax] = ObliczZ5(a, ...
    e, i, OMEGA, omega, m, A_D, A_SRP, C_D, C_P, t_date, Const.rhomax, Const, 0);

fprintf('\n=== ANALIZA CZUŁOŚCI ===');

Scenarios = {
    Const.rho0,   Z1,    Z2,    Z5,    t_nom,    h_nom,    a_nom,    (e_nom ...
    ),    'Średnia', 2;
    Const.rhomin, Z1min, Z2min, Z5min, t_nommin, h_nommin, a_nommin, (e_nommin ...
    ), 'Niska', 1;
    Const.rhomax, Z1max, Z2max, Z5max, t_nommax, h_nommax, a_nommax, (e_nommax ...
    ), 'Wysoka', 3 
};

order = [2, 1, 3]; 

for j = 1:3
    k = order(j);
    
    curr_rho  = Scenarios{k, 1};
    curr_desc = Scenarios{k, 9};
    plot_idx  = Scenarios{k, 10};
    
    fprintf('\n---> Rysowanie wykresu %d/3: %s aktywność słoneczna...', ...
        plot_idx, curr_desc);

    Z_struct_curr.Z1 = Scenarios{k, 2};
    Z_struct_curr.Z2 = Scenarios{k, 3};
    Z_struct_curr.Z3 = Z3;
    Z_struct_curr.Z4 = Z4;
    Z_struct_curr.Z5 = Scenarios{k, 4};

    Adapt_struct_curr.t_nom = Scenarios{k, 5};
    Adapt_struct_curr.h_nom = Scenarios{k, 6};
    Adapt_struct_curr.a_nom = Scenarios{k, 7};
    Adapt_struct_curr.e_nom = Scenarios{k, 8};

    WysokoscPerigeum(curr_rho, sat, ID_Satelity, Const, a, e, i, OMEGA, omega, ...
  m, A_D, A_SRP, C_D, C_P, t_date, Z_struct_curr, Adapt_struct_curr, curr_desc, plot_idx);
             
    if k == 1
       % AnalizaCzulosciParametrow(curr_rho, sat, ID_Satelity, Const, a, e, ...
   %i, OMEGA, omega, m, A_D, A_SRP, C_D, C_P, t_date, Z_struct_curr, Adapt_struct_curr);
    end
     
end

%WysokoscPerigeumZbiorcza(Scenarios, Z3, sat, ID_Satelity, Const, t_date);

fprintf('\n UKOŃCZONO CAŁĄ ANALIZĘ \n');

function WysokoscPerigeumZbiorcza(Scenarios, Z3_val, sat, ID_Satelity, Const, t_date)
    
    max_lifetime = 0;
    for k = 1:size(Scenarios, 1)
        vals = [Scenarios{k, 2}, Scenarios{k, 3}, Scenarios{k, 4}, Z3_val];
        max_lifetime = max(max_lifetime, max(vals(isfinite(vals))));
    end
    
    styles_rho = {'-', '--', ':'};
    labels_rho = {'Średnia', 'Niska', 'Wysoka'};
    
    fig_h = figure('Name', 'Analiza Zbiorcza - Wysokość Perigeum', ...
        'Units', 'normalized', 'Position', [0.05 0.05 0.9 0.85]); 
    hold on; grid on;
    
    title({['Analiza czułości zmiany wysokości perigeum dla kroku ' ...
        'adaptacyjnego dla różnej aktywności słonecznej.']; ...
        sprintf('Satelita %s (ID: %d)', sat.Nazwa, ID_Satelity)});
        
    ylabel('Wysokość perigeum [km]'); 
    xlabel('Dni od zakończenia manewrowości');
    
    yline(Const.hpp/1000, 'g:', 'LineWidth', 2.0, ...
        'Label', 'Granica deorbitacji', ...
        'LabelVerticalAlignment', 'bottom', ...
        'HandleVisibility', 'off'); 
        
    ylim([0 inf]);
    
    legend_entries = [];
    legend_txt = {};

    % RYSOWANIE KRZYWYCH ADAPTACYJNYCH (Czerwone)
    for k = 1:size(Scenarios, 1)
        style = styles_rho{k};
        desc  = labels_rho{k};
        
        t_ad = Scenarios{k, 5}; 
        h_ad = Scenarios{k, 6}; 
        
        if ~isempty(t_ad)
            d_ad = t_ad * 365.25; 
            p = plot(d_ad, h_ad, style, 'Color', 'r', 'LineWidth', 2.0);
            legend_entries = [legend_entries, p];
            legend_txt{end+1} = sprintf('%s aktywność słoneczna', desc);
        end
    end

    % SKALOWANIE CZASU (Dni -> Lata, jeśli > 3 lata)
    time_scale = 1;
    time_label = 'Dni od zakończenia manewrowości';
    if max_lifetime > 3
        time_scale = 365.25;
        time_label = 'Lata od zakończenia manewrowości';
    end
    
    if time_scale > 1
        all_lines = findobj(gca, 'Type', 'line');
        for L = all_lines'
            if max(L.XData) > (max_lifetime * 0.5) 
                L.XData = L.XData / time_scale;
            end
        end
        xlabel(time_label);
    end
    
    % DANE RZECZYWISTE (TLE) - Czarne
    if isfield(sat.Orbita, 'vec')
        d_real = days(sat.Orbita.vec.Czas - t_date) / time_scale;
        pr = plot(d_real, sat.Orbita.vec.h_per, 'k-', 'LineWidth', 2.5); 
        pr.Color = [0 0 0 0.8]; 
        
        legend_entries = [legend_entries, pr];
        legend_txt{end+1} = 'Rzeczywiste dane';
    end
    
    % LINIE ŻYWOTNOŚCI (Pionowe)
    v_color = [0.5 0.5 0.5];
    
    for k = 1:size(Scenarios, 1)
        if k == 1
            continue; 
        end
        
        zs = [Scenarios{k, 2}, Scenarios{k, 3}, Scenarios{k, 4}];
        names = {'Z1', 'Z2', 'Z5'};
        style = styles_rho{k};
        
        for z_idx = 1:length(zs)
            val = zs(z_idx);
            nm  = names{z_idx};
            
            if isfinite(val) && val > 0
                pos = (val * 365.25) / time_scale; 
                lbl_suffix = '';
                if k==2, lbl_suffix='_{min}'; elseif k==3, lbl_suffix='_{max}'; end
                
                h_align = 'left'; 
                if strcmp(nm, 'Z2')
                    h_align = 'right';
                else
                    h_align = 'left';
                end
                
                xline(pos, style, 'Color', v_color, 'LineWidth', 1.0, ...
                    'Label', [nm, lbl_suffix], 'LabelVerticalAlignment', ...
                    'bottom', 'LabelHorizontalAlignment', h_align, ...
                    'LabelOrientation', 'aligned', 'HandleVisibility', 'off');
            end
        end
    end
    
    if isfinite(Z3_val) && Z3_val > 0
        pos_z3 = (Z3_val * 365.25) / time_scale;
        xline(pos_z3, '-.', 'Color', v_color, 'LineWidth', 1.0, 'Label', 'Z3', ...
            'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', ...
            'left', 'LabelOrientation', 'aligned', 'HandleVisibility', 'off');
    end
    
    legend(legend_entries, legend_txt, 'Location', 'northeast');
    current_xlim = xlim;
    max_x_visible = (max_lifetime * 365.25) / time_scale;
    if max_x_visible > current_xlim(2)
        xlim([0, max_x_visible * 1.1]);
    end
end

function WysokoscPerigeum(rho_val, sat, ID_Satelity, Const, a, e, i, OMEGA, omega, ...
m, A_D, A_SRP, C_D, C_P, t_date, Z_struct, Adapt_struct, activity_label, subplot_idx)
    
    Z1 = Z_struct.Z1; Z2 = Z_struct.Z2; Z4 = Z_struct.Z4; Z5 = Z_struct.Z5;
    sim_days = max(400, (max([1.5*Z1, 1.5*Z2, Z4/20, 400/365]) * 365));
    lifetimes = [Z1, Z2]; 
    valid_Z = lifetimes(lifetimes > 0 & isfinite(lifetimes));
    if isempty(valid_Z), valid_Z = 1; end
    L_ref_days = min(valid_Z) * 365;
    dt_candidates = [0.1, 1, 10, 50, 100, 200]; 
    max_steps_limit = 300000; 
    dt_vec = dt_candidates( (sim_days ./ dt_candidates) <= max_steps_limit );
    styles = {'g', 'b', 'k', 'm', 'c', 'y'};
    legend_labels = {};
    plots_h = []; 
    fig_name = 'Czułość - Wysokość Perigeum (Multi-View)';
 
    if subplot_idx == 1
        fig_h = figure('Name', fig_name, 'Units', 'normalized', ['' ...
            'Position'], [0.05 0.15 0.9 0.6]);
        sgtitle(sprintf(['Analiza czułości zmiany wysokości perigeum dla satelity ' ...
        '%s (ID: %d)'], sat.Nazwa, ID_Satelity), 'FontWeight', 'bold', 'FontSize', 14);
    else
        fig_h = findobj('Type', 'figure', 'Name', fig_name);
        if isempty(fig_h)
            fig_h = figure('Name', fig_name, 'Units', 'normalized', ...
                'Position', [0.05 0.15 0.9 0.6]);
        else
            figure(fig_h); 
        end
    end
    
    subplot(1, 3, subplot_idx);
    hold on; grid on;
    title(sprintf('%s aktywność słoneczna', activity_label), 'FontWeight', 'bold');
    ylabel('Wysokość perigeum [km]'); 
    xlabel('Dni od zakończenia manewrowości');
    yline(Const.hpp/1000, 'r--', 'LineWidth', 1.5, 'Label', '78km', ...
        'HandleVisibility', 'off'); 
    ylim([0 inf]);
    
    for idx = 1:length(dt_vec)
        dt_val = dt_vec(idx);
        sim_days_sensitivity = max(sim_days, L_ref_days * 2.5);
        
        [t_sim, h_sim, ~, ~, ~, ~, ~, ~] = Symulacja(a, e, i, OMEGA, omega, m, ...
            A_D, A_SRP, C_D, C_P, t_date, rho_val, Const, dt_val, sim_days_sensitivity);
        
        if ~isempty(t_sim)
            d_sim = days(t_sim - t_sim(1));
            curr_style = styles{mod(idx-1, length(styles)) + 1};
            
            p = plot(d_sim, h_sim/1000, curr_style, 'LineWidth', 1.2);
            plots_h = [plots_h, p]; 
        end
        legend_labels{end+1} = sprintf('dt = %.4g d', dt_val);
    end
    
    time_scale = 1;
    time_label = 'Dni od zakończenia manewrowości';
    if exist('Z5', 'var') && Z5 > 3
        time_scale = 365.25;
        time_label = 'Lata od zakończenia manewrowości';
    end
    
    if time_scale > 1
        for ph = plots_h
            if isvalid(ph) && max(ph.XData) > (Z5 * 1.5) 
                ph.XData = ph.XData / time_scale;
            end
        end
        xlabel(gca, time_label);
    end
    
    % --- DANE RZECZYWISTE ---
    if isfield(sat.Orbita, 'vec')
        d_real = days(sat.Orbita.vec.Czas - t_date) / time_scale;
        p_real = plot(d_real, sat.Orbita.vec.h_per, '-.', 'Color', ...
            'k', 'LineWidth', 0.8);
        plots_h = [plots_h, p_real];
        legend_labels{end+1} = 'Rzeczywiste wartości';
    end
    
    % --- METODA ADAPTACYJNA ---
    legend_labels_main = legend_labels;
    if isfield(Adapt_struct, 't_nom') && ~isempty(Adapt_struct.t_nom)
        d_nom = (Adapt_struct.t_nom * 365.25) / time_scale;
        p_nom = plot(d_nom, Adapt_struct.h_nom, 'r-', 'LineWidth', 2);
        plots_h = [plots_h, p_nom];
        legend_labels_main{end+1} = 'Krok adaptacyjny';
    end
    
    % --- LINIE PIONOWE ---
    line_labels = {'Z1 (King-Hele)', 'Z2 (Simple)', 'Z3 (Empir.)','Z5 (Num.)'};
    vals_Z = [Z1, Z2, Z_struct.Z3, Z5]; 
    line_color = [0.4 0.4 0.4 0.6]; 
    max_val_to_show = 0;
    
    for k = 1:length(vals_Z)
        L_yrs = vals_Z(k); 
        if isfinite(L_yrs) && L_yrs > 0
            pos_val = (L_yrs * 365.25) / time_scale;
            xline(pos_val, '-.', line_labels{k}, 'Color', line_color, ...
                'LineWidth', 1.2, 'LabelVerticalAlignment', 'bottom', ...
                'LabelHorizontalAlignment', 'right', 'HandleVisibility', 'off');
            max_val_to_show = max(max_val_to_show, pos_val);
        end
    end
    
    if max_val_to_show > xlim(gca).*(1)
       xlim([0, max(xlim(gca).*(1), max_val_to_show * 1.1)]);
    end
    
    if subplot_idx == 3
        lgd = legend(plots_h, legend_labels_main, 'Orientation', 'horizontal');
        set(lgd, 'Units', 'normalized');
        set(lgd, 'Position', [0.1, 0.02, 0.8, 0.05]); 
    end
end

function AnalizaCzulosciParametrow(rho_val, sat, ID_Satelity, Const, a, e, ...
    i, OMEGA, omega, m, A_D, A_SRP, C_D, C_P, t_date, Z_struct, Adapt_struct)
    
    Z1 = Z_struct.Z1; Z2 = Z_struct.Z2; Z4 = Z_struct.Z4; Z5 = Z_struct.Z5;
    
    sim_days = max(400, (max([1.5*Z1, 1.5*Z2, Z4/20, 400/365]) * 365));
    lifetimes = [Z1, Z2];
    valid_Z = lifetimes(lifetimes > 0 & isfinite(lifetimes));
    if isempty(valid_Z), valid_Z = 1; end
    L_ref_days = min(valid_Z) * 365;
    dt_candidates = [0.1, 1, 10, 50, 100, 200];
    max_steps_limit = 300000;
    dt_vec = dt_candidates( (sim_days ./ dt_candidates) <= max_steps_limit );
    styles = {'g', 'b', 'k', 'm', 'c', 'y'};
    legend_labels = {};
    plots_legend_handles = [];
    
    fig_elem = figure('Name', 'Czułość - Elementy Orbitalne', 'Units', ...
        'normalized', 'Position', [0.1 0.15 0.8 0.7]);
    sgtitle(sprintf('Analiza czułości parametrów orbitalnych satelity %s (ID: %d)', ...
        sat.Nazwa, ID_Satelity), 'FontWeight', 'bold', 'FontSize', 14);
    
    ax_a = subplot(2,3,1); hold on; title('Półoś wielka a [km]'); 
    grid on; xlabel('Dni'); ylim(ax_a, [6400 inf]);
    ax_e = subplot(2,3,2); hold on; title('Mimośród e [-]'); 
    grid on; xlabel('Dni');
    ax_i = subplot(2,3,3); hold on; title('Inklinacja i [deg]'); 
    grid on; xlabel('Dni');
    ax_w = subplot(2,3,4); hold on; title(['Argument perigeum \omega [deg]']); 
    grid on; xlabel('Dni');
    ax_O = subplot(2,3,5); hold on; title('RAAN \Omega [deg]'); 
    grid on; xlabel('Dni');
    ax_M = subplot(2,3,6); hold on; title('Anomalia Średnia M [deg]'); 
    grid on; xlabel('Dni');
    
    for idx = 1:length(dt_vec)
        dt_val = dt_vec(idx);
        sim_days_sensitivity = max(sim_days, L_ref_days * 2.5);
        
        [t_sim, ~, a_sim, e_sim, i_sim, O_sim, w_sim, M_sim] = Symulacja(a, ...
            e, i, OMEGA, omega, m, A_D, A_SRP, C_D, C_P, t_date, rho_val, ...
            Const, dt_val, sim_days_sensitivity);
        
        if ~isempty(t_sim)
            d_sim = days(t_sim - t_sim(1));
            curr_style = styles{mod(idx-1, length(styles)) + 1};
            p = plot(ax_a, d_sim, a_sim/1000, curr_style, 'LineWidth', 1.2);
            plot(ax_e, d_sim, e_sim, curr_style, 'LineWidth', 1.2);
            plot(ax_i, d_sim, i_sim, curr_style, 'LineWidth', 1.2);
            plot(ax_w, d_sim, w_sim, curr_style, 'LineWidth', 1.2);
            plot(ax_O, d_sim, O_sim, curr_style, 'LineWidth', 1.2);
            plot(ax_M, d_sim, M_sim, curr_style, 'LineWidth', 1.2);
            plots_legend_handles = [plots_legend_handles, p];
        end
        legend_labels{end+1} = sprintf('dt = %.4g d', dt_val);
    end
    
    time_scale = 1;
    time_label = 'Dni';
    if exist('Z5', 'var') && Z5 > 3
        time_scale = 365.25;
        time_label = 'Lata';
    end
    
    all_axes = [ax_a, ax_e, ax_i, ax_w, ax_O, ax_M];
    if time_scale > 1
        for ax = all_axes
            lines = findobj(ax, 'Type', 'line');
            for L = lines'
                if max(L.XData) > (Z5 * 1.5)
                    L.XData = L.XData / time_scale;
                end
            end
            xlabel(ax, time_label);
        end
    end
    
    % --- DANE RZECZYWISTE ---
    if isfield(sat.Orbita, 'vec')
        d_real = days(sat.Orbita.vec.Czas - t_date) / time_scale;
        real_style = '-.'; real_color = 'k'; real_width = 0.8;
        
        p_real = plot(ax_a, d_real, sat.Orbita.vec.semiMajorAxis, ...
            real_style, 'Color', real_color, 'LineWidth', real_width);
        plot(ax_e, d_real, sat.Orbita.vec.ecc, real_style, 'Color', ...
            real_color, 'LineWidth', real_width);
        plot(ax_i, d_real, sat.Orbita.vec.inc, real_style, 'Color', ...
            real_color, 'LineWidth', real_width);
        plot(ax_w, d_real, sat.Orbita.vec.argPer, real_style, 'Color', ...
            real_color, 'LineWidth', real_width);
        plot(ax_O, d_real, sat.Orbita.vec.raan, real_style, 'Color', ...
            real_color, 'LineWidth', real_width);
        plot(ax_M, d_real, zeros(size(d_real)), real_style, 'Color', ...
            real_color, 'LineWidth', real_width);
        
        plots_legend_handles = [plots_legend_handles, p_real];
        legend_labels{end+1} = 'Rzeczywiste';
    end
    
    % --- METODA ADAPTACYJNA ---
    if isfield(Adapt_struct, 't_nom') && ~isempty(Adapt_struct.t_nom)
         d_nom = (Adapt_struct.t_nom * 365.25) / time_scale;
         
         if isfield(Adapt_struct, 'a_nom') && isfield(Adapt_struct, 'e_nom')
             adapt_style = 'r-'; width = 1.5;
             p_nom = plot(ax_a, d_nom, Adapt_struct.a_nom, adapt_style, ...
                 'LineWidth', width);
             plot(ax_e, d_nom, Adapt_struct.e_nom, adapt_style, ...
                 'LineWidth', width);
             
             plots_legend_handles = [plots_legend_handles, p_nom];
             legend_labels{end+1} = 'Krok adaptacyjny';
         end
    end
    
    ylim(ax_w, [0 360]); yticks(ax_w, 0:90:360);
    ylim(ax_O, [0 360]); yticks(ax_O, 0:90:360);
    ylim(ax_M, [0 360]); yticks(ax_M, 0:90:360);
    
    lgd = legend(plots_legend_handles, legend_labels, 'Orientation', 'horizontal');
    set(lgd, 'Units', 'normalized');
    set(lgd, 'Position', [0.1, 0.02, 0.8, 0.05]);
   
end

function [t_out, h_out, a_out, e_out, i_out, O_out, w_out, M_out] = Symulacja( ...
    a0, e0, i0, O0, w0, m, AD, ASRP, CD, CP, t0, rho_val, C, dt, max_days)
    est_steps = ceil(max_days/dt) + 10;
    if isempty(t0.TimeZone)
        t_hist = NaT(est_steps, 1);
    else
        t_hist = NaT(est_steps, 1, 'TimeZone', t0.TimeZone);
    end
    
    h_hist = zeros(est_steps, 1); a_hist = zeros(est_steps, 1); 
    e_hist = zeros(est_steps, 1); i_hist = zeros(est_steps, 1);
    O_hist = zeros(est_steps, 1); w_hist = zeros(est_steps, 1);
    M_hist = zeros(est_steps, 1);
    ca = a0; ce = e0; ci = i0; cO = O0; cw = w0; ct = t0; cM = 0;
    step = 1; elapsed = 0;
    
    while elapsed <= max_days
        hp = ca * (1 - ce) - C.R_E;
        
        if step > length(h_hist)
            block = 1000;
            h_hist = [h_hist; zeros(block,1)];
            a_hist = [a_hist; zeros(block,1)];
            e_hist = [e_hist; zeros(block,1)];
            i_hist = [i_hist; zeros(block,1)];
            O_hist = [O_hist; zeros(block,1)];
            w_hist = [w_hist; zeros(block,1)];
            M_hist = [M_hist; zeros(block,1)];
            
            if isempty(ct.TimeZone)
                t_hist = [t_hist; NaT(block,1)];
            else
                t_hist = [t_hist; NaT(block,1,'TimeZone',ct.TimeZone)];
            end
        end
        
        t_hist(step) = ct; h_hist(step) = hp; a_hist(step) = ca;
        e_hist(step) = ce; i_hist(step) = rad2deg(ci); O_hist(step) = rad2deg(cO);
        w_hist(step) = rad2deg(cw); M_hist(step) = rad2deg(mod(cM, 2*pi));
        
        if hp < C.hpp, break; end
        if ca <= C.R_E, break; end
        
        [res, ~] = ObliczPerturbacje(ca, ce, ci, cO, cw, m, AD, ASRP, CD, ...
            CP, ct, rho_val, C);
        
        ca = ca + res.da_dtS * dt; ce = ce + res.de_dtS * dt;
        ci = ci + deg2rad(res.di_dtS * dt); cO = cO + deg2rad(res.dOMEGA_dtS * dt);
        cw = cw + deg2rad(res.domega_dtS * dt); cM = cM + deg2rad(res.dM_dtS * dt);
        
        if ce < 1e-12, ce = 1e-12; end
        if ce >= 1, break; end
        
        ct = ct + days(dt);
        elapsed = elapsed + dt;
        step = step + 1;
    end
    
    idx = 1:step;
    t_out = t_hist(idx); h_out = h_hist(idx);
    a_out = a_hist(idx); e_out = e_hist(idx);
    i_out = mod(i_hist(idx), 360); O_out = mod(O_hist(idx), 360);
    w_out = mod(w_hist(idx), 360); M_out = mod(M_hist(idx), 360);
end

function [rho_val, H_i] = Gestosc(r_query_meters, Const, rho_input)
    
    R_E = Const.R_E; h0 = Const.h0; H = Const.H;
    h_query = (r_query_meters - R_E) / 1000;
    
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

function satelliteData = Dane_Full(satID, Const)
  
    if nargin < 1 || isempty(satID)
        fprintf('Wprowadzanie interaktywne.\n');
        satID = input('Podaj ID satelity (SATCAT, np. 25544): ');
    end
    
    fprintf('--- POBIERANIE DANYCH DLA SATELITY ID: %d ---\n', satID);
    
    satName = "Nieznany"; 
    line1 = ""; line2 = "";
    found_tle_local = false;
    
    %% Orbital Data (TLE)
    local_tle_file = 'active_tle_catalog.txt';
    
    % Sprawdzenie bazy lokalnej
    if exist(local_tle_file, 'file')
        try
      
            fileText = fileread(local_tle_file); % Wczytujemy cały plik jako tekst
            
            % Szukamy wzorca dla Linii 1: "1 25544U" (spacja + ID + litera)
            % ID musi mieć 5 cyfr, więc używamy %05d
            pattern = sprintf('1 %05d', satID);
            idx = strfind(fileText, pattern);
            
            if ~isempty(idx)
                % Musimy wyciągnąć linię przed (Nazwa) i linię po (Linia 2)
                % Dzielimy fragment tekstu na linie
                % Bierzemy margines +/- 200 znaków wokół znaleziska
                startPos = max(1, idx(1) - 100);
                endPos = min(length(fileText), idx(1) + 200);
                chunk = fileText(startPos:endPos);
                chunkLines = splitlines(chunk);
                
                % Szukamy linii w wycinku
                localIdx = find(contains(chunkLines, pattern), 1);
                
                if ~isempty(localIdx) && localIdx > 1 && length(chunkLines) > localIdx
                    satName = strtrim(chunkLines{localIdx-1});
                    line1 = chunkLines{localIdx};
                    line2 = chunkLines{localIdx+1};
                    
                    % Dodatkowa weryfikacja czy Linia 2 pasuje do ID
                    if startsWith(strtrim(line2), sprintf('2 %05d', satID))
                        found_tle_local = true;
                        fprintf('   Znaleziono dane orbitalne.\n');
                    end
                end
            end
        catch
            fprintf('   Błąd odczytu pliku lokalnego. Poszukiwanie online.\n');
        end
    else
        fprintf('   Brak lokalnej bazy TLE (%s).\n', local_tle_file);
    end
    
    % Pobieranie z Internetu (jeśli brak lokalnie)
    if ~found_tle_local
        fprintf('1. Pobieranie danych orbitalnych z CelesTrak (ONLINE)...\n');
        url_tle = sprintf( ...
            'https://celestrak.org/NORAD/elements/gp.php?CATNR=%d&FORMAT=TLE', ...
            satID);
        
        try
            options = weboptions('Timeout', 30, 'UserAgent', 'Mozilla/5.0');
            tle_text = webread(url_tle, options);
            lines = splitlines(strtrim(tle_text));
            
            if length(lines) >= 3
                satName = strtrim(lines{1});
                line1 = lines{2};
                line2 = lines{3};
            else
                error('Otrzymano puste lub błędne TLE.');
            end
        catch
            error('Błąd pobierania TLE (możliwy ban IP lub złe ID).');
        end
    end
    
    % --- PARSOWANIE DANYCH TLE ---
    try
        epochStr = line1(19:32); 
        yy = str2double(epochStr(1:2));
        if yy < 57, year = 2000 + yy; else, year = 1900 + yy; end
        doy = str2double(epochStr(3:end));
        t_utc = datetime(year, 1, 1, 'TimeZone', 'UTC') + days(doy - 1);
        t_pl = t_utc; t_pl.TimeZone = 'Europe/Warsaw';
        
        inc = str2double(line2(9:16));
        raan = str2double(line2(18:25));
        ecc = str2double(['0.' line2(27:33)]);
        argPer = str2double(line2(35:42));
        meanMotion = str2double(line2(53:63));
        
        mu = Const.mi;
        n_rad_s = meanMotion * (2 * pi) / 86400;
        a_m = (mu / n_rad_s^2)^(1/3);
        a_km = a_m / 1000;
        h_per = a_km * (1 - ecc) - Const.R_E/1000;
        h_apo = a_km * (1 + ecc) - Const.R_E/1000;
    catch
        error('Błąd parsowania linii TLE. Dane mogą być uszkodzone.');
    end
    
    %% Dane Fizyczne (GCAT)
    fprintf('2. Analiza katalogu fizycznego GCAT...\n');
    url_gcat = 'https://planet4589.org/space/gcat/tsv/cat/satcat.tsv';
    local_gcat_file = 'satcat_physical.tsv';
    needDownload = true;
    if exist(local_gcat_file, 'file')
        fprintf('   Używam lokalnego katalogu %s.\n', local_gcat_file);
        needDownload = false;
    else
        fprintf('   Brak pliku %s.\n', local_gcat_file);
    end
    
    if needDownload
        fprintf('   Pobieranie katalogu (ok. 15-20 MB) na dysk...\n');
        try
            options = weboptions('Timeout', 60, 'UserAgent', 'Mozilla/5.0');
            websave(local_gcat_file, url_gcat, options);
            fprintf('   Pobrano pomyślnie.\n');
        catch
            fprintf('   Błąd pobierania katalogu!\n');
        end
    end
    
    mass = NaN; len = NaN; dia = NaN; span = NaN; source = "Brak";
    
    % Skanowanie pliku w poszukiwaniu naszego satelity
    if exist(local_gcat_file, 'file')
        fid = fopen(local_gcat_file, 'r');
        
        % 1. Znalezienie nagłówka i mapowanie kolumn
        idx_ID = 2; idx_Mass = 16; idx_Len = 22; idx_Dia = 24; 
        idx_Span = 26; 
        headerFound = false;
        
        % Skanujemy max 100 pierwszych linii w poszukiwaniu nagłówka
        for i=1:100
            line = fgetl(fid);
            if ~ischar(line), break; end
            if contains(line, 'Satcat') && contains(line, 'Mass')
                headers = split(string(line), sprintf('\t'));
                for k = 1:length(headers)
                    hName = strtrim(headers(k));
                    if hName == "Satcat", idx_ID = k; end
                    if hName == "Mass", idx_Mass = k; end
                    if hName == "Length", idx_Len = k; end
                    if hName == "Diameter", idx_Dia = k; end
                    if hName == "Span", idx_Span = k; end
                end
                headerFound = true;
                break;
            end
        end
        
        % 2. Szukanie konkretnego ID (Skanowanie liniowe)
        % Nie wczytujemy całego pliku do RAM, tylko lecimy linijka po linijce
        if headerFound
            while ~feof(fid)
                line = fgetl(fid);
                if ~ischar(line), break; end
                
                % Szybki check czy linia zawiera ID
                if contains(line, num2str(satID))
                    parts = split(string(line), sprintf('\t'));
                    
                    % Sprawdzamy czy to właściwa kolumna
                    if length(parts) >= idx_ID
                        valID = str2double(parts(idx_ID));
                        if valID == satID
                            cleanVal = @(val) str2double(replace(val, "-", "NaN"));
                            
                            if length(parts) >= idx_Mass, mass = cleanVal(parts( ...
                                    idx_Mass)); end
                            if length(parts) >= idx_Len, len = cleanVal(parts( ...
                                    idx_Len)); end
                            if length(parts) >= idx_Dia, dia = cleanVal(parts( ...
                                    idx_Dia)); end
                            if length(parts) >= idx_Span, span = cleanVal(parts( ...
                                    idx_Span)); end
                            
                            if ~isnan(mass) || ~isnan(len)
                                source = "Katalog Jonathana McDowella (Lokalny)";
                            else
                                source = "Katalog GCAT";
                            end
                            break;
                        end
                    end
                end
            end
        else
            fprintf('   Nie udało się zmapować kolumn w pliku GCAT.\n');
        end
        fclose(fid);
    else
        fprintf('   Brak pliku katalogu. Pomijam dane fizyczne.\n');
    end

    %% Wyświetlenie Wyników
    fprintf('\n=========================================\n');
    fprintf(' SATELITA: %s (ID: %d)\n', satName, satID);
    fprintf('=========================================\n');
    fprintf('Czas danych (PL):       %s\n', datestr(t_pl, 'yyyy-mm-dd HH:MM:SS'));
    fprintf('-----------------------------------------\n');
    
    if ~isnan(mass)
        fprintf('Masa:                   %.2f kg\n', mass);
    else
        fprintf('Masa:                   Brak danych w katalogu\n');
    end
    
    gotDims = false;
    if ~isnan(len), fprintf('Długość (Length):       %.2f m\n', len); 
        gotDims = true; end
    if ~isnan(dia), fprintf('Średnica (Diameter):    %.2f m\n', dia); 
        gotDims = true; end
    if ~isnan(span), fprintf('Rozpiętość (Span):      %.2f m\n', span); 
        gotDims = true; end
    
    if ~gotDims
        fprintf('Wymiary:                Brak danych w katalogu\n');
    end
    
    fprintf('Źródło danych fiz.:     %s\n', source);
    fprintf('-----------------------------------------\n');
    fprintf('Ekscentryczność:        %.7f\n', ecc);
    fprintf('Inklinacja:             %.4f deg\n', inc);
    fprintf('Wys. Perigeum:          %.2f km\n', h_per);
    fprintf('Wys. Apogeum:           %.2f km\n', h_apo);
    fprintf('RAAN:                   %.4f deg\n', raan);
    fprintf('Arg. Perigeum:          %.4f deg\n', argPer);
    fprintf('Obiegi dzienne:         %.8f\n', meanMotion);
    fprintf('=========================================\n');

    %% Zapis do struktury wyjściowej
    satelliteData.Nazwa = satName; satelliteData.Masa = mass;
    satelliteData.Wymiary_Raw = [len, dia, span]; satelliteData.Czas_PL = t_pl;
    
    satelliteData.Orbita.ecc = ecc; satelliteData.Orbita.inc = inc;
    satelliteData.Orbita.raan = raan; satelliteData.Orbita.argPer = argPer;
    satelliteData.Orbita.meanMotion = meanMotion;
    satelliteData.Orbita.h_per = h_per; satelliteData.Orbita.h_apo = h_apo;
end

function satelliteData = Dane_Excel(satID, excelFileName, Const)

    %% Wczytanie danych z Excela
    sheetName = num2str(satID); 
       
    % Wczytujemy z zachowaniem oryginalnych nagłówków
        opts = detectImportOptions(excelFileName, 'Sheet', sheetName);
        opts.VariableNamingRule = 'preserve';
        rawTable = readtable(excelFileName, opts);

    %% Przetwarzanie danych
    satNameStr = string(rawTable.OBJECT_NAME{1});
    
    % Masa i Wymiary
        mass_vec = rawTable.Masa; len_vec = rawTable.Length;
        dia_vec = rawTable.Diameter; span_vec = rawTable.Span;

    t_pl_vec = datetime(rawTable.EPOCH, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSSSSS');
    
    ecc_vec = rawTable.ECCENTRICITY; inc_vec = rawTable.INCLINATION;
    raan_vec = rawTable.RA_OF_ASC_NODE; argPer_vec = rawTable.ARG_OF_PERICENTER;
    meanMotion_vec = rawTable.MEAN_MOTION; n_rad_s = meanMotion_vec .* (2*pi / 86400);
    a_vec = ((Const.mi ./ n_rad_s.^2).^(1/3)) / 1000;
    r_per_vec = a_vec .* (1 - ecc_vec); r_apo_vec = a_vec .* (1 + ecc_vec);
    h_per_vec = r_per_vec - Const.R_E / 1000; h_apo_vec = r_apo_vec - Const.R_E / 1000;

    %% Wyświetlenie
    fprintf('\n=========================================\n');
    fprintf(' SATELITA: %s (ID: %d)\n', satNameStr, satID);
    fprintf('=========================================\n');
    fprintf('Zakres danych:          %s  ->  %s\n', datestr(min( ...
        t_pl_vec)), datestr(max(t_pl_vec)));
    fprintf('-----------------------------------------\n');
    fprintf('Masa:                   %.2f kg\n', mass_vec(1));
    fprintf('Długość (Length):       %.2f m\n', len_vec(1));
    fprintf('Średnica (Diameter):    %.2f m\n', dia_vec(1));
    fprintf('Rozpiętość (Span):      %.2f m\n', span_vec(1));
    fprintf('-----------------------------------------\n');
    fprintf('Ekscentryczność:        %.7f\n', ecc_vec(1));
    fprintf('Inklinacja:             %.4f deg\n', inc_vec(1));
    fprintf('Wys. Perigeum:          %.2f km\n', h_per_vec(1));
    fprintf('Wys. Apogeum:           %.2f km\n', h_apo_vec(1));
    fprintf('RAAN:                   %.4f deg\n', raan_vec(1));
    fprintf('Arg. Perigeum:          %.4f deg\n', argPer_vec(1));
    fprintf('Obiegi dzienne:         %.8f\n', meanMotion_vec(1));
    fprintf('=========================================\n');

    %% Zapis do struktury wyjściowej
    satelliteData.Nazwa = satNameStr; satelliteData.Masa = mass_vec(1);
    satelliteData.Wymiary_Raw = [len_vec(1), dia_vec(1), span_vec(1)];
    satelliteData.Czas_PL = t_pl_vec(1);
    
    satelliteData.Orbita.ecc = ecc_vec(1); satelliteData.Orbita.inc = inc_vec(1);
    satelliteData.Orbita.raan = raan_vec(1); satelliteData.Orbita.argPer = argPer_vec(1);
    satelliteData.Orbita.meanMotion = meanMotion_vec(1);
    satelliteData.Orbita.h_per = h_per_vec(1); satelliteData.Orbita.h_apo = h_apo_vec(1);
    
    satelliteData.Orbita.vec.Czas = t_pl_vec; satelliteData.Orbita.vec.ecc = ecc_vec;
    satelliteData.Orbita.vec.inc = inc_vec; satelliteData.Orbita.vec.raan = raan_vec;
    satelliteData.Orbita.vec.argPer = argPer_vec;
    satelliteData.Orbita.vec.meanMotion = meanMotion_vec;
    satelliteData.Orbita.vec.h_per = h_per_vec; satelliteData.Orbita.vec.h_apo = h_apo_vec;
    satelliteData.Orbita.vec.semiMajorAxis = a_vec;
end

function [Zywotnosc_1, Zywotnosc_2, Zywotnosc_3, Zywotnosc_4] = ObliczZywotnosci( ...
    a, e, n, r_p, h_p, C_B, T_S, avg, rho_val, Const, ObliczZ3_handle, wynik)
    
    mi = Const.mi; R_E=Const.R_E; h0=Const.h0; H=Const.H; 
    rpp = Const.r_pp; dadt = avg.da_dtS; dedt = avg.de_dtS; 
    
    % METODA 1 - King-Hele
    [rhop0, Hp0_km] = Gestosc(r_p, Const, rho_val); 
    Hp0 = Hp0_km * 1000;
    
    B = n / C_B * rhop0 * a * e * besseli(1, a * e / Hp0) * exp(-e ...
        * (1 + a / Hp0));
    
    if B ~= 0
        Zywotnosc_1s = e ^ 2 * (1 - 11 / 6 * e + 29 / 16 * e ^ 2 + 7 / 8 * ...
            Hp0 / a) / 2 / B; 
    else
        Zywotnosc_1s = Inf;
    end
    Zywotnosc_1 = Zywotnosc_1s / 3600 / 24 / 365;
    
    % METODA 2 - Uproszczona
    if e > 0.005
        Zywotnosc_2 = NaN;
    else
        Zywotnosc_2s = Hp0 * C_B / rhop0 / sqrt(mi * a) * (1 - ...
            exp(-(a - R_E) / Hp0) * (1 + (a - R_E) / 2 / a)); 
        Zywotnosc_2 = Zywotnosc_2s / 3600 / 24 / 365;
    end
 
    % METODA 3 - Empiryczna
    if h_p >= 150000 && h_p <= 550000
        NC = ObliczZ3_handle(e, h_p / 1000);
        Zywotnosc_3s = NC * C_B * T_S; 
        Zywotnosc_3 = Zywotnosc_3s / 3600 / 24 / 365;
    else
        Zywotnosc_3 = NaN;
    end
    
    % METODA 4 - Autorska, Uproszczona
    dadt_srednia = mean(dadt); dedt_srednia = mean(dedt);
    
    a_delta = dadt_srednia * dedt_srednia;
    b_delta = a * dedt_srednia + dadt_srednia * (1 - e);
    c_delta = r_p - rpp;
    Delta = b_delta ^ 2 - 4 * a_delta * c_delta;
    
    if Delta >= 0
        Zywotnosc_4x = (-b_delta - sqrt(Delta)) / 2 / a_delta;
        Zywotnosc_4s = (-b_delta + sqrt(Delta)) / 2 / a_delta;
        Zywotnosc_4 = max(Zywotnosc_4x, Zywotnosc_4s) / 365;
    else
        Zywotnosc_4 = NaN;
    end

    % WYŚWIETLANIE
    if wynik == 1
        fprintf('\n--- WYNIKI OBLICZEŃ ŻYWOTNOŚCI ---\n');
        if Zywotnosc_1 < 1
            fprintf('Metoda 1 (Dziki Wzór):      %.1f dni\n', Zywotnosc_1 * 365);
        else
            fprintf('Metoda 1 (Dziki Wzór):      %.4f lat\n', Zywotnosc_1);
        end
        if e > 0.005
            fprintf(['Metoda 2 (Uproszczona):     Mimośrodowość e nie jest ' ...
                'pomijalna. Metoda poza zakresem stosowalności (e>0,005).\n']);
        elseif Zywotnosc_2 < 1
            fprintf('Metoda 2 (Uproszczona):     %.1f dni\n', Zywotnosc_2 * 365);
        else
            fprintf('Metoda 2 (Uproszczona):     %.4f lat\n', Zywotnosc_2);
        end
        if h_p >= 150000 && h_p <= 550000
            if Zywotnosc_3 < 1
                fprintf('Metoda 3 (Empiryczna):      %.1f dni\n', Zywotnosc_3 * 365);
            else
                fprintf('Metoda 3 (Empiryczna):      %.4f lat\n', Zywotnosc_3);
            end
        else
            fprintf(['Metoda 3 (Empiryczna):      Wysokość h_p poza ' ...
                'zakresem stosowalności (150-550 km).\n']);
        end
    end
end

function WykresyPerturbacji(vec, sat, ID_Satelity)
    
    theta_deg = vec.theta_deg; dadt = vec.dadt; dadtatm = vec.dadtatm;
    dadtSRP = vec.dadtSRP; dedt = vec.dedt; dedtatm = vec.dedtatm;
    dedtSRP = vec.dedtSRP; dedtniejed = vec.geo_simple.de;
    dedtniejedC = vec.geo_complex.de; dedtC = vec.dedtC;
    didt = vec.didt; didtSRP = vec.didtSRP; didtniejed = vec.geo_simple.di;
    didtniejedC = vec.geo_complex.di; didtC = vec.didtC; dMdt = vec.dMdt;
    dMdtSRP = vec.dMdtSRP; dMdtniejed = vec.geo_simple.dM; 
    dMdtniejedC = vec.geo_complex.dM; dMdtC = vec.dMdtC; domegadt = vec.domegadt;
    domegadtatm = vec.domegadtatm; domegadtSRP = vec.domegadtSRP;
    domeganiejed = vec.geo_simple.domega; domegasun = vec.sun_simple.domega;
    domegamoon = vec.moon_simple.domega; domegadtniejedC = vec.geo_complex.domega;
    domegadtC = vec.domegadtC; dOMEGAdt = vec.dOMEGAdt; dOMEGAdtSRP = vec.dOMEGAdtSRP;
    dOMEGAniejed = vec.geo_simple.dOMEGA; dOMEGAsun = vec.sun_simple.dOMEGA; 
    dOMEGAmoon = vec.moon_simple.dOMEGA; dOMEGAdtniejedC = vec.geo_complex.dOMEGA;
    dOMEGAdtC = vec.dOMEGAdtC;

    c_SRP  = 'r'; c_Atm  = 'b'; c_Geo  = 'g'; c_Sun  = [0.93, 0.69, 0.13];        
    c_Moon = [0.45, 0.45, 0.45]; c_Tot  = 'k';                       
    s_Comp = '--'; s_Tot  = '-'; s_C = '-.'; w_Comp = 1.2; w_Tot  = 1.5;  
    x_pos = 350;

    fprintf('Generowanie wykresów perturbacji...\n');

    %% Wykres 1: Półoś wielka (a)
    figure;
    hold on;
    plot(theta_deg, dadtatm, 'Color', c_Atm, 'LineStyle', s_Comp, 'LineWidth', w_Comp);
    plot(theta_deg, dadtSRP, 'Color', c_SRP, 'LineStyle', s_Comp, 'LineWidth', w_Comp);
    plot(theta_deg, dadt, 'Color', c_Tot, 'LineStyle', s_Tot, 'LineWidth', w_Tot);
    yline(0, 'Color', [0.3 0.3 0.3], 'LineWidth', 0.5, 'HandleVisibility', 'off'); 
    
    xlabel('\theta [deg]');
    xlim([0 360]);
    xticks(0:30:360);
    ylabel('da/dt [ m / day ]');
    title('Zmiana półosi wielkiej w czasie');
    legend('Opór', 'SRP', 'Łącznie');
    grid on;

    %% Wykres 2: Mimośród (e)
    figure;
    hold on;
    plot(theta_deg, dedtatm, 'Color', c_Atm, 'LineStyle', s_Comp, 'LineWidth', w_Comp);
    plot(theta_deg, dedtSRP, 'Color', c_SRP, 'LineStyle', s_Comp, 'LineWidth', w_Comp);
    plot(theta_deg, ones(size(theta_deg)) * dedtniejed, 'Color', ...
        c_Geo, 'LineStyle', s_Comp, 'LineWidth', w_Comp);
    plot(theta_deg, dedt, 'Color', c_Tot, 'LineStyle', s_Tot, 'LineWidth', w_Tot);
    plot(theta_deg, ones(size(theta_deg)) * dedtniejedC, 'Color', ...
        c_Geo, 'LineStyle', s_C, 'LineWidth', w_Comp);
    plot(theta_deg, dedtC, 'Color', c_Tot, 'LineStyle', s_C, 'LineWidth', w_Tot);
    yline(0, 'Color', [0.3 0.3 0.3], 'LineWidth', 0.5, 'HandleVisibility', 'off'); 
    text(x_pos, dedtniejed, ' Model J2', 'Color', c_Geo, 'VerticalAlignment', ...
        'bottom', 'HorizontalAlignment', 'right', 'FontSize', 9, 'FontWeight', 'bold');
    text(x_pos, dedtniejedC, ' Model J2-J4', 'Color', c_Geo, 'VerticalAlignment', ...
        'top', 'HorizontalAlignment', 'right', 'FontSize', 9, 'FontWeight', 'bold');

    xlabel('\theta [deg]');
    xlim([0 360]);
    xticks(0:30:360);
    ylabel('de/dt [ 1 / day ]');
    title('Zmiana mimośrodu w czasie');
    legend('Opór', 'SRP', 'Niejednorodność (model J2)', ...
        'Łącznie (model J2)', 'Niejednorodność (model J2, J3, J4)', ...
        'Łącznie (model J2, J3, J4)' );
    grid on;

    %% Wykres 3: Inklinacja (i)
    figure;
    hold on;
    plot(theta_deg, didtSRP, 'Color', c_SRP, 'LineStyle', s_Comp, 'LineWidth', w_Comp);
    plot(theta_deg, ones(size(theta_deg)) * didtniejed, 'Color', ...
        c_Geo, 'LineStyle', s_Comp, 'LineWidth', w_Comp);
    plot(theta_deg, didt, 'Color', c_Tot, 'LineStyle', s_Tot, 'LineWidth', w_Tot);
    plot(theta_deg, ones(size(theta_deg)) * didtniejedC, 'Color', ...
        c_Geo, 'LineStyle', s_C, 'LineWidth', w_Comp);
    plot(theta_deg, didtC, 'Color', c_Tot, 'LineStyle', s_C, 'LineWidth', w_Tot);
    yline(0, 'Color', [0.3 0.3 0.3], 'LineWidth', 0.5, 'HandleVisibility', 'off');
    text(x_pos, didtniejed, ' Model J2', 'Color', c_Geo, 'VerticalAlignment', ...
        'bottom', 'HorizontalAlignment', 'right', 'FontSize', 9, 'FontWeight', 'bold');
    text(x_pos, didtniejedC, ' Model J2-J4', 'Color', c_Geo, 'VerticalAlignment', ...
        'top', 'HorizontalAlignment', 'right', 'FontSize', 9, 'FontWeight', 'bold');

    xlabel('\theta [deg]');
    xlim([0 360]);
    xticks(0:30:360);
    ylabel('di/dt [ deg / day ]');
    title('Zmiana inklinacji w czasie');
    legend('SRP', 'Niejednorodność (model J2)', 'Łącznie (model J2)', ...
        'Niejednorodność (model J2, J3, J4)', 'Łącznie (model J2, J3, J4)');
    grid on;

    %% Wykres 4: Anomalia średnia (M)
    figure;
    hold on;
    plot(theta_deg, dMdtSRP, 'Color', c_SRP, 'LineStyle', s_Comp, 'LineWidth', w_Comp);
    plot(theta_deg, ones(size(theta_deg)) * dMdtniejed, 'Color', c_Geo, ...
        'LineStyle', s_Comp, 'LineWidth', w_Comp);
    plot(theta_deg, dMdt, 'Color', c_Tot, 'LineStyle', s_Tot, 'LineWidth', w_Tot);
    plot(theta_deg, ones(size(theta_deg)) * dMdtniejedC, 'Color', c_Geo, ...
        'LineStyle', s_C, 'LineWidth', w_Comp);
    plot(theta_deg, dMdtC, 'Color', c_Tot, 'LineStyle', s_C, 'LineWidth', w_Tot);

    yline(0, 'Color', [0.3 0.3 0.3], 'LineWidth', 0.5, 'HandleVisibility', 'off');
    
    text(x_pos, dMdtniejed, ' Model J2', 'Color', c_Geo, 'VerticalAlignment', ...
        'bottom', 'HorizontalAlignment', 'right', 'FontSize', 9, 'FontWeight', 'bold');
    text(x_pos, dMdtniejedC, ' Model J2-J4', 'Color', c_Geo, 'VerticalAlignment', ...
        'top', 'HorizontalAlignment', 'right', 'FontSize', 9, 'FontWeight', 'bold');

    xlabel('\theta [deg]');
    xlim([0 360]);
    xticks(0:30:360);
    ylabel('dM/dt - n [ deg / day ]');
    title('Zmiana anomalii średniej w czasie');
    legend('SRP', 'Niejednorodność (model J2)', 'Łącznie (model J2)', ...
        'Niejednorodność (model J2, J3, J4)', 'Łącznie (model J2, J3, J4)');
    grid on;

    %% Wykres 5: Argument pericentrum (omega)
    figure;
    hold on;
    plot(theta_deg, domegadtatm, 'Color', c_Atm, 'LineStyle', s_Comp, 'LineWidth', w_Comp);
    plot(theta_deg, domegadtSRP, 'Color', c_SRP, 'LineStyle', s_Comp, 'LineWidth', w_Comp);
    plot(theta_deg, ones(size(theta_deg)) * domegasun, 'Color', ...
        c_Sun, 'LineStyle', s_Comp, 'LineWidth', w_Comp);
    plot(theta_deg, ones(size(theta_deg)) * domegamoon, 'Color', ...
        c_Moon, 'LineStyle', s_Comp, 'LineWidth', w_Comp);
    plot(theta_deg, ones(size(theta_deg)) * domeganiejed, 'Color', ...
        c_Geo, 'LineStyle', s_Comp, 'LineWidth', w_Comp);
    plot(theta_deg, domegadt, 'Color', c_Tot, 'LineStyle', s_Tot, 'LineWidth', w_Tot);
    plot(theta_deg, ones(size(theta_deg)) * domegadtniejedC, 'Color', ...
        c_Geo, 'LineStyle', s_C, 'LineWidth', w_Comp);
    plot(theta_deg, domegadtC, 'Color', c_Tot, 'LineStyle', s_C, 'LineWidth', w_Tot);
    yline(0, 'Color', [0.3 0.3 0.3], 'LineWidth', 0.5, 'HandleVisibility', 'off');
    text(x_pos, domeganiejed, ' Model J2', 'Color', c_Geo, 'VerticalAlignment', ...
        'bottom', 'HorizontalAlignment', 'right', 'FontSize', 9, 'FontWeight', 'bold');
    text(x_pos, domegadtniejedC, ' Model J2-J4', 'Color', c_Geo, 'VerticalAlignment', ...
        'top', 'HorizontalAlignment', 'right', 'FontSize', 9, 'FontWeight', 'bold');

    xlabel('\theta [deg]');
    xlim([0 360]);
    xticks(0:30:360);
    ylabel('d\omega/dt [ deg / day ]');
    title('Zmiana argumentu pericentrum w czasie');
    legend('Opór', 'SRP', 'Grawitacja Słońca', 'Grawitacja Księżyca', ...
        'Niejednorodność (model J2)', 'Łącznie (model J2)', ...
        'Niejednorodność (model J2, J3, J4)', 'Łącznie (model J2, J3, J4)');
    grid on;

    %% Wykres 6: Węzeł wstępujący (OMEGA)
    figure;
    hold on;
    plot(theta_deg, dOMEGAdtSRP, 'Color', c_SRP, 'LineStyle', s_Comp, 'LineWidth', w_Comp);
    plot(theta_deg, ones(size(theta_deg)) * dOMEGAsun, 'Color', ...
        c_Sun, 'LineStyle', s_Comp, 'LineWidth', w_Comp);
    plot(theta_deg, ones(size(theta_deg)) * dOMEGAmoon, 'Color', ...
        c_Moon, 'LineStyle', s_Comp, 'LineWidth', w_Comp);
    plot(theta_deg, ones(size(theta_deg)) * dOMEGAniejed, 'Color', ...
        c_Geo, 'LineStyle', s_Comp, 'LineWidth', w_Comp);
    plot(theta_deg, dOMEGAdt, 'Color', c_Tot, 'LineStyle', s_Tot, 'LineWidth', w_Tot);
    plot(theta_deg, ones(size(theta_deg)) * dOMEGAdtniejedC, 'Color', ...
        c_Geo, 'LineStyle', s_C, 'LineWidth', w_Comp);
    plot(theta_deg, dOMEGAdtC, 'Color', c_Tot, 'LineStyle', s_C, 'LineWidth', w_Tot);
    yline(0, 'Color', [0.3 0.3 0.3], 'LineWidth', 0.5, 'HandleVisibility', 'off');

    text(x_pos, dOMEGAniejed, ' Model J2', 'Color', c_Geo, 'VerticalAlignment', ...
        'bottom', 'HorizontalAlignment', 'right', 'FontSize', 9, 'FontWeight', 'bold');
    text(x_pos, dOMEGAdtniejedC, ' Model J2-J4', 'Color', c_Geo, 'VerticalAlignment', ...
        'top', 'HorizontalAlignment', 'right', 'FontSize', 9, 'FontWeight', 'bold');
    
    xlabel('\theta [deg]');
    xlim([0 360]);
    xticks(0:30:360);
    ylabel('d\Omega/dt [ deg / day ]');
    title('Zmiana węzła wstępującego w czasie');
    legend('SRP', 'Grawitacja Słońca', 'Grawitacja Księżyca', ...
        'Niejednorodność (model J2)', 'Łącznie (model J2)', ...
        'Niejednorodność (model J2, J3, J4)', 'Łącznie (model J2, J3, J4)');
    grid on;

    %% Jedno okno z 6 wykresami (TOTAL)
    figure;
        sgtitle(sprintf(['Całkowite perturbacje parametrów orbitalnych satelity ' ...
        '%s (ID: %d)'], sat.Nazwa, ID_Satelity), 'FontWeight', 'bold', 'FontSize', 14);
    subplot(2, 3, 1);
    plot(theta_deg, dadt);
    title('da/dt total [ m / day ]');
    xlabel('\theta [deg]');
    xlim([0 360]);
    xticks(0:30:360);
    grid on;
    yline(mean(dadt), '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1, ...
        'HandleVisibility', 'off'); 
    if min(dadt) <= 0 && max(dadt) >= 0
        yline(0, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1, ...
         'HandleVisibility', 'off');
    end

    subplot(2, 3, 2);
    plot(theta_deg, dedt);
    title('de/dt total [ - / day ]');
    xlabel('\theta [deg]');
    xlim([0 360]);
    xticks(0:30:360);
    grid on;
    yline(mean(dedt), '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1, ...
        'HandleVisibility', 'off');
    if min(dedt) <= 0 && max(dedt) >= 0
        yline(0, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1, ...
            'HandleVisibility', 'off');
    end

    subplot(2, 3, 3);
    plot(theta_deg, didt);
    title('di/dt total [ deg / day ]');
    xlabel('\theta [deg]');
    xlim([0 360]);
    xticks(0:30:360);
    grid on;
    yline(mean(didt), '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1, ...
        'HandleVisibility', 'off');
    if min(didt) <= 0 && max(didt) >= 0
        yline(0, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1, ...
            'HandleVisibility', 'off');
    end

    subplot(2, 3, 4);
    plot(theta_deg, dOMEGAdt);
    title('d\Omega/dt total [ deg / day ]');
    xlabel('\theta [deg]');
    xlim([0 360]);
    xticks(0:30:360);
    grid on;
    yline(mean(dOMEGAdt), '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1, ...
        'HandleVisibility', 'off');
    if min(dOMEGAdt) <= 0 && max(dOMEGAdt) >= 0
        yline(0, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1, ...
            'HandleVisibility', 'off');
    end

    subplot(2, 3, 5);
    plot(theta_deg, domegadt);
    title('d\omega/dt total [ deg / day ]');
    xlabel('\theta [deg]');
    xlim([0 360]);
    xticks(0:30:360);
    grid on;
    yline(mean(domegadt), '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1, ...
        'HandleVisibility', 'off');
    if min(domegadt) <= 0 && max(domegadt) >= 0
        yline(0, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1, ...
            'HandleVisibility', 'off');
    end

    subplot(2, 3, 6);
    plot(theta_deg, dMdt);
    title('dM/dt - n total [ deg / day ]');
    xlabel('\theta [deg]');
    xlim([0 360]);
    xticks(0:30:360);
    grid on;
    yline(mean(dMdt), '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1, ...
        'HandleVisibility', 'off');
    if min(dMdt) <= 0 && max(dMdt) >= 0
        yline(0, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1, ...
            'HandleVisibility', 'off');
    end
end

function [avg, vec] = ObliczPerturbacje(a, e, i, OMEGA, omega, m, ...
    A_D, A_SRP, C_D, C_P, t_date, rho_vec, Const)
    
    mi = Const.mi; R_E = Const.R_E; J_2 = Const.J_2; c = Const.c;
    epsilon_deg = Const.epsilon_deg; a_S = Const.a_S; r_S = Const.r_S;
    
    if isempty(t_date.TimeZone), t_ref = datetime(year(t_date), 7, 4);
    else, t_ref = datetime(year(t_date), 7, 4, 'TimeZone', t_date.TimeZone); end
    D = 2 * pi / 365 * days(t_date - t_ref);
    day_of_year = day(t_date, 'dayofyear');
    lambda = deg2rad(mod((day_of_year - 80) * 360 / 365, 360));
    epsilon = deg2rad(epsilon_deg);
    f_SRP = 1358 / (1.0004 + 0.0334 * cos(D)) / c * A_SRP / m * C_P * (a_S / r_S) ^ 2;
    n = sqrt(mi / a ^ 3); n_deg_day = n * (180/pi) * 86400;

    % Model Uproszczony Niejednorodności
    dOMEGA_S = (-3 / 2 * J_2 * n * cos(i) / (1 - e ^ 2) ^ 2 * (R_E / a) ^ ...
        2) * 86400 * 180 / pi;
    domega_S = (-3 / 4 * n * J_2 * (1 - 5 * cos(i) ^ 2) / (1 - e ^ ...
        2) ^ 2 * (R_E / a) ^ 2) * 86400 * 180 / pi;
    dM_perturb_S = (3 / 4 * n * J_2 * (3 * cos(i) ^ 2 - 1) / (1 - e ^ ...
        2) ^ (3 / 2) * (R_E / a) ^ 2) * 86400 * 180 / pi;
    dM_total_S = dM_perturb_S;
    
    % Model Złożony Niejednorodności
    [de_C, di_C, dOMEGA_C, domega_C, dM_perturb_C] = ObliczNiejednorodnosc( ...
        a, e, i, omega, n, Const);
    dM_total_C = dM_perturb_C;

    % Perturbacje Słońca/Księżyca
    n_obr = n * 3600 * 24 / 2 / pi;
    dOMEGAmoon = -3.4e-3 * cos(i) / n_obr; domegamoon = 1.7e-3 * (5 * ...
        cos(i)^2 - 1) / n_obr;
    dOMEGAsun = -1.5e-3 * cos(i) / n_obr; domegasun = 0.8e-3 * (5 * ...
        cos(i)^2 - 1) / n_obr;
    i_sun = deg2rad(Const.i_sundeg); i_moon = deg2rad(Const.i_moondeg);
    dOMEGA3bodysun = (-3/8*Const.n_sun^2/n*(1+1.5*e^2)/sqrt(1-e^2)*(3* ...
        cos(i_sun)^2-1)*cos(i))*86400*180/pi;
    domega3bodysun = (3/8*Const.n_sun^2/n*(1-1.5*sin(i_sun)^2)/sqrt(1-e^ ...
        2)*(5*cos(i)^2-1+e^2))*86400*180/pi;
    dOMEGA3bodymoon = (-3/8*Const.n_moon^2/n*(1+1.5*e^2)/sqrt(1-e^2)*(3* ...
        cos(i_moon)^2-1)*cos(i))*86400*180/pi;
    domega3bodymoon = (3/8*Const.n_moon^2/n*(1-1.5*sin(i_moon)^2)/sqrt(1-e^ ...
        2)*(5*cos(i)^2-1+e^2))*86400*180/pi;

    % ATMOSFERA i SRP
    theta_deg = 0:10:359; len_th = length(theta_deg);
    dadt = zeros(1, len_th); dadtatm = zeros(1, len_th); 
    dadtSRP = zeros(1, len_th); dedt = zeros(1, len_th); 
    dedtatm = zeros(1, len_th); dedtSRP = zeros(1, len_th);
    didtSRP = zeros(1, len_th); dOMEGAdtSRP = zeros(1, len_th); 
    domegadtSRP = zeros(1, len_th); dMdtSRP = zeros(1, len_th);

    for k = 1:len_th
        theta_val = deg2rad(theta_deg(k));
        u_val = omega + theta_val;
        r_val = a * (1 - e^2) / (1 + e * cos(theta_val));
        v_k = sqrt(mi * (2/r_val - 1/a));
        E_anom = 2 * atan(sqrt((1 - e)/(1 + e)) * tan(theta_val / 2));
        [rho_val,~] = Gestosc(r_val, Const, rho_vec); 
        F_D_val = 0.5 * rho_val * A_D * v_k ^ 2 * C_D;
        
        dadtatm(k) = (- 2 * F_D_val / m / n / sqrt(1 - e ^ 2) * sqrt(e ^ ...
            2 + 1 + 2 * e * cos(theta_val))) * 86400;
        dedtatm(k) = (- F_D_val / m / n / a / sqrt(e ^ 2 + 1 + 2 * e * ...
            cos(theta_val)) * sqrt(1 - e ^ 2) * (2 * e + 2 * ...
            cos(theta_val))) * 86400;
        domegadtatm(k) = (- 2 * F_D_val * sin(theta_val) / m / n / a / e / ...
            sqrt(e ^ 2 + 1 + 2 * e * cos(theta_val)) * sqrt(1 - e ^ ...
            2)) * 86400 * 180 / pi;
        
        % SRP
        R_val = f_SRP * ( - cos(i/2)^2 * cos(epsilon/2)^2 * ...
            cos(lambda - (u_val) - OMEGA) - sin(i/2)^2 * ...
            sin(epsilon/2)^2 * cos(lambda - (u_val) + OMEGA) - 0.5 * ...
            sin(i) * sin(epsilon) * (cos(lambda - (u_val)) - ...
            cos(-lambda - (u_val))) - sin(i/2)^2 * cos(epsilon/2)^2 * ...
            cos(-lambda - (u_val) + OMEGA) - cos(i/2)^2 * sin(epsilon/2)^2 * ...
            cos(-lambda - (u_val) - OMEGA));
        T_val = f_SRP * ( - cos(i/2)^2 * cos(epsilon/2)^2 * sin(lambda - ...
            (u_val) - OMEGA) - sin(i/2)^2 * sin(epsilon/2)^2 * sin(lambda - ...
            (u_val) + OMEGA) - 0.5 * sin(i) * sin(epsilon) * (sin(lambda - ...
            (u_val)) - sin(-lambda - (u_val))) - sin(i/2)^2 * ...
            cos(epsilon/2)^2 * sin(-lambda - (u_val) + OMEGA) - ...
            cos(i/2)^2 * sin(epsilon/2)^2 * sin(-lambda - (u_val) - OMEGA));
        W_val = f_SRP * (sin(i) * cos(epsilon/2)^2 * sin(lambda - OMEGA) - ...
            sin(i) * sin(epsilon/2)^2 * sin(lambda + OMEGA) - cos(i) * ...
            sin(epsilon) * sin(lambda));
        
        dadtSRP(k) = (2 / n / sqrt(1 - e ^ 2) * (e * sin(theta_val) * ...
            R_val + (1 + e * cos(theta_val))* T_val)) * 86400;
        dedtSRP(k) = (sqrt(1 - e ^ 2) / n / a * (sin(theta_val) * R_val + ...
            (cos(E_anom) + cos(theta_val)) * T_val)) * 86400;
        didtSRP(k) = (r_val / a ^ 2 / n / sqrt(1 - e ^ 2) * ...
            cos(theta_val + omega) * W_val) * 86400 * 180 / pi;
        dOMEGAdtSRP(k) = (r_val / a ^ 2 / n / sqrt(1 - e ^ 2) * ...
            sin(theta_val + omega) / sin(i) * W_val) * 86400 * 180 / pi;
        domegadtSRP(k) = (sqrt(1 - e ^ 2) / n / a / e * (- ...
            cos(theta_val) * R_val + (1 + 1 / (1 + e * cos(theta_val))) * ...
            sin(theta_val) * T_val - cos(i) * (dOMEGAdtSRP(k) / 86400 / ...
            180 * pi))) * 86400 * 180 / pi;
        dMdtSRP(k) = ((1 - e ^ 2) / n / a / e * ((-2 * e / (1 + e * ...
            cos(theta_val)) + cos(theta_val)) * R_val - (1 + 1 / (1 + e * ...
            cos(theta_val))) * T_val * sin(theta_val))) * 86400 * 180 / pi;
        
        % Sumowanie
        dadt_S(k) = dadtatm(k) + dadtSRP(k); dedt_S(k) = dedtatm(k) + dedtSRP(k);
        didt_S(k) = didtSRP(k); dMdt_S(k) = dM_total_S + dMdtSRP(k);
        dOMEGAdt_S(k) = dOMEGAdtSRP(k) + dOMEGA_S + dOMEGA3bodysun + dOMEGA3bodymoon;
        domegadt_S(k) = (domegadtatm(k) + domegadtSRP(k) + domega_S + ...
            domega3bodysun + domega3bodymoon);
        
        % Sumowanie Complex
        dadt_C(k) = dadtatm(k) + dadtSRP(k); dedt_C(k) = dedtatm(k) + dedtSRP(k) + de_C;
        didt_C(k) = didtSRP(k) + di_C; dMdt_C(k) = dM_total_C + dMdtSRP(k);
        dOMEGAdt_C(k) = dOMEGAdtSRP(k) + dOMEGA_C + dOMEGA3bodysun + dOMEGA3bodymoon;
        domegadt_C(k) = (domegadtatm(k) + domegadtSRP(k) + domega_C + ...
            domega3bodysun + domega3bodymoon);   
    end
    
    avg.da_dtS = mean(dadt_S); avg.de_dtS = mean(dedt_S); avg.di_dtS = mean(didt_S); 
    avg.dOMEGA_dtS = mean(dOMEGAdt_S); avg.domega_dtS = mean(domegadt_S); 
    avg.dM_dtS = mean(dMdt_S); avg.dM_dtC = mean(dMdt_C);
    
    avg.da_dtC = mean(dadt_C); avg.de_dtC = mean(dedt_C); avg.di_dtC = mean(didt_C); 
    avg.dOMEGA_dtC = mean(dOMEGAdt_C); avg.domega_dtC = mean(domegadt_C); 

    % Wektory dla wykresów
    vec.theta_deg = theta_deg; vec.dadt = dadt_S; vec.dadtatm = dadtatm; 
    vec.dadtSRP = dadtSRP; vec.dedt = dedt_S; vec.dedtatm = dedtatm; 
    vec.dedtSRP = dedtSRP; vec.didt = didt_S; vec.didtSRP = didtSRP; 
    vec.dMdt = dMdt_S; vec.dMdtSRP = dMdtSRP; vec.domegadt = domegadt_S; 
    vec.domegadtatm = domegadtatm; vec.domegadtSRP = domegadtSRP;
    vec.dOMEGAdt = dOMEGAdt_S; vec.dOMEGAdtSRP = dOMEGAdtSRP;
    vec.geo_simple.de = 0; vec.geo_simple.di = 0; 
    vec.geo_simple.dOMEGA = dOMEGA_S; vec.geo_simple.domega = domega_S; 
    vec.geo_simple.dM = dM_total_S; vec.geo_complex.de = de_C; 
    vec.geo_complex.di = di_C; vec.geo_complex.dOMEGA = dOMEGA_C; 
    vec.geo_complex.domega = domega_C; vec.geo_complex.dM = dM_total_C;
    vec.dOMEGAdtC= dOMEGAdt_C; vec.domegadtC= domegadt_C; vec.dMdtC= dMdt_C; 
    vec.didtC= didt_C; vec.dedtC= dedt_C; vec.dadtC= dadt_C;  
    vec.sun_simple.domega = domegasun; vec.sun_simple.dOMEGA = dOMEGAsun;
    vec.moon_simple.domega = domegamoon; vec.moon_simple.dOMEGA = dOMEGAmoon;
end

function X = ObliczZ3(e_query, h_query, coeff1, deg1, mx1, sx1, my1, ...
    sy1, coeff2, deg2, mx2, sx2, my2, sy2)

X = zeros(size(e_query));

for k = 1:numel(e_query)
    e = e_query(k);
    h = h_query(k);

    if e <= 0.1
        xn = (h - mx1)/sx1;
        yn = (e - my1)/sy1;
        A = Wielomian([xn],[yn], deg1);
        log10N = A * coeff1(:);
    else
        xn = (h - mx2)/sx2;
        yn = (e - my2)/sy2;
        A = Wielomian([xn],[yn], deg2);
        log10N = A * coeff2(:);
    end
    X(k) = 10.^log10N;
end
end

function A = Wielomian(x, y, deg)
count = 0;
for d = 0:deg
    for i = 0:d
        j = d - i;
        count = count + 1;
        A(:,count) = (x(:).^i) .* (y(:).^j);
    end
end
end

function [dedt, didt, dOmegadt, domegadt, dMdt_perturb] = ObliczNiejednorodnosc( ...
    a, e, i, omega, n, Const)

    R_E = Const.R_E; J_2 = Const.J_2; J_3 = Const.J_3; J_4 = Const.J_4;
    
    p = a * (1 - e^2); RE_p = R_E / p; sin_i = sin(i); sin2_i = sin(i)^2;
    sin4_i = sin(i)^4; cos_i = cos(i); sin_2i = sin(2*i); sin_omega = sin(omega);
    cos_omega = cos(omega); sin_2omega = sin(2*omega); cos_2omega = cos(2*omega);
    cos_4omega = cos(4*omega); sqrt_1_e2 = sqrt(1 - e^2);
    
    % de/dt
    term1 = -3/32 * n * J_2^2 * (RE_p)^4 * sin2_i * (14 - 15* ...
        sin2_i) * e * (1 - e^2) * sin_2omega;
    term2 = -3/8 * n * J_3 * (RE_p)^3 * sin_i * (4 - 5* ...
        sin2_i) * (1 - e^2) * cos_omega;
    term3 = -15/32 * n * J_4 * (RE_p)^4 * sin2_i * (6 - 7* ...
        sin2_i) * e * (1 - e^2) * sin_2omega;
    dedt_rad_s = term1 + term2 + term3;

    % di/dt
    term1 = 3/64 * n * J_2^2 * (RE_p)^4 * sin_2i * (14 - 15* ...
        sin2_i) * e^2 * sin_2omega;
    term2 = 3/8 * n * J_3 * (RE_p)^3 * cos_i * (4 - 5* ...
        sin2_i) * e * cos_omega;
    term3 = 15/64 * n * J_4 * (RE_p)^4 * sin_2i * (6 - 7* ...
        sin2_i) * e^2 * sin_2omega;
    didt_rad_s = term1 + term2 + term3;

    % dOmega/dt
    term1 = -3/2 * n * J_2 * (RE_p)^2 * cos_i;
    term2 = +15/16 * n * J_4 * (RE_p)^4 * cos_i * (4 - 7* ...
        sin2_i) * (1 + 1.5*e^2);
    term3 = -3/16 * n * J_2^2 * (RE_p)^4 * cos_i * (7 - 15* ...
        sin2_i) * e^2 * cos_2omega;
    term4 = -3/2 * n * J_3 * (RE_p)^3 * cot(i) * e * sin_omega * (15/4* ...
        sin2_i - 1);
    term5 = -15/16 * n * J_4 * (RE_p)^4 * cos_i * (3 - 7* ...
        sin2_i) * e^2 * cos_2omega;
    bracket = (9/4 + 3/2 * sqrt_1_e2) - sin2_i * (5/2 + 9/4 * ...
        sqrt_1_e2) + 1/4 * (1 + 5/4 * sin2_i) * e^2;   
    term6 = -3/2 * n * J_2^2 * (RE_p)^4 * cos_i * bracket;
    dOmegadt_rad_s = term1 + term2 + term3 + term4 + term5 + term6;

    % domega/dt
    term1 = 3/4 * n * J_2 * (RE_p)^2 * (4 - 5*sin2_i);
    bracket1 = 12 - 103/4 * sin2_i + 215/16 * sin4_i + (7/4 - 9/8 * ...
        sin2_i - 45/32 * sin4_i) * e^2;
    term2 = 3/4 * n * J_2^2 * (RE_p)^4 * bracket1;
    bracket2 = 3/2 * (1 - 3/2 * sin2_i) * (4 - 5*sin2_i) * sqrt_1_e2;
    term3 = 3/4 * n * J_2^2 * (RE_p)^4 * bracket2;
    bracket3 = (16 - 62 * sin2_i + 49 * sin4_i) + 3/4 * (24 - 84 * ...
        sin2_i + 63 * sin4_i) * e^2;
    term4 = -15/32 * n * J_4 * (RE_p)^4 * bracket3;
    term5 = 3/64 * n * J_2^2 * (RE_p)^4 * cos_2omega * (-2 * (14 - 15* ...
        sin2_i));
    term6 = 3/64 * n * J_2^2 * (RE_p)^4 * cos_2omega * e^2 * (28 - 158* ...
        sin2_i + 135*sin4_i);
    
    if e > 1e-9 && abs(sin_i) > 1e-9
        term7 = 3/8 * n * J_3 * (RE_p)^3 * (sin_omega / (e * ...
            sin_i)) * (4 - 5*sin2_i) * (sin2_i - e^2 * cos(i)^2);
        term8 = 3/4 * n * J_3 * (RE_p)^3 * (sin_omega / (e * ...
            sin_i)) * sin2_i * e^2 * (13 - 15*sin2_i);
    else
        term7 = 0; term8 = 0;
    end
    
    term9 = -3/16 * n * J_4 * (RE_p)^4 * (3 * sin2_i * (6 - 7* ...
        sin2_i) + (-18 + 105*sin2_i - 189/2 * sin4_i) * e^2) * cos_2omega;      
    domegadt_rad_s = (term1 + term2 + term3 + term4 + term5 + ...
        term6 + term7 + term8 + term9);

    % dM/dt - n
    term1 = 3/2 * n * J_2 * (RE_p)^2 * (1 - 3/2 * sin2_i) * sqrt_1_e2;
    term2 = -15/8 * n * J_2^2 * (RE_p)^4 * sqrt(1-e) * (-1 + 5/2 * ...
        sin2_i - 13/8 * sin4_i);
    term3 = -15/16 * n * J_2^2 * (RE_p)^4 * sqrt(1-e) * e^2 * (-1 + ...
        sin2_i + 5/8 * sin4_i);
    term4 = -3/2 * n * J_2^2 * (RE_p)^4 * (1 - e^2) * (1 - 3/2 * sin2_i);
    term5 = -45/128 * n * J_4 * (RE_p)^4 * sqrt_1_e2 * e^2 * (8 - 40* ...
        sin2_i + 35*sin4_i);
    term6 = -9/64 * n * J_2^2 * (RE_p)^4 * sin2_i * sqrt_1_e2 * e^2 * (14 - ...
        15*sin2_i) * cos_2omega;
    term7 = 3/32 * n * J_2^2 * (RE_p)^4 * sin2_i * (1 - e^2)^(3/2) * (14 - ...
        15*sin2_i) * cos_2omega;
    term8 = -3/8 * n * J_3 * (RE_p)^3 * sin_i * sin_omega * sqrt_1_e2 * (4 - ...
        5*sin2_i) * (1 - 4*e^2)/e;
    if e < 1e-9, term8 = 0; end
    term9 = 15/64 * n * J_4 * (RE_p)^4 * sin2_i * cos_2omega * sqrt_1_e2 * (6 - ...
        7*sin2_i) * (2 - 5*e^2);
    term10 = 9/8 * n * J_2^2 * (RE_p)^4 * (1/sqrt_1_e2) * (3 - 15/2 * ...
        sin2_i + 47/8 * sin4_i);
    term11 = 9/8 * n * J_2^2 * (RE_p)^4 * (e^2/sqrt_1_e2) * (3/2 - 5* ...
        sin2_i + 117/16 * sin4_i);
    term12 = 9/64 * n * J_2^2 * (RE_p)^4 * (e^4/sqrt_1_e2) * (1 + 5* ...
        sin2_i - 101/8 * sin4_i);
    term13 = 81/1024 * n * J_2^2 * (RE_p)^4 * (e^4 * ...
        cos_4omega / sqrt_1_e2) * sin4_i;
    bracket_term14 = (70 - 123*sin2_i) * e^2 + 2 * (28 - 33*sin2_i) * e^4;
    term14 = 9/192 * n * J_2^2 * (RE_p)^4 * (cos_2omega / ...
        sqrt_1_e2) * sin2_i * bracket_term14;
    dMdt_rad_s = (term1 + term2 + term3 + term4 + term5 + term6 + term7 + ...
        term8 + term9 + term10 + term11 + term12 + term13 + term14);

    conv_factor = (180/pi) * (24 * 3600);
  
    dedt = dedt_rad_s * (24 * 3600); didt = didt_rad_s * conv_factor;
    dOmegadt = dOmegadt_rad_s * conv_factor; domegadt = domegadt_rad_s * conv_factor;
    dMdt_perturb = dMdt_rad_s * conv_factor;
end

function [Z_num, t_out, h_out, a_out, e_out] = ObliczZ5( ...
    a0, e0, i0, O0, w0, m, AD, ASRP, CD, CP, t0, rho_vec, C, wynik)
    
    ca = a0; ce = e0; ci = i0; cO = O0; cw = w0; ct = t0; elapsed_days = 0;
    
    max_iter = 1000000; min_h = C.hpp; max_lifetime_years = 60000; 
    
    hist_t = zeros(max_iter, 1); hist_hp = zeros(max_iter, 1); 
    hist_a = zeros(max_iter, 1); hist_e = zeros(max_iter, 1);
    rec_idx = 0; dt = 1.0;

    for iter = 1:max_iter
        hp = ca * (1 - ce) - C.R_E;
        
        rec_idx = rec_idx + 1; hist_t(rec_idx) = elapsed_days / 365.25; 
        hist_hp(rec_idx) = hp / 1000; hist_a(rec_idx) = ca / 1000;             
        hist_e(rec_idx) = ce;                    
        
        if hp < min_h, break; end
        if isnan(hp) || isinf(hp)
            fprintf('BŁĄD: Wykryto NaN w iteracji %d.\n', iter); break;
        end
        if elapsed_days > max_lifetime_years * 365.25
             Z_num = max_lifetime_years; break;
        end
        
        [res, ~] = ObliczPerturbacje(ca, ce, ci, cO, cw, m, AD, ASRP, CD, ...
            CP, ct, rho_vec, C);
        
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

    if wynik == 1
    if Z_num < 1
        fprintf('Metoda 4 (Numeryczna):      %.1f dni\n', elapsed_days);
    else
        fprintf('Metoda 4 (Numeryczna):      %.4f lat\n', Z_num);
    end
    end
    
    if rec_idx > 0
        t_out = hist_t(1:rec_idx); h_out = hist_hp(1:rec_idx);
        a_out = hist_a(1:rec_idx); e_out = hist_e(1:rec_idx);
    else
        t_out = []; h_out = []; a_out = []; e_out = [];        
    end
end

function PorownajModele(vec)
    
    theta = vec.theta_deg;
    
    c_Geo  = 'g'; s_Sim = '--'; s_Com = '-'; w_Lin = 1.5;
    
    figure('Name', 'Porównanie Modeli Geopotencjału', 'Units', ...
        'normalized', 'Position', [0.2 0.2 0.6 0.6]);
    
    de_S = vec.geo_simple.de; di_S = vec.geo_simple.di; dO_S = vec.geo_simple.dOMEGA;
    dw_S = vec.geo_simple.domega; dM_S = vec.geo_simple.dM;
    de_C = vec.geo_complex.de; di_C = vec.geo_complex.di; dO_C = vec.geo_complex.dOMEGA;
    dw_C = vec.geo_complex.domega; dM_C = vec.geo_complex.dM;
    
    % de/dt
    subplot(2,3,1); hold on; grid on;
    h1 = plot(theta, ones(size(theta)) * de_S, 'Color', ...
        c_Geo, 'LineStyle', s_Sim, 'LineWidth', w_Lin);
    h2 = plot(theta, ones(size(theta)) * de_C, 'Color', ...
        c_Geo, 'LineStyle', s_Com, 'LineWidth', w_Lin);
    title('de/dt [1/day]'); xlabel('\theta [deg]');
    
    % di/dt
    subplot(2,3,2); hold on; grid on;
    plot(theta, ones(size(theta)) * di_S, 'Color', ...
        c_Geo, 'LineStyle', s_Sim, 'LineWidth', w_Lin);
    plot(theta, ones(size(theta)) * di_C, 'Color', ...
        c_Geo, 'LineStyle', s_Com, 'LineWidth', w_Lin);
    title('di/dt [deg/day]'); xlabel('\theta [deg]');
    
    % dOmega/dt
    subplot(2,3,3); hold on; grid on;
    plot(theta, ones(size(theta)) * dO_S, 'Color', ...
        c_Geo, 'LineStyle', s_Sim, 'LineWidth', w_Lin);
    plot(theta, ones(size(theta)) * dO_C, 'Color', ...
        c_Geo, 'LineStyle', s_Com, 'LineWidth', w_Lin);
    title('d\Omega/dt [deg/day]'); xlabel('\theta [deg]');
    
    % domega/dt
    subplot(2,3,4); hold on; grid on;
    plot(theta, ones(size(theta)) * dw_S, 'Color', ...
        c_Geo, 'LineStyle', s_Sim, 'LineWidth', w_Lin);
    plot(theta, ones(size(theta)) * dw_C, 'Color', ...
        c_Geo, 'LineStyle', s_Com, 'LineWidth', w_Lin);
    title('d\omega/dt [deg/day]'); xlabel('\theta [deg]');
    
    % dM/dt
    subplot(2,3,5); hold on; grid on;
    plot(theta, ones(size(theta)) * dM_S, 'Color', c_Geo, ...
        'LineStyle', s_Sim, 'LineWidth', w_Lin);
    plot(theta, ones(size(theta)) * dM_C, 'Color', ...
        c_Geo, 'LineStyle', s_Com, 'LineWidth', w_Lin);
    title('dM/dt - n [deg/day]'); xlabel('\theta [deg]');
    
    sgtitle('Wpływ Geopotencjału: Model Uproszczony vs Dokładny');
    
    lgd = legend([h1, h2], 'Geo Simple (J2)', 'Geo Complex (J2^2, J3...)', ...
        'Orientation', 'horizontal');
    
    set(lgd, 'Units', 'normalized');
    set(lgd, 'Position', [0.25, 0.02, 0.5, 0.05]); 
    
end
