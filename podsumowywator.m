clc; clear; close all;

%% 1. Konfiguracja plików i kolorów
filenames = {'ONEWEB.txt', 'STARLINK.txt', 'KUIPER.txt', 'GUOWANG.txt', 'QIANFAN.txt'};
colors = [
    0.00, 0.45, 0.74; % Niebieski
    0.85, 0.33, 0.10; % Pomarańczowy
    0.93, 0.69, 0.13; % Żółty
    0.49, 0.18, 0.56; % Fioletowy
    0.47, 0.67, 0.19  % Zielony
];

all_data_raw = struct();

%% 2. Pancerny import danych (Obsługa 5 kolumn numerycznych)
fprintf('Wczytywanie danych...\n');
valid_idx = 0;
for i = 1:length(filenames)
    fname = filenames{i};
    if isfile(fname)
        fid = fopen(fname, 'r');
        fgetl(fid); % Pominięcie nagłówka
        
        hp_tmp = []; z1_tmp = []; z2_tmp = []; z5_tmp = []; e_tmp = [];
        
        while ~feof(fid)
            line = fgetl(fid);
            if ischar(line) && ~isempty(strtrim(line))
                parts = strsplit(strtrim(line));
                if length(parts) >= 5
                    % Pobieramy 5 ostatnich wartości numerycznych
                    vals = str2double(parts(end-4:end));
                    
                    if all(~isnan(vals))
                        % MAPOWANIE KOLUMN (Jeśli Z5 dalej jest 0, zamień vals(4) z vals(5))
                        hp_tmp(end+1) = vals(1); 
                        z1_tmp(end+1) = vals(2);
                        z2_tmp(end+1) = vals(3);
                        z5_tmp(end+1) = vals(4);
                        e_tmp(end+1)  = vals(5); % e jako ostatnia kolumna
                    end
                end
            end
        end
        fclose(fid);
        
        if ~isempty(hp_tmp)
            valid_idx = valid_idx + 1;
            all_data_raw(valid_idx).name = strrep(fname, '.txt', '');
            all_data_raw(valid_idx).hp = hp_tmp';
            all_data_raw(valid_idx).z1 = z1_tmp';
            all_data_raw(valid_idx).z2 = z2_tmp';
            all_data_raw(valid_idx).z5 = z5_tmp';
            all_data_raw(valid_idx).color = colors(mod(i-1, size(colors,1))+1, :);
        end
    end
end

%% 3. Generowanie Wykresów (Dwa okna: Wszystko oraz < 750km)
for case_idx = 1:2
    all_data = all_data_raw;
    if case_idx == 1
        filter_label = 'Wszystkie dane';
        max_hp_limit = inf;
    else
        max_hp_limit = 750;
        filter_label = 'Filtr: hp <= 750 km';
        % Aplikacja filtra wysokości
        for i = 1:length(all_data)
            m = all_data(i).hp <= max_hp_limit;
            all_data(i).hp = all_data(i).hp(m);
            all_data(i).z1 = all_data(i).z1(m);
            all_data(i).z2 = all_data(i).z2(m);
            all_data(i).z5 = all_data(i).z5(m);
        end
    end
    
    % Zliczanie satelitów
    total_sat_count = 0;
    all_Y = []; all_HP = [];
    for i = 1:length(all_data)
        total_sat_count = total_sat_count + length(all_data(i).hp);
        all_Y = [all_Y; all_data(i).z1; all_data(i).z2; all_data(i).z5];
        all_HP = [all_HP; all_data(i).hp];
    end
    
    if total_sat_count == 0, continue; end
    
    % Ograniczenia osi i histogramu
    common_ylim = [0, max(all_Y) * 1.1];
    bin_width = 100;
    edges = floor(min(all_HP)/bin_width)*bin_width : bin_width : ceil(max(all_HP)/bin_width)*bin_width;

    figure('Name', [filter_label ' (N=' num2str(total_sat_count) ')'], 'Units', 'normalized', ...
           'Position', [0.05 + (case_idx-1)*0.05, 0.1, 0.85, 0.75], 'Color', 'w');
    tlo = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'loose');
    
    titles = {'Metoda 1 (King-Hele)', 'Metoda 2 (Uproszczona)', 'Metoda 3 (Numeryczna)'};
    fields = {'z1', 'z2', 'z5'};
    h_scatters = []; h_trend = [];

    for m = 1:3
        ax = nexttile; hold on; grid on;
        X_global = []; Y_global = [];
        
        % --- PRAWA OŚ (Histogram % co 100km) ---
        yyaxis right
        h_hist = histogram(all_HP, edges, 'Normalization', 'probability', ...
            'FaceColor', [0.8 0.8 0.8], 'FaceAlpha', 0.4, 'EdgeColor', [0.6 0.6 0.6]);
        ylabel('Udział populacji [%]');
        ax.YAxis(2).TickLabelFormat = '%.0f%%'; 
        curr_ticks = get(ax, 'YTick');
        set(ax, 'YTickLabel', string(curr_ticks * 100) + "%");
        ax.YAxis(2).Color = [0.4 0.4 0.4];
        
        % --- LEWA OŚ (Żywotność) ---
        yyaxis left
        for i = 1:length(all_data)
            if isempty(all_data(i).hp), continue; end
            h = scatter(all_data(i).hp, all_data(i).(fields{m}), 22, all_data(i).color, 'filled', 'MarkerFaceAlpha', 0.5);
            if m == 1, h_scatters(end+1) = h; end
            X_global = [X_global; all_data(i).hp];
            Y_global = [Y_global; all_data(i).(fields{m})];
        end
        
        % Linia trendu (Cienka, przerywana)
        v = Y_global > 0 & ~isnan(X_global);
        if any(v)
            p = polyfit(X_global(v), log(Y_global(v)), 1);
            x_range = linspace(min(X_global), max(X_global), 100);
            y_trend_val = exp(p(2)) * exp(p(1) * x_range);
            h_trend = plot(x_range, y_trend_val, 'k--', 'LineWidth', 0.8);
        end
        
        ylim(common_ylim);
        xlabel('Wysokość perygeum [km]');
        if m == 1, ylabel('Żywotność [Lata]'); end
        title(titles{m}, 'FontSize', 12, 'FontWeight', 'bold');
        ax.YAxis(1).Color = 'k';
    end

    % --- WSPÓLNA LEGENDA NA DOLE ---
    leg_entries = [h_scatters];
    leg_names = {all_data(1:length(h_scatters)).name};
    if ~isempty(h_trend)
        leg_entries = [leg_entries, h_trend];
        leg_names = [leg_names, 'Linia trendu'];
    end
    
    lgd = legend(leg_entries, leg_names);
    lgd.Layout.Tile = 'south';
    lgd.Orientation = 'horizontal';
    lgd.Interpreter = 'none';
    lgd.FontSize = 11;

    sgtitle(sprintf('Analiza żywotności LEO - %s | Populacja: %d satelitów', ...
            filter_label, total_sat_count), 'FontSize', 16, 'FontWeight', 'bold');
end