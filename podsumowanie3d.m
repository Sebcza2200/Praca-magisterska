clc; clear; close all;

%% 1. Konfiguracja plików i kolorów
filenames = {'ONEWEB.txt', 'STARLINK.txt', 'KUIPER.txt', 'GUOWANG.txt', 'QIANFAN.txt'};
colors = lines(length(filenames)); 

all_data_raw = struct();

%% 2. Pancerny import danych (5 kolumn: hp, Z1, Z2, Z5, e)
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
                    vals = str2double(parts(end-4:end));
                    if all(~isnan(vals))
                        % Mapowanie: hp, Z1, Z2, Z5, e
                        hp_tmp(end+1) = vals(1); 
                        z1_tmp(end+1) = vals(2);
                        z2_tmp(end+1) = vals(3);
                        z5_tmp(end+1) = vals(4);
                        e_tmp(end+1)  = vals(5); 
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
            all_data_raw(valid_idx).e  = e_tmp';
            all_data_raw(valid_idx).color = colors(valid_idx, :);
        end
    end
end

%% 3. GENEROWANIE WYKRESÓW 2D (Dwa okna: Pełne i Filtrowane)
for case_idx = 1:2
    all_data = all_data_raw;
    filter_txt = 'Wszystkie dane';
    if case_idx == 2
        max_hp = 750;
        filter_txt = ['hp <= ' num2str(max_hp) ' km'];
        for i = 1:length(all_data)
            m = all_data(i).hp <= max_hp;
            all_data(i).hp = all_data(i).hp(m);
            all_data(i).z1 = all_data(i).z1(m);
            all_data(i).z2 = all_data(i).z2(m);
            all_data(i).z5 = all_data(i).z5(m);
            all_data(i).e  = all_data(i).e(m);
        end
    end
    
    total_sats = sum(arrayfun(@(x) length(x.hp), all_data));
    if total_sats == 0, continue; end
    
    all_Z_case = [vertcat(all_data.z1); vertcat(all_data.z2); vertcat(all_data.z5)];
    all_HP_case = vertcat(all_data.hp);
    common_ylim = [0, max(all_Z_case) * 1.1];
    edges = floor(min(all_HP_case)/100)*100 : 100 : ceil(max(all_HP_case)/100)*100;

    figure('Name', ['2D: ' filter_txt], 'Units', 'normalized', 'Position', [0.05, 0.1, 0.8, 0.7], 'Color', 'w');
    tlo = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'loose');
    fields = {'z1', 'z2', 'z5'}; titles = {'Z1 (King-Hele)', 'Z2 (Uproszczona)', 'Z5 (Numeryczna)'};
    h_scat = [];

    for m = 1:3
        ax = nexttile; hold on; grid on;
        yyaxis right; histogram(all_HP_case, edges, 'Normalization', 'probability', 'FaceColor', [0.8 0.8 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        ylabel('Udział'); ax.YAxis(2).TickLabelFormat = '%.0f%%'; 
        curr = get(ax, 'YTick'); set(ax, 'YTickLabel', string(curr * 100) + "%");
        yyaxis left;
        X_curr = []; Y_curr = [];
        for i = 1:length(all_data)
            if isempty(all_data(i).hp), continue; end
            h = scatter(all_data(i).hp, all_data(i).(fields{m}), 20, all_data(i).color, 'filled', 'MarkerFaceAlpha', 0.5);
            if m == 1, h_scat(end+1) = h; end
            X_curr = [X_curr; all_data(i).hp]; Y_curr = [Y_curr; all_data(i).(fields{m})];
        end
        v = Y_curr > 0; if sum(v) > 2
            p = polyfit(X_curr(v), log(Y_curr(v)), 1); xr = linspace(min(X_curr), max(X_curr), 100);
            h_tr = plot(xr, exp(p(2))*exp(p(1)*xr), 'k--', 'LineWidth', 0.8);
        end
        ylim(common_ylim); title(titles{m});
    end
    lg = legend([h_scat, h_tr], [{all_data.name}, 'Linia trendu']); lg.Layout.Tile = 'south'; lg.Orientation = 'horizontal';
    sgtitle(['Analiza 2D | ' filter_txt ' | Populacja: ' num2str(total_sats)], 'FontWeight', 'bold');
end

%% 4. ANALIZY 3D (Dwa okna: Wszystko oraz < 750km)
for d3_case = 1:2
    all_data_3d = all_data_raw;
    filter_label = 'Wszystkie dane';
    if d3_case == 2
        max_hp = 750;
        filter_label = ['hp <= ' num2str(max_hp) ' km'];
        for i = 1:length(all_data_3d)
            m = all_data_3d(i).hp <= max_hp;
            all_data_3d(i).hp = all_data_3d(i).hp(m);
            all_data_3d(i).z1 = all_data_3d(i).z1(m);
            all_data_3d(i).z2 = all_data_3d(i).z2(m);
            all_data_3d(i).z5 = all_data_3d(i).z5(m);
            all_data_3d(i).e  = all_data_3d(i).e(m);
        end
    end
    
    total_3d_sats = sum(arrayfun(@(x) length(x.hp), all_data_3d));
    if total_3d_sats == 0, continue; end
    
    % Wyznaczanie limitu Z dla bieżącego okna 3D
    all_Z_3d = [vertcat(all_data_3d.z1); vertcat(all_data_3d.z2); vertcat(all_data_3d.z5)];
    z_max_3d = max(all_Z_3d) * 1.05;

    figure('Name', ['3D: ' filter_label], 'Units', 'normalized', 'Position', [0.1, 0.1, 0.85, 0.75], 'Color', 'w');
    tlo3d = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'loose');
    z_fields = {'z1', 'z2', 'z5'}; z_titles = {'Żywotność Z1', 'Żywotność Z2', 'Żywotność Z5'};
    h_scat3d = [];

    for m = 1:3
        ax3d = nexttile; hold on; grid on;
        for i = 1:length(all_data_3d)
            if isempty(all_data_3d(i).hp), continue; end
            h = scatter3(all_data_3d(i).hp, all_data_3d(i).e, all_data_3d(i).(z_fields{m}), ...
                         25, all_data_3d(i).color, 'filled', 'MarkerFaceAlpha', 0.6);
            if m == 1, h_scat3d(end+1) = h; end
        end
        xlabel('hp [km]'); ylabel('e [-]'); zlabel('Żywotność [Lata]');
        title(z_titles{m}); zlim([0, z_max_3d]); view(45, 30);
    end

    lg3d = legend(h_scat3d, {all_data_3d(1:length(h_scat3d)).name});
    lg3d.Layout.Tile = 'south'; lg3d.Orientation = 'horizontal'; lg3d.Interpreter = 'none';
    sgtitle(['Analiza 3D | ' filter_label ' | Populacja: ' num2str(total_3d_sats)], 'FontSize', 16, 'FontWeight', 'bold');
end