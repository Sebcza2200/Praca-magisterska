%% combined_surface_fit_to_lines.m
% Dopasowanie powierzchni 3. rzędu nie do punktów pomiarowych,
% a do linii ekstrapolowanych N(h) dla poszczególnych e.
close all force;
clearvars;
clc;

%% ---------- USTAWIENIA ----------
filename = 'wykres_zywotnosc.txt';
deg1 = 3;
deg2 = 3;
deg3 = 5;

nx = 140; ny = 140;
ridge_lambda = 1e-9;
big_weight = 1e8;
tol_e_zero = 1e-8;
n_h_samples = 200;   % ile próbek w h dla każdej krzywej ekstrapolowanej

%% ---------- WCZYTANIE DANYCH ----------
opts = detectImportOptions(filename, 'FileType', 'text');
opts.VariableNamingRule = 'preserve';
opts.Delimiter = '\t';
T = readtable(filename, opts);

x = T.("h [km]");
y = T.("e [-]");
z = T.("N/C_B");

ok = z > 0 & isfinite(z);
x = x(ok); y = y(ok); z = z(ok);

%% ---------- 1D DOPASOWANIA (krzywe ekstrapolacyjne) ----------
% dla każdej unikalnej e dopasuj wielomian 2 rzędu do ln(N)
data = readtable(filename, 'Delimiter', '\t', 'VariableNamingRule', 'preserve');
e_vals = data.('e [-]');
h_vals = data.('h [km]');
NC_vals = data.('N/C_B');
unique_e = unique(e_vals);

colors = lines(length(unique_e));
coeffs_all = [];

for i = 1:length(unique_e)
    e_current = unique_e(i);
    idx = e_vals == e_current;
    h = h_vals(idx);
    N = NC_vals(idx);
    valid = N > 0 & isfinite(h);
    h = h(valid); N = N(valid);
    if isempty(h)
        continue;
    end
    [h, ord] = sort(h); N = N(ord);
    lnN = log(N); % natural log
    if numel(h) >= 3
        p = polyfit(h, lnN, 2); % ln(N) = a*h^2 + b*h + c
    elseif numel(h) == 2
        p_lin = polyfit(h, lnN, 1);
        p = [0, p_lin(1), p_lin(2)];
    else
        p = [0,0,mean(lnN)];
    end
    coeffs_all = [coeffs_all; e_current, p];
end

% (opcjonalnie) rysowanie linii 1D - oryginalne krzywe
figure('Name','Krzywe ekstrapolowane','Color','w');
hold on; grid on; set(gca,'YScale','log');
xlabel('h_p [km]'); ylabel('N_{life}/C_B');
title('Ekstrapolowane krzywe (skala log)');

% --- naniesienie danych pomiarowych (wszystkich e) ---
scatter(h_vals, NC_vals, 25, 'ko', 'filled', 'DisplayName', 'Dane oryginalne', ...
    'MarkerFaceAlpha', 0.6, 'MarkerEdgeAlpha', 0.6);

% --- przygotowanie legendy i kolory
legendEntries = {'Dane oryginalne'};

for i = 1:size(coeffs_all,1)
    e_current = coeffs_all(i,1);
    p = coeffs_all(i,2:4);
    h_full = linspace(170,800,400);

    % jeśli e_current == 0.6, zamiast rysować faktyczne N_full,
    % rysujemy przesuniętą krzywą e=0.4 v2
    if abs(e_current - 0.6) < 1e-6
        % --- krzywa e=0.4
        idx_e04 = find(abs(coeffs_all(:,1)-0.4)<1e-6,1);
        p_e04 = coeffs_all(idx_e04,2:4);
        lnN_e04 = polyval(p_e04, h_full);
        % lnN e=0.6 faktyczne (dla przesunięcia)
        idx_e06 = find(abs(coeffs_all(:,1)-0.6)<1e-6,1);
        p_e06 = coeffs_all(idx_e06,2:4);
        lnN_e06 = polyval(p_e06, h_full);
        % przesunięcie w logspace
        shift_ln = mean(lnN_e06 - lnN_e04)+0.1;
        lnN_v2 = lnN_e04 + shift_ln;
        N_full = exp(lnN_v2); % krzywa „oszukana”
    else
        lnN = polyval(p, h_full);
        N_full = exp(lnN);
        N_full = cummax(N_full);
    end

    % --- obliczenie lokalnego R²
    mask_e = abs(e_vals - e_current) < 1e-6;
    if sum(mask_e) >= 3
        h_local = h_vals(mask_e);
        NC_local = NC_vals(mask_e);
        lnN_pred_local = polyval(p, h_local);
        N_pred_local = exp(lnN_pred_local);
        R2_local = 1 - sum((NC_local - N_pred_local).^2) / sum((NC_local - mean(NC_local)).^2);
    else
        R2_local = NaN;
    end

    % --- rysowanie linii
    plot(h_full, N_full, 'Color', colors(i,:), 'LineWidth', 1.4, ...
        'DisplayName', sprintf('e=%.2f (R²=%.3f)', e_current, R2_local));

    % --- dodanie wpisu do legendy
    legendEntries{end+1} = sprintf('e=%.2f (R²=%.3f)', e_current, R2_local);
end

% --- znajdź indeks e=0.6
idx_e06 = find(abs(coeffs_all(:,1)-0.6)<1e-6,1);

% --- współczynniki krzywej e=0.4
idx_e04 = find(abs(coeffs_all(:,1)-0.4)<1e-6,1);
p_e04 = coeffs_all(idx_e04,2:4);
lnN_e04 = polyval(p_e04, h_full);  % h_full musi być tym samym wektorem co w dalszej części

% --- obliczenie średniej różnicy w logspace dla przesunięcia
p_e06 = coeffs_all(idx_e06,2:4);
lnN_e06 = polyval(p_e06, h_full);
shift_ln = mean(lnN_e06 - lnN_e04)+0.1;

% --- nowa krzywa (oszukana)
lnN_v2 = lnN_e04 + shift_ln;

% --- dopasowanie wielomianu kwadratowego do lnN_v2 (żeby pasowało do formatu coeffs_all)
p_v2 = polyfit(h_full, lnN_v2, 2);

% --- podmień w coeffs_all
coeffs_all(idx_e06,2:4) = p_v2;

% --- legenda
legend(legendEntries, 'Location', 'bestoutside', 'Interpreter', 'none');


%% ========== NOWE: tworzenie próbek z linii i dopasowanie powierzchni ==========

% Ustal zakres h (dla próbkowania linii)
h_min = 100; h_max = 700;
h_sample_vec = linspace(h_min, h_max, n_h_samples);

% --- Model 1: użyj tylko krzywych z e <= 0.1 ---
idx1 = (y >= 0) & (y <= 0.1);
x1 = x(idx1); y1 = y(idx1); z1 = z(idx1);

rows1 = find(coeffs_all(:,1) <= 0.1);
if isempty(rows1)
    warning('Brak krzywych dla e<=0.1 — model 1 nie zostanie policzony.');
else
    H_samps1 = [];
    E_samps1 = [];
    Log10N_samps1 = [];
    W_samps1 = [];
    for k = 1:length(rows1)
        e_val = coeffs_all(rows1(k),1);
        p_ln = coeffs_all(rows1(k),2:4); % coefficients for ln(N)
        lnN_samples = polyval(p_ln, h_sample_vec); % natural log
        % convert to log10
        log10N_samples = lnN_samples / log(10);
        H_samps1 = [H_samps1(:); h_sample_vec(:)];
        E_samps1 = [E_samps1(:); e_val*ones(size(h_sample_vec(:)))];
        Log10N_samps1 = [Log10N_samps1(:); log10N_samples(:)];
        % weights: default 1, boost near e=0
        w = ones(size(h_sample_vec(:)));
        if abs(e_val) <= tol_e_zero
            w = w * big_weight;
        end
        W_samps1 = [W_samps1(:); w(:)];
    end

    idx_orig_zero = find(abs(y - 0) <= tol_e_zero);
    if ~isempty(idx_orig_zero)
        H_samps1 = [H_samps1; x(idx_orig_zero)];
        E_samps1 = [E_samps1; y(idx_orig_zero)];
        Log10N_samps1 = [Log10N_samps1; log10(z(idx_orig_zero))];
        W_samps1 = [W_samps1; big_weight*ones(length(idx_orig_zero),1)];
    end

    % Normalizacja (użyjemy tej samej normalizacji jak poprzednio jeśli x1,y1 istniały)
    % Jeśli x1/y1 puste, użyj normalizacji próbek:
    if exist('x1','var') && ~isempty(x1)
        mx1 = mean(x1); sx1 = std(x1); if sx1==0, sx1=1; end
        my1 = mean(y1); sy1 = std(y1); if sy1==0, sy1=1; end
    else
        mx1 = mean(H_samps1); sx1 = std(H_samps1); if sx1==0, sx1=1; end
        my1 = mean(E_samps1); sy1 = std(E_samps1); if sy1==0, sy1=1; end
    end

    xn1 = (H_samps1 - mx1) / sx1;
    yn1 = (E_samps1 - my1) / sy1;

    A_surf1 = Wielomian(xn1, yn1, deg1);
    Wdiag1 = spdiags(W_samps1,0,length(W_samps1),length(W_samps1));
    M_surf1 = (A_surf1' * Wdiag1 * A_surf1) + ridge_lambda * eye(size(A_surf1,2));
    b_surf1 = A_surf1' * Wdiag1 * Log10N_samps1;
    coeff1 = M_surf1 \ b_surf1; % nowe coeff dla powierzchni (log10 space)

    % R^2 (ważony) względem próbek linii
    log10_fit1 = A_surf1 * coeff1;
    R2_surf1 = 1 - sum(W_samps1 .* (Log10N_samps1 - log10_fit1).^2) / sum(W_samps1 .* (Log10N_samps1 - mean(Log10N_samps1)).^2);
    fprintf('\nMODEL 1: deg=%d, R² (względem lini) = %.8f\n', deg1, R2_surf1);
% --- R² względem oryginalnych danych pomiarowych (nie linii)
if exist('x1','var') && ~isempty(x1)
    log10N_pred1 = Wielomian((x1 - mx1)/sx1, (y1 - my1)/sy1, deg1) * coeff1;
    log10N_true1 = log10(z1);
    R2_data1 = 1 - sum((log10N_true1 - log10N_pred1).^2) / sum((log10N_true1 - mean(log10N_true1)).^2);
    fprintf('R² (względem danych) = %.8f\n', R2_data1);
else
    R2_data1 = NaN;
end

end

% --- Model 2: użyj krzywych z e >= 0.1 ---
idx2 = (y >= 0.1);
x2 = x(idx2); y2 = y(idx2); z2 = z(idx2);

rows2 = find(coeffs_all(:,1) >= 0.1);
if isempty(rows2)
    warning('Brak krzywych dla e>=0.1 — model 2 nie zostanie policzony.');
else
    H_samps2 = [];
    E_samps2 = [];
    Log10N_samps2 = [];
    W_samps2 = [];
    for k = 1:length(rows2)
        e_val = coeffs_all(rows2(k),1);
        p_ln = coeffs_all(rows2(k),2:4);
        lnN_samples = polyval(p_ln, h_sample_vec);
        log10N_samples = lnN_samples / log(10);
        H_samps2 = [H_samps2(:); h_sample_vec(:)];
        E_samps2 = [E_samps2(:); e_val*ones(size(h_sample_vec(:)))];
        Log10N_samps2 = [Log10N_samps2(:); log10N_samples(:)];
        w = ones(size(h_sample_vec(:)));
        % przydatne: jeśli chcesz preferować zbieżność w okolicach e=0 również w modelu2,
        % możesz dać lekkie podniesienie wagi pobliskim e (ale tu e>=0.1 więc zwykle brak)
        W_samps2 = [W_samps2(:); w(:)];
    end

    % Normalizacja dla modelu 2 (użyjemy x2,y2 jeśli istnieją)
    if exist('x2','var') && ~isempty(x2)
        mx2 = mean(x2); sx2 = std(x2); if sx2==0, sx2=1; end
        my2 = mean(y2); sy2 = std(y2); if sy2==0, sy2=1; end
    else
        mx2 = mean(H_samps2); sx2 = std(H_samps2); if sx2==0, sx2=1; end
        my2 = mean(E_samps2); sy2 = std(E_samps2); if sy2==0, sy2=1; end
    end

    xn2 = (H_samps2 - mx2) / sx2;
    yn2 = (E_samps2 - my2) / sy2;

    A_surf2 = Wielomian(xn2, yn2, deg2);
    Wdiag2 = spdiags(W_samps2,0,length(W_samps2),length(W_samps2));
    M_surf2 = (A_surf2' * Wdiag2 * A_surf2) + ridge_lambda * eye(size(A_surf2,2));
    b_surf2 = A_surf2' * Wdiag2 * Log10N_samps2;
    coeff2 = M_surf2 \ b_surf2; % nowe coeff dla powierzchni (log10 space)

    % R^2 (ważony)
    log10_fit2 = A_surf2 * coeff2;
    R2_surf2 = 1 - sum(W_samps2 .* (Log10N_samps2 - log10_fit2).^2) / sum(W_samps2 .* (Log10N_samps2 - mean(Log10N_samps2)).^2);
    fprintf('MODEL 2: deg=%d, R² (względem lini) = %.8f\n', deg2, R2_surf2);
% --- R² względem oryginalnych danych pomiarowych (nie linii)
if exist('x2','var') && ~isempty(x2)
    log10N_pred2 = Wielomian((x2 - mx2)/sx2, (y2 - my2)/sy2, deg2) * coeff2;
    log10N_true2 = log10(z2);
    R2_data2 = 1 - sum((log10N_true2 - log10N_pred2).^2) / sum((log10N_true2 - mean(log10N_true2)).^2);
    fprintf('R² (względem danych) = %.8f\n', R2_data2);
else
    R2_data2 = NaN;
end

end

%% === MODEL 3: e ∈ [0, 0.6] ===
idx3 = (y >= 0);
x3 = x(idx3); y3 = y(idx3); z3 = z(idx3);

% --- filtrowanie danych ---
mask3 = (e_vals >= 0) & (e_vals <= 0.6);
h3 = h_vals(mask3);
e3 = e_vals(mask3);
NC3 = NC_vals(mask3);

% --- przekształcenie logarytmiczne ---
logNC3 = log(NC3);

% --- dopasowanie powierzchni deg3-rzędu z centrowaniem i skalowaniem ---
fitType3 = sprintf('poly%d%d', deg3, deg3);                % dynamiczny typ
fitoptions3 = fitoptions(fitType3, 'Normalize', 'on');     % automatyczne centrowanie i skalowanie
[fitresult3, gof3] = fit([h3, e3], logNC3, fitType3, fitoptions3);

% --- przewidywane wartości ---
logNC3_pred = feval(fitresult3, h3, e3);
NC3_pred = exp(logNC3_pred);

% --- R² względem danych ---
R2_data3 = 1 - sum((NC3 - NC3_pred).^2) / sum((NC3 - mean(NC3)).^2);

% --- R² względem powierzchni ---
[hq3, eq3] = meshgrid(linspace(min(h3), max(h3), 50), linspace(min(e3), max(e3), 50));
logNCq3 = feval(fitresult3, hq3, eq3);
NCq3 = exp(logNCq3);
logNCfit3 = log(NC3_pred);
R2_surf3 = 1 - sum((logNC3 - logNCfit3).^2) / sum((logNC3 - mean(logNC3)).^2);

% --- wypisanie R² do command window ---
fprintf('MODEL 3: deg=%d, R² (względem lini) = %.8f\n', deg3, R2_surf3);
fprintf('R² (względem danych) = %.8f\n', R2_data3);



%% ============================================================
% ZAPIS RÓWNAŃ MODELI I LINII EKSTRAPOLACYJNYCH DO PLIKÓW
% ============================================================

output_folder = fullfile(pwd, 'wyniki_równania');
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

%% --- Zapis modeli powierzchniowych ---
fid = fopen(fullfile(output_folder, 'modele_powierzchni.txt'), 'w');

fprintf(fid, '============================================================\n');
fprintf(fid, 'RÓWNANIA MODELI POWIERZCHNIOWYCH\n');
fprintf(fid, '============================================================\n\n');

% Model 1
if exist('coeff1','var')
    fprintf(fid, 'MODEL 1: e ∈ [0, 0.1]\n');
    fprintf(fid, 'Normalizacja: h_n = (h_p - %.6f)/%.6f,  e_n = (e - %.6f)/%.6f\n', mx1, sx1, my1, sy1);
    fprintf(fid, 'Stopień wielomianu: %d\n', deg1);
    fprintf(fid, '\nWzór:\n');
    fprintf(fid, 'N_{life} = 10^(Σ_{i+j≤3} c_{ij} * h_n^i * e_n^j)\n\n');
    fprintf(fid, 'Współczynniki c_{ij} (porządek: d=0..3, i=0..d, j=d-i):\n');

    count = 0;
    for d = 0:deg1
        for i = 0:d
            j = d - i;
            count = count + 1;
            fprintf(fid, 'c(%d,%d) = %.12e\n', i, j, coeff1(count));
        end
    end
    fprintf(fid, '\nR² (względem linii) = %.8f\n', exist('R2_surf1','var')*R2_surf1);
    fprintf(fid, '\nR² (względem danych) = %.8f\n', exist('R2_data1','var')*R2_data1);
    fprintf(fid, '\n------------------------------------------------------------\n\n');
end

% Model 2
if exist('coeff2','var')
    fprintf(fid, 'MODEL 2: e ≥ 0.1\n');
    fprintf(fid, 'Normalizacja: h_n = (h_p - %.6f)/%.6f,  e_n = (e - %.6f)/%.6f\n', mx2, sx2, my2, sy2);
    fprintf(fid, 'Stopień wielomianu: %d\n', deg2);
    fprintf(fid, '\nWzór:\n');
    fprintf(fid, 'N_{life} = 10^(Σ_{i+j≤3} c_{ij} * h_n^i * e_n^j)\n\n');
    fprintf(fid, 'Współczynniki c_{ij} (porządek: d=0..3, i=0..d, j=d-i):\n');

    count = 0;
    for d = 0:deg2
        for i = 0:d
            j = d - i;
            count = count + 1;
            fprintf(fid, 'c(%d,%d) = %.10e\n', i, j, coeff2(count));
        end
    end
    fprintf(fid, '\nR² (względem linii) = %.8f\n', exist('R2_surf2','var')*R2_surf2);
    fprintf(fid, '\nR² (względem danych) = %.8f\n', exist('R2_data2','var')*R2_data2);
    fprintf(fid, '\n------------------------------------------------------------\n\n');
end

% --- Model 3 ---
fprintf(fid, 'MODEL 3: pełny zakres e\n');
eqn3 = formula(fitresult3);
coeffs3 = coeffvalues(fitresult3);
coeffnames3 = coeffnames(fitresult3);
eqn_text3 = char(eqn3);
for i = 1:length(coeffs3)
    eqn_text3 = strrep(eqn_text3, coeffnames3{i}, sprintf('%.12e', coeffs3(i)));
end
eqn_text3 = sprintf('z = exp(%s)\n', eqn_text3);
fprintf(fid, '%s\n', eqn_text3);

    fprintf(fid, '\nR² (względem linii) = %.8f\n', exist('R2_surf3','var')*R2_surf3);
    fprintf(fid, '\nR² (względem danych) = %.8f\n', exist('R2_data3','var')*R2_data3);
    fprintf(fid, '\n------------------------------------------------------------\n\n');

fclose(fid);

%% --- Zapis ekstrapolowanych krzywych N(h) (oryginał) ---
fid = fopen(fullfile(output_folder, 'linie_ekstrapolacyjne.txt'), 'w');
fprintf(fid, 'Współczynniki do równań ekstrapolacji danych, postaci:\n [ln(N) = a*h^2 + b*h + c]\n');
fprintf(fid, 'e\t a (h^2)\t b (h)\t c (const) \n');
for i = 1:size(coeffs_all,1)
    fprintf(fid, '%.3f\t %.12e\t %.12e\t %.12e\n', coeffs_all(i,1), coeffs_all(i,2), coeffs_all(i,3), coeffs_all(i,4));
end
fclose(fid);

%% ---------- POWIERZCHNIE 3D + Naniesienie krzywych i danych ----------
% Użyjemy teraz coeff1/coeff2 (powierzchnie dopasowane do linii)
h_min_plot = 100; h_max_plot = 700;
hq1_plot = linspace(h_min_plot, h_max_plot, nx);
if exist('y1','var') && ~isempty(y1)
    eq1_plot = linspace(min(y1), max(y1), ny);
else
    eq1_plot = linspace(0,0.1,ny);
end
[Xq1_plot, Yq1_plot] = meshgrid(hq1_plot, eq1_plot);
Aq1_plot = Wielomian((Xq1_plot - mx1)/sx1, (Yq1_plot - my1)/sy1, deg1);
Zq1_plot = reshape(10.^(Aq1_plot * coeff1), size(Xq1_plot));

hq2_plot = linspace(h_min_plot, h_max_plot, nx);
if exist('y2','var') && ~isempty(y2)
    eq2_plot = linspace(min(y2), max(y2), ny);
else
    eq2_plot = linspace(0.1, max(unique_e), ny);
end
[Xq2_plot, Yq2_plot] = meshgrid(hq2_plot, eq2_plot);
Aq2_plot = Wielomian((Xq2_plot - mx2)/sx2, (Yq2_plot - my2)/sy2, deg2);
Zq2_plot = reshape(10.^(Aq2_plot * coeff2), size(Xq2_plot));

figure('Name','Modele powierzchniowe','Position',[100 100 1200 500]);

% MODEL 1 subplot
subplot(1,2,1);
colormap(turbo);
s1 = surf(Xq1_plot, Yq1_plot, Zq1_plot, log10(Zq1_plot), 'EdgeColor', 'none', 'DisplayName', 'Model interpolacji powierzchniowej');
set(gca,'ZScale','log'); view(45,30); grid on; hold on;
cb = colorbar;
cb.Label.String = 'log_{10}(N_{life}/C_B)';
cb.Label.FontSize = 12; cb.Label.FontWeight = 'bold';
xlabel('h_p [km]'); ylabel('e [-]'); zlabel('N_{life}/C_B');
title(sprintf('MODEL 1: e∈[0, 0.1], R² (wzgl. linii)=%.4f, (wzgl. punktów)=%.4f', ...
    exist('R2_surf1','var')*R2_surf1, exist('R2_data1','var')*R2_data1));

% nałóż krzywe (ekstrapolacje 1D) i dane pomiarowe
for i = 1:size(coeffs_all,1)
    e_current = coeffs_all(i,1);
    if e_current <= 0.1
        p_ln = coeffs_all(i,2:4);
        lnN_vals = polyval(p_ln, hq1_plot);
        N_vals = exp(lnN_vals);
        N_vals = cummax(N_vals);
        plot3(hq1_plot, e_current*ones(size(hq1_plot)), N_vals, ...
    'Color', [0.3 0.3 0.3], 'LineWidth', 2.2, ...
    'DisplayName', 'Model ekstrapolacji liniowej');
                % --- naniesienie danych pomiarowych ---
scatter3(x1, y1, z1, 25, 'ko', 'filled', 'DisplayName', 'Dane oryginalne');
    end
end

legend({'Model interpolacji powierzchniowej', 'Model ekstrapolacji liniowej', 'Dane oryginalne'}, 'Location','best');

% MODEL 2 subplot
subplot(1,2,2);
colormap(turbo);
s2 = surf(Xq2_plot, Yq2_plot, Zq2_plot, log10(Zq2_plot), 'EdgeColor', 'none', 'DisplayName', 'Model interpolacji powierzchniowej');
set(gca,'ZScale','log'); view(45,30); grid on; hold on;
cb = colorbar;
cb.Label.String = 'log_{10}(N_{life}/C_B)';
cb.Label.FontSize = 12; cb.Label.FontWeight = 'bold';
xlabel('h_p [km]'); ylabel('e [-]'); zlabel('N_{life}/C_B');
title(sprintf('MODEL 2: e∈(0.1, 0.6], R² (wzgl. linii)=%.4f, (wzgl. punktów)=%.4f', ...
    exist('R2_surf2','var')*R2_surf2, exist('R2_data2','var')*R2_data2));

for i = 1:size(coeffs_all,1)
    e_current = coeffs_all(i,1);
    if e_current >= 0.1
        p_ln = coeffs_all(i,2:4);
        lnN_vals = polyval(p_ln, hq2_plot);
        N_vals = exp(lnN_vals);
        N_vals = cummax(N_vals);
plot3(hq2_plot, e_current*ones(size(hq2_plot)), N_vals, ...
    'Color', [0.3 0.3 0.3], 'LineWidth', 2.2, ...
    'DisplayName', 'Model ekstrapolacji liniowej');
            % --- naniesienie danych pomiarowych ---
scatter3(x2, y2, z2, 25, 'ko', 'filled', 'DisplayName', 'Dane oryginalne');
    end
end


legend({'Model interpolacji powierzchniowej', 'Model ekstrapolacji liniowej', 'Dane oryginalne'}, 'Location','best');


%% ---------- MODEL 3 ----------
hq3_plot = linspace(h_min_plot, h_max_plot, nx);
if exist('e3','var') && ~isempty(e3)
    eq3_plot = linspace(min(e3), max(e3), ny);
else
    eq3_plot = linspace(0, 0.6, ny);
end
[Xq3_plot, Yq3_plot] = meshgrid(hq3_plot, eq3_plot);

% Oblicz Z (fitresult3 przewiduje ln(N))
logNCq3 = feval(fitresult3, Xq3_plot, Yq3_plot);
Zq3_plot = reshape(exp(logNCq3), size(Xq3_plot));

figure('Name','MODEL 3: e∈[0,0.6]','Position',[100 100 700 600]);
colormap(turbo);

% --- rysowanie powierzchni (z kolorem log10) ---
s3 = surf(Xq3_plot, Yq3_plot, Zq3_plot, log10(Zq3_plot), ...
    'EdgeColor', 'none', 'DisplayName', 'Model interpolacji powierzchniowej');
set(gca,'ZScale','log'); view(45,30); grid on; hold on;

cb = colorbar;
cb.Label.String = 'log_{10}(N_{life}/C_B)';
cb.Label.FontSize = 12;
cb.Label.FontWeight = 'bold';

xlabel('h_p [km]');
ylabel('e [-]');
zlabel('N_{life}/C_B');
title(sprintf('MODEL 3: e∈[0, 0.6], R²(wzgl. linii)=%.4f, R²(wzgl. punktów)=%.4f', ...
    R2_surf3, R2_data3));

% --- placeholdery do legendy ---
h_line = plot3(nan, nan, nan, 'Color', [0.3 0.3 0.3], 'LineWidth', 2.2, ...
    'DisplayName', 'Model ekstrapolacji liniowej');
h_data = scatter3(nan, nan, nan, 25, 'ko', 'filled', ...
    'DisplayName', 'Dane oryginalne');

% --- ekstrapolowane linie 1D ---
for i = 1:size(coeffs_all,1)
    e_current = coeffs_all(i,1);
    if e_current >= 0 && e_current <= 0.6
        p_ln = coeffs_all(i,2:4);
        lnN_vals = polyval(p_ln, hq3_plot);
        N_vals = exp(lnN_vals);
        N_vals = cummax(N_vals);
        plot3(hq3_plot, e_current*ones(size(hq3_plot)), N_vals, ...
            'Color', [0.3 0.3 0.3], 'LineWidth', 2.2, 'HandleVisibility','off');
    end
end

% --- naniesienie danych pomiarowych (na końcu, żeby były na wierzchu) ---
if exist('h3','var') && ~isempty(h3)
    scatter3(x3, y3, z3, 25, 'ko', 'filled', 'HandleVisibility','off');
end

legend([s3, h_line, h_data], ...
    {'Model interpolacji powierzchniowej', 'Model ekstrapolacji liniowej', 'Dane oryginalne'}, ...
    'Location', 'best', 'Interpreter', 'none');

function A = Wielomian(x, y, deg)
% Buduje macierz wielomianów Σ_{i+j<=deg} x^i y^j
count = 0;
for d = 0:deg
    for i = 0:d
        j = d - i;
        count = count + 1;
        A(:,count) = (x(:).^i) .* (y(:).^j);
    end
end
end

