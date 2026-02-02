% Wyznaczanie N_{life} / C_B 
close all force;
clearvars;
clc;

%% ---------- USTAWIENIA ----------
filename = 'wykres_zywotnosc.txt';
deg1 = 3;   % stopień powierzchni dla e <= 0.1
deg2 = 3;   % stopień powierzchni dla e > 0.1
ridge_lambda = 1e-9;
big_weight = 1e8;
tol_e_zero = 1e-8;
n_h_samples = 200;   % ile próbek w h dla każdej krzywej ekstrapolowanej
h_min = 100; h_max = 700;
h_sample_vec = linspace(h_min, h_max, n_h_samples);

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
data = readtable(filename, 'Delimiter', '\t', 'VariableNamingRule', 'preserve');
e_vals = data.('e [-]');
h_vals = data.('h [km]');
NC_vals = data.('N/C_B');
unique_e = unique(e_vals);

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
        p = polyfit(h, lnN, 2);
    elseif numel(h) == 2
        p_lin = polyfit(h, lnN, 1);
        p = [0, p_lin(1), p_lin(2)];
    else
        p = [0,0,mean(lnN)];
    end
    coeffs_all = [coeffs_all; e_current, p];
end

%% ---------- Korekcja przesunięcia e=0.6 względem e=0.4 ----------
idx_e04 = find(abs(coeffs_all(:,1)-0.4)<1e-6,1);
idx_e06 = find(abs(coeffs_all(:,1)-0.6)<1e-6,1);
if ~isempty(idx_e04) && ~isempty(idx_e06)
    p_e04 = coeffs_all(idx_e04,2:4);
    lnN_e04 = polyval(p_e04, h_sample_vec);

    p_e06 = coeffs_all(idx_e06,2:4);
    lnN_e06 = polyval(p_e06, h_sample_vec);

    shift_ln = mean(lnN_e06 - lnN_e04) + 0.1;
    lnN_v2 = lnN_e04 + shift_ln;

    p_v2 = polyfit(h_sample_vec, lnN_v2, 2);
    coeffs_all(idx_e06,2:4) = p_v2;
end

%% ---------- Tworzenie próbek dla modeli powierzchni ----------
% Model 1: e <= 0.1
rows1 = find(coeffs_all(:,1) <= 0.1);
H_samps1 = []; E_samps1 = []; Log10N_samps1 = []; W_samps1 = [];

for k = 1:length(rows1)
    e_val = coeffs_all(rows1(k),1);
    p_ln = coeffs_all(rows1(k),2:4);
    lnN_samples = polyval(p_ln, h_sample_vec);
    log10N_samples = lnN_samples / log(10);

    H_samps1 = [H_samps1(:); h_sample_vec(:)];
    E_samps1 = [E_samps1(:); e_val*ones(size(h_sample_vec(:)))];
    Log10N_samps1 = [Log10N_samps1(:); log10N_samples(:)];

    w = ones(size(h_sample_vec(:)));
    if abs(e_val) <= tol_e_zero
        w = w * big_weight;
    end
    W_samps1 = [W_samps1(:); w(:)];
end

% Model 2: e > 0.1
rows2 = find(coeffs_all(:,1) > 0.1);
H_samps2 = []; E_samps2 = []; Log10N_samps2 = []; W_samps2 = [];

for k = 1:length(rows2)
    e_val = coeffs_all(rows2(k),1);
    p_ln = coeffs_all(rows2(k),2:4);
    lnN_samples = polyval(p_ln, h_sample_vec);
    log10N_samples = lnN_samples / log(10);

    H_samps2 = [H_samps2(:); h_sample_vec(:)];
    E_samps2 = [E_samps2(:); e_val*ones(size(h_sample_vec(:)))];
    Log10N_samps2 = [Log10N_samps2(:); log10N_samples(:)];
    W_samps2 = [W_samps2(:); ones(size(h_sample_vec(:)))];
end

%% ---------- Normalizacja i dopasowanie powierzchni ----------
% Model 1
mx1 = mean(H_samps1); sx1 = std(H_samps1); if sx1==0, sx1=1; end
my1 = mean(E_samps1); sy1 = std(E_samps1); if sy1==0, sy1=1; end

xn1 = (H_samps1 - mx1) / sx1;
yn1 = (E_samps1 - my1) / sy1;

A_surf1 = Wielomian(xn1, yn1, deg1);
Wdiag1 = spdiags(W_samps1,0,length(W_samps1),length(W_samps1));
M_surf1 = (A_surf1' * Wdiag1 * A_surf1) + ridge_lambda * eye(size(A_surf1,2));
b_surf1 = A_surf1' * Wdiag1 * Log10N_samps1;
coeff1 = M_surf1 \ b_surf1;

% Model 2
mx2 = mean(H_samps2); sx2 = std(H_samps2); if sx2==0, sx2=1; end
my2 = mean(E_samps2); sy2 = std(E_samps2); if sy2==0, sy2=1; end

xn2 = (H_samps2 - mx2) / sx2;
yn2 = (E_samps2 - my2) / sy2;

A_surf2 = Wielomian(xn2, yn2, deg2);
Wdiag2 = spdiags(W_samps2,0,length(W_samps2),length(W_samps2));
M_surf2 = (A_surf2' * Wdiag2 * A_surf2) + ridge_lambda * eye(size(A_surf2,2));
b_surf2 = A_surf2' * Wdiag2 * Log10N_samps2;
coeff2 = M_surf2 \ b_surf2;

%% ---------- FUNKCJA ObliczZ3 ----------
ObliczZ3_handle = @(e_query, h_query) ObliczZ3(e_query, h_query, ...
    coeff1, deg1, mx1, sx1, my1, sy1, ...
    coeff2, deg2, mx2, sx2, my2, sy2);

%% ----------------- FUNKCJE POMOCNICZE -----------------
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

function X = ObliczZ3(e_query, h_query, ...
    coeff1, deg1, mx1, sx1, my1, sy1, ...
    coeff2, deg2, mx2, sx2, my2, sy2)

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

save('coeffs_X.mat', 'coeff1', 'deg1', 'mx1', 'sx1', 'my1', 'sy1', ...
                     'coeff2', 'deg2', 'mx2', 'sx2', 'my2', 'sy2');


% Po uruchomieniu skryptu:
X_val = ObliczZ3_handle([0.5 0.4 0.2], [1000 1000 1000])