function GenerateSpaceDatabase()
    % GenerateSpaceDatabase - Generuje plik Excel ze wszystkimi satelitami
    % WERSJA LEO ONLY + PODZIAŁ NA 4 ARKUSZE (Starlink/OneWeb/Konstelacje/Inne)
    
    clc;
    disp('==========================================================');
    disp('GENERATOR BAZY DANYCH OBIEKTÓW KOSMICZNYCH (OFFLINE MODE)');
    disp('==========================================================');

    % Opcje połączenia (udajemy przeglądarkę + długi timeout)
    netOpts = weboptions('Timeout', 60, 'UserAgent', 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36');

    %% KROK 1: Pobranie Masowego Katalogu Fizycznego (GCAT)
    url_gcat = 'https://planet4589.org/space/gcat/tsv/cat/satcat.tsv';
    local_gcat_file = 'satcat_physical.tsv'; % Stała nazwa
    
    disp('1. Sprawdzanie katalogu fizycznego GCAT...');
    
    file_exists = false;
    if exist(local_gcat_file, 'file')
        disp('   Plik satcat_physical.tsv już istnieje. Pomijam pobieranie.');
        file_exists = true;
    else
        disp('   Pobieranie katalogu (ok. 15-20 MB)...');
        try
            websave(local_gcat_file, url_gcat, netOpts);
            fprintf('   Pobrano plik: %s\n', local_gcat_file);
            file_exists = true;
        catch ME
            fprintf('   Błąd pobierania GCAT: %s\n', ME.message);
        end
    end
    
    if file_exists
        disp('   Parsowanie danych fizycznych (niskopoziomowe)...');
        
        fid = fopen(local_gcat_file, 'r');
        headerFound = false;
        
        % Domyślne indeksy
        idx_ID = 2; idx_Mass = 16; idx_Len = 22; idx_Dia = 24; idx_Span = 26;
        lineCount = 0;
        dataStartLine = 1;
        
        while ~feof(fid)
            line = fgetl(fid);
            lineCount = lineCount + 1;
            if contains(line, 'Satcat') && contains(line, 'Mass')
                headerFound = true;
                dataStartLine = lineCount + 1;
                headers = split(string(line), sprintf('\t'));
                for k = 1:length(headers)
                    hName = strtrim(headers(k));
                    if hName == "Satcat", idx_ID = k; end
                    if hName == "Mass", idx_Mass = k; end
                    if hName == "Length", idx_Len = k; end
                    if hName == "Diameter", idx_Dia = k; end
                    if hName == "Span", idx_Span = k; end
                end
                break;
            end
            if lineCount > 100, break; end 
        end
        fclose(fid);
        
        opts = detectImportOptions(local_gcat_file, 'FileType', 'text', 'Delimiter', '\t', 'CommentStyle', '');
        opts.DataLine = dataStartLine;
        opts.VariableNamesLine = lineCount;
        opts.PreserveVariableNames = true;
        opts.VariableTypes = repmat({'string'}, 1, length(opts.VariableNames));
        
        rawTable = readtable(local_gcat_file, opts);
        
        gcatTable = table();
        convertCol = @(tbl, idx) str2double(replace(replace(table2array(tbl(:,idx)), "-", "NaN"), " ", ""));
        
        gcatTable.ID = convertCol(rawTable, idx_ID);
        
        numCols = width(rawTable);
        if idx_Mass <= numCols, gcatTable.Mass = convertCol(rawTable, idx_Mass); else, gcatTable.Mass = NaN(height(gcatTable),1); end
        if idx_Len <= numCols, gcatTable.Length = convertCol(rawTable, idx_Len); else, gcatTable.Length = NaN(height(gcatTable),1); end
        if idx_Dia <= numCols, gcatTable.Diameter = convertCol(rawTable, idx_Dia); else, gcatTable.Diameter = NaN(height(gcatTable),1); end
        if idx_Span <= numCols, gcatTable.Span = convertCol(rawTable, idx_Span); else, gcatTable.Span = NaN(height(gcatTable),1); end
        
        gcatTable = gcatTable(~isnan(gcatTable.ID), :);
        if ~isempty(gcatTable)
            [~, uniqueIdx] = unique(gcatTable.ID, 'stable');
            gcatTable = gcatTable(uniqueIdx, :);
        end
        
        fprintf('   Załadowano %d obiektów z katalogu fizycznego.\n', height(gcatTable));
        % NIE USUWA PLIKU!
    end

    %% KROK 2: Pobranie Masowego Katalogu Orbitalnego (TLE)
    disp('2. Sprawdzanie katalogu TLE (Active)...');
    url_tle_all = 'https://celestrak.org/NORAD/elements/gp.php?GROUP=active&FORMAT=TLE';
    local_tle_file = 'active_tle_catalog.txt'; % Stała nazwa
    
    % Tutaj zawsze próbujemy pobrać nowy, chyba że już jest świeży (np. z dzisiaj)
    need_download = true;
    if exist(local_tle_file, 'file')
        finfo = dir(local_tle_file);
        fileAgeDays = now - finfo.datenum;
        if fileAgeDays < 1 % Młodszy niż 1 dzień
            disp('   Plik active_tle_catalog.txt jest świeży. Pomijam pobieranie.');
            need_download = false;
        else
            disp('   Plik TLE jest stary. Próba aktualizacji...');
        end
    end
    
    if need_download
        success = false;
        try
            for attempt = 1:3
                try
                    websave(local_tle_file, url_tle_all, netOpts);
                    % Weryfikacja
                    fileInfo = dir(local_tle_file);
                    if fileInfo.bytes < 1000
                        content = fileread(local_tle_file);
                        if contains(content, '503') || contains(content, 'Service Unavailable')
                            error('Błąd 503.');
                        end
                    end
                    success = true;
                    break;
                catch
                    fprintf('   [Próba %d/3 nieudana] Czekam 5s...\n', attempt);
                    pause(5);
                end
            end
        catch
        end
        
        % MANUAL FALLBACK
        if ~success
            fprintf('\n!!! OSTRZEŻENIE: Nie udało się pobrać TLE automatycznie. !!!\n');
            if exist(local_tle_file, 'file')
                fprintf('   Używam starej wersji pliku, którą znalazłem na dysku.\n');
            else
                fprintf('   Brak pliku na dysku. Pobierz go ręcznie:\n');
                fprintf('   Link: %s\n', url_tle_all);
                fprintf('   Zapisz jako: %s w folderze roboczym.\n', local_tle_file);
                input('   Naciśnij ENTER gdy plik będzie gotowy...', 's');
            end
        end
    end
    
    if ~exist(local_tle_file, 'file')
        error('Brak pliku TLE. Nie można kontynuować.');
    end
    
    tle_text = fileread(local_tle_file);
    lines = splitlines(strtrim(tle_text));
    numSats = floor(length(lines) / 3);
    fprintf('   Dostępne dane orbitalne dla %d obiektów.\n', numSats);
    
    % NIE PARSUJEMY WSZYSTKIEGO TUTAJ DLA EXCELA, JEŚLI CHCESZ TYLKO POBRAĆ.
    % Ale skoro to funkcja "GenerateSpaceDatabase", to tworzymy Excela.
    
    disp('   Przetwarzanie danych do Excela...');
    
    % Prealokacja
    ids = zeros(numSats, 1);
    names = strings(numSats, 1);
    eccs = zeros(numSats, 1);
    incs = zeros(numSats, 1);
    ras = zeros(numSats, 1);
    argps = zeros(numSats, 1);
    mmas = zeros(numSats, 1);
    h_pers = zeros(numSats, 1);
    h_apos = zeros(numSats, 1);
    tle_yys = zeros(numSats, 1);
    tle_doys = zeros(numSats, 1);
    
    mu = 3.986004418e14;
    
    for i = 1:numSats
        idx = (i-1)*3 + 1;
        l0 = lines{idx}; l1 = lines{idx+1}; l2 = lines{idx+2};
        
        names(i) = strtrim(l0);
        ids(i) = str2double(l2(3:7));
        epochStr = l1(19:32);
        tle_yys(i) = str2double(epochStr(1:2));
        tle_doys(i) = str2double(epochStr(3:end));
        
        incs(i) = str2double(l2(9:16));
        ras(i) = str2double(l2(18:25));
        eccs(i) = str2double(['0.' l2(27:33)]);
        argps(i) = str2double(l2(35:42));
        mm = str2double(l2(53:63));
        mmas(i) = mm;
        
        n_rad_s = mm * (2 * pi) / 86400;
        if n_rad_s > 0
            a_km = ((mu / n_rad_s^2)^(1/3)) / 1000;
            h_pers(i) = a_km * (1 - eccs(i)) - 6378.137;
            h_apos(i) = a_km * (1 + eccs(i)) - 6378.137;
        else
            h_pers(i) = NaN; h_apos(i) = NaN;
        end
    end
    
    % Konwersja czasu na PL
    full_years = zeros(numSats, 1);
    mask2000 = tle_yys < 57;
    full_years(mask2000) = 2000 + tle_yys(mask2000);
    full_years(~mask2000) = 1900 + tle_yys(~mask2000);
    dt_utc = datetime(full_years, 1, 1, 'TimeZone', 'UTC') + days(tle_doys - 1);
    dt_pl = dt_utc; dt_pl.TimeZone = 'Europe/Warsaw';
    epochs_pl = string(dt_pl, 'yyyy-MM-dd HH:mm:ss');
    
    % --- POPRAWKA: Dodanie wszystkich parametrów do tabeli ---
    tleTable = table(ids, names, h_pers, h_apos, incs, eccs, ras, argps, mmas, epochs_pl, ...
        'VariableNames', {'ID', 'Nazwa_TLE', 'Perygeum_km', 'Apogeum_km', ...
        'Inklinacja_deg', 'Ekscentrycznosc', 'RAAN_deg', 'Arg_Perygeum_deg', ...
        'Obiegi_na_dzien', 'Czas_PL'});
        
    if ~isempty(tleTable)
        [~, uniqueTleIdx] = unique(tleTable.ID, 'stable');
        tleTable = tleTable(uniqueTleIdx, :);
    end
    
    % Nie usuwamy pliku TLE w trybie offline-first
    % delete(local_tle_file); 

    %% KROK 3: Łączenie i Zapis
    disp('3. Łączenie i zapis Excela (Baza_Satelitow_LEO_2025.xlsx)...');
    
    if exist('gcatTable', 'var')
        finalTable = outerjoin(tleTable, gcatTable, 'Keys', 'ID', 'Type', 'left', 'MergeKeys', true);
    else
        finalTable = tleTable;
    end
    
    if ismember('Name', finalTable.Properties.VariableNames), finalTable.Name = []; end
    
    % Filtrowanie LEO
    rowsLEO = (finalTable.Perygeum_km <= 2000) & (finalTable.Apogeum_km <= 2000);
    finalTable = finalTable(rowsLEO, :);
    
    % Podział na arkusze
    namesUpper = upper(finalTable.Nazwa_TLE);
    maskStar = contains(namesUpper, 'STARLINK');
    maskOneWeb = contains(namesUpper, 'ONEWEB');
    maskConst = contains(namesUpper, {'HULIANWANG', 'QIANFAN', 'KUIPER'});
    maskRest = ~maskStar & ~maskOneWeb & ~maskConst;
    
    filename = 'Baza_Satelitow_LEO_2025.xlsx';
    if exist(filename, 'file'), delete(filename); end
    
    if any(maskStar), writetable(finalTable(maskStar,:), filename, 'Sheet', 'Starlinki'); end
    if any(maskOneWeb), writetable(finalTable(maskOneWeb,:), filename, 'Sheet', 'OneWeb'); end
    if any(maskConst), writetable(finalTable(maskConst,:), filename, 'Sheet', 'Inne_Konstelacje'); end
    if any(maskRest), writetable(finalTable(maskRest,:), filename, 'Sheet', 'Pozostale'); end
    
    disp('==========================================================');
    disp('SUKCES! Baza została podzielona i zapisana.');
end