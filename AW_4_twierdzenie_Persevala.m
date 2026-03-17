%% Sekcja 4:Twierdzenie Parsevala


% Tp = 0.001; N = 2000; n = 0:N-1;
% x = sin(2*pi*5*n*Tp) + 0.5*sin(2*pi*10*n*Tp) + 0.25*sin(2*pi*30*n*Tp);
E_czas = Tp * sum(x.^2);

% wzor (14)
PHI_xx = 1 / (N*Tp) * ( abs(Tp * fft(x)) ).^2

% energia (15b)
E_widmo = sum( PHI_xx);

% energia (15a)
E_N = Tp * sum(x.^2 );

fprintf('Energia czasu: %.6f\n', E_czas); 
fprintf('Energia Parseval (całe widmo): %.6f\n', E_widmo);
fprintf('Energia Parseval (metoda druga): %.6f\n', E_N);
fprintf('Błąd względny: %.2e\n', abs(E_czas - E_widmo)/E_czas); 

% PERIODOGRAM:
figure;
f = (0:N/2-1) * (1/Tp)/N;  % f [Hz]
subplot(1,2,1)
stem(f, PHI_xx(1:N/2));
title('Periodogram PSD \Phi_{xx}(f) - tylko dodatnie częstotliwości'); 
xlabel('Częstotliwość f [Hz]'); ylabel('\Phi_{xx}(f)');
grid on; xlim([0 50]);

% Zakładam, że w workspace są wektory v i e o tej samej długości N i ta sama Tp.
% Sprawdzenie i przypisanie N, Tp jeśli nie istnieją:
if ~exist('Tp','var')
    error('Nie znaleziono zmiennej Tp w workspace.');
end
if ~exist('v','var') || ~exist('e','var')
    error('Nie znaleziono sygnałów v i/lub e w workspace.');
end

% Upewnij się, że wektory kolumnowe:
v = v(:).'; e = e(:).';
Nv = length(v); Ne = length(e);
if Nv ~= N
    warning('Długość v (%d) różna od N (%d). Używam Nv=%d dla obliczeń.', Nv, N, Nv);
end
if Ne ~= N
    warning('Długość e (%d) różna od N (%d). Używam Ne=%d dla obliczeń.', Ne, N, Ne);
end

% Funkcja pomocnicza do obliczeń Parsevala i periodogramu
compute_and_plot = @(sig, signame, sigN) ...
    deal( Tp * sum(sig.^2), ...                                  % energia w czasie
          1/(sigN*Tp) * (abs(Tp * fft(sig))).^2 );                % PHI_xx

% Obliczenia dla v
[E_v_time, PHI_v] = compute_and_plot(v, 'v', Nv);
E_v_spectrum = sum(PHI_v);
E_v_N = Tp * sum(v.^2);

fprintf('\n--- Sygnał v ---\n');
fprintf('Energia czasu (v): %.6f\n', E_v_time);
fprintf('Energia Parseval (v, całe widmo): %.6f\n', E_v_spectrum);
fprintf('Energia Parseval (v, metoda druga): %.6f\n', E_v_N);
fprintf('Błąd względny (v): %.2e\n', abs(E_v_time - E_v_spectrum)/E_v_time);

% Rysunek periodogramu dla v (po prawej stronie subplotu)
subplot(1,2,2);
f_v = (0:floor(Nv/2)-1) * (1/Tp)/Nv;
stem(f_v, PHI_v(1:floor(Nv/2)));
title('Periodogram \Phi_{vv}(f) - tylko dodatnie częstotliwości');
xlabel('Częstotliwość f [Hz]'); ylabel('\Phi_{vv}(f)');
grid on; xlim([0 50]);

% Nowy figure dla e
figure;
[E_e_time, PHI_e] = compute_and_plot(e, 'e', Ne);
E_e_spectrum = sum(PHI_e);
E_e_N = Tp * sum(e.^2);

fprintf('\n--- Sygnał e ---\n');
fprintf('Energia czasu (e): %.6f\n', E_e_time);
fprintf('Energia Parseval (e, całe widmo): %.6f\n', E_e_spectrum);
fprintf('Energia Parseval (e, metoda druga): %.6f\n', E_e_N);
fprintf('Błąd względny (e): %.2e\n', abs(E_e_time - E_e_spectrum)/E_e_time);

% Rysunek periodogramu dla e (dodatnie częstotliwości)
f_e = (0:floor(Ne/2)-1) * (1/Tp)/Ne;
stem(f_e, PHI_e(1:floor(Ne/2)));
title('Periodogram \Phi_{ee}(f) - tylko dodatnie częstotliwości');
xlabel('Częstotliwość f [Hz]'); ylabel('\Phi_{ee}(f)');
grid on; xlim([0 50]);


% % Odpowiedź: periodogram zwykle prezentuje się tylko dla dodatnich częstotliwości,
% % ponieważ sygnał rzeczywisty ma widmo symetryczne wokół zera i dodatnie
% % składowe zawierają pełną informację o gęstości mocy. Poniżej krótka demonstracja:
% %
% % - Jeśli chcemy tylko dodatnie częstotliwości: wyświetlamy indeksy 1:floor(N/2).
% % - Jeśli chcemy całe widmo (oba bieguny): wyświetlamy pełny wektor PHI_xx
% %   i zakres częstotliwości od -Fs/2 do Fs/2.
% %
% % Poniższy kod rysuje oba warianty dla sygnału x z początku skryptu (jeśli istnieje).
% 
% if exist('x','var') && exist('Tp','var')
%     Fs = 1/Tp;
%     Nfft = length(x);
%     PHI_full = 1/(Nfft*Tp) * (abs(Tp * fft(x))).^2;
%     % dodatnie częstotliwości
%     posIdx = 1:floor(Nfft/2);
%     f_pos = (posIdx-1) * Fs / Nfft;
%     figure;
%     subplot(2,1,1);
%     stem(f_pos, PHI_full(posIdx), 'filled');
%     title('Periodogram - tylko dodatnie częstotliwości');
%     xlabel('Częstotliwość [Hz]'); ylabel('\Phi_{xx}(f)'); xlim([0 min(50, Fs/2)]);
%     grid on;
%     % całe widmo z osią od -Fs/2 do Fs/2
%     PHI_shift = fftshift(PHI_full);
%     f_full = (-floor(Nfft/2):ceil(Nfft/2)-1) * Fs / Nfft;
%     subplot(2,1,2);
%     stem(f_full, PHI_shift, 'filled');
%     title('Periodogram - całe widmo (ujemne i dodatnie częstotliwości)');
%     xlabel('Częstotliwość [Hz]'); ylabel('\Phi_{xx}(f)'); xlim([-min(50,Fs/2) min(50,Fs/2)]);
%     grid on;
% else
%     % Jeśli x nie istnieje, informacja tekstowa wystarczy.
%     fprintf('\nDomyślna praktyka: pokazuje się tylko dodatnie częstotliwości dla sygnałów rzeczywistych.\n');
%     fprintf('Aby zobaczyć pełne widmo (ujemne + dodatnie), można użyć fftshift i oś od -Fs/2 do Fs/2.\n');
% end