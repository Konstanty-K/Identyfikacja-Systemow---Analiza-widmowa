%% Sekcja 5: Metoda korelogramowa (okno Hanninga, Mw=N/5)
rng(77,'twister'); % stałe ziarno generatora liczb losowych
Tp = 0.001;
N=2000;
Mw = round(N/6);  % Mw=400 próbek
N = 2*Mw+1;
tau = (-Mw) : Mw;
n = 0:N-1;
sigma = 0.8;
x = sin(2*pi*5*n*Tp) + 0.5*sin(2*pi*10*n*Tp) + 0.25*sin(2*pi*30*n*Tp);
e = sigma * randn(1,N);
H = tf([0.1], [1, -0.9], Tp);
v = lsim(H, e, n*Tp)';
% DLA porównania z mniejszym okneml modyfikujemy Mw

x_DFT = Tp * fft(x);
e_DFT = Tp * fft(e);
v_DFT = Tp * fft(v);
% rectangle window (rectangular/boxcar) of length Mw (centered on 0..2*Mw)
% w_p should be length 2*Mw+1 to match rxx indexing (which is 2*Mw+1)
%w_p = ones(1, 2*Mw+1);  % prostokątne okno o wartości 1
w_p = zeros(2*Mw+1, 1);
for i = 1 : 2*Mw+1
    %tau_i = i - Mw;
    %if abs(tau_i) <= Mw
    w_p(i) = 1;  % długość 2*Mw+1
    %end
end
%w_p = ones(size(tau)); 

% autocorrelation 
rxx = xcorr(x, x, Mw, 'biased');  % długość 2*Mw+1, odpowiada indeksom -Mw: Mw
ree = xcorr(e, e, Mw, 'biased');
rvv = xcorr(v, v, Mw, 'biased');

corr_x = rxx;  
corr_e = ree; 
corr_v = rvv; 

PHI_s_xx_p = Tp * sum(w_p .* corr_x .* x_DFT);
PHI_s_ee_p = Tp * sum(w_p .* corr_e .* e_DFT);
PHI_s_vv_p = Tp * sum(w_p .* corr_v .* v_DFT);

% Hanning window (okno Hanninga) of length 2*Mw+1, 
w_h = zeros(N, 1);
for i = 1 : 2*Mw+1
    %tau_i = i - Mw;
    %if abs(tau_i) <= Mw
    w_h(i) = 0.5 * (1 + cos(pi*i/Mw)); %tau_i/Mw));  % długość 2*Mw+1
    %end
end

PHI_s_xx_h = Tp * sum(w_h .* corr_x .* x_DFT); % 3 sinusy
PHI_s_ee_h = Tp * sum(w_h .* corr_e .* e_DFT); % szum biały
PHI_s_vv_h = Tp * sum(w_h .* corr_v .* v_DFT); % szum kolorowy
%%
k = -Mw:Mw;  % przesunizenia τ
lags = k;

figure;

% Convert DFT frequency bins to Hz
fs = 1/Tp;
f = (0:N-1)*(fs/N);

% Plot PHI_s_xx_h (rectangular vs Hanning)
subplot(3,2,1);
plot(f, abs(PHI_s_xx_p), 'b-');
xlim([0 fs/2]); title('PHI\_s\_xx (rectangular)');
xlabel('Frequency (Hz)'); ylabel('\Phi_{s,xx}(f)'); grid on;

subplot(3,2,2);
plot(f, abs(PHI_s_xx_h), 'r-');
xlim([0 fs/2]); title('PHI\_s\_xx (Hanning)');
xlabel('Frequency (Hz)'); ylabel('\Phi_{s,xx}(f)'); grid on;

% Plot PHI_s_ee_h (rectangular vs Hanning)
subplot(3,2,3);
plot(f, abs(PHI_s_ee_p), 'b-');
xlim([0 fs/2]); title('PHI\_s\_ee (rectangular)');
xlabel('Frequency (Hz)'); ylabel('\Phi_{s,ee}(f)'); grid on;

subplot(3,2,4);
plot(f, abs(PHI_s_ee_h), 'r-');
xlim([0 fs/2]); title('PHI\_s\_ee (Hanning)');
xlabel('Frequency (Hz)'); ylabel('\Phi_{s,ee}(f)'); grid on;

% Plot PHI_s_vv_h (rectangular vs Hanning)
subplot(3,2,5);
plot(f, abs(PHI_s_vv_p), 'b-');
xlim([0 fs/2]); title('PHI\_s\_vv (rectangular)');
xlabel('Frequency (Hz)'); ylabel('\Phi_{s,vv}(f)'); grid on;

subplot(3,2,6);
plot(f, abs(PHI_s_vv_h), 'r-');
xlim([0 fs/2]); title('PHI\_s\_vv (Hanning)');
xlabel('Frequency (Hz)'); ylabel('\Phi_{s,vv}(f)'); grid on;

sgtitle('Porównanie PSD: prostokątne (niebieski) vs Hanning (czerwony)');
