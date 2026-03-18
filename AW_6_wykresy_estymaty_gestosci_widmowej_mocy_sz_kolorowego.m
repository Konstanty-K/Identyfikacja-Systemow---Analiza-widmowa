%% PSD szumu kolorowego v(nTp)
Tp=0.001; N=2000; n=0:N-1;
sigma=0.8; e=sigma*randn(1,N);

H = tf([0.1],[1 -0.9],Tp);
v = lsim(H, e, n*Tp)';

f = (0:N/2-1)*(1/Tp)/N;

% 1. Periodogram (14) - POPRAWIONY
Vk = Tp * fft(v);                       % X_N = Tp * fft(v)
Phi_p = (1/(N*Tp)) * abs(Vk).^2;       % = (Tp/N) * abs(fft(v))^2

% 2. Korelogram (16) - BEZ PĘTLI
Mw = round(N/5);
tau = -Mw:Mw;                           % wektor przesunięć τ
r_vv = xcorr(v, Mw, 'biased');          % r_vv(τ), dł. 2*Mw+1
w = 0.5*(1 + cos(pi*tau/Mw));           % okno Hanninga (18): 1+cos !
Phi_k = Tp * abs(fft(r_vv.*w, N));      % DFT sumy z wzoru (16)

figure;
subplot(2,1,1);
plot(f, Phi_p(1:N/2),'b-','LineWidth',1.5);
title('Periodogram (14)');
xlabel('f [Hz]'); ylabel('\Phi_{vv}'); grid on; xlim([0 250]);

subplot(2,1,2);
plot(f, Phi_k(1:N/2),'r-','LineWidth',1.5);
title('Korelogram (16)');
xlabel('f [Hz]'); ylabel('\Phi_{vv}'); grid on; xlim([0 250]);

