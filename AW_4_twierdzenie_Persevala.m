%% Sekcja 4:Twierdzenie Parsevala


% Tp = 0.001; N = 2000; n = 0:N-1;
% x = sin(2*pi*5*n*Tp) + 0.5*sin(2*pi*10*n*Tp) + 0.25*sin(2*pi*30*n*Tp);
E_czas = Tp * sum(x.^2);

Xk = fft(x);
Sxx_periodo = (Tp/N) * abs(Xk).^2;


E_widmo = sum(Sxx_periodo);

fprintf('Energia czasu: %.6f\n', E_czas); 
fprintf('Energia Parseval (całe widmo): %.6f\n', E_widmo);
fprintf('Błąd względny: %.2e\n', abs(E_czas - E_widmo)/E_czas); 

figure;
f = (0:N/2-1) * (1/Tp)/N;  % f [Hz]
subplot(1,2,1)
stem(f, Sxx_periodo(1:N/2));
title('Periodogram PSD \Phi_{xx}(f) - tylko dodatnie częstotliwości'); 
xlabel('Częstotliwość f [Hz]'); ylabel('\Phi_{xx}(f)');
grid on; xlim([0 50]);