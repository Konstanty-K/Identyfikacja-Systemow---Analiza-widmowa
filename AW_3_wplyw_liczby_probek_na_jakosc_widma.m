%% Sekcja 3: wpływ liczby N na widmo sygnału

Tp = 0.001; fs = 1/Tp;
Ns = [1000, 200, 100];
okresy_5Hz = Ns * 5 * Tp;  % liczba okresów składowej 5 Hz w oknie

figure('Position',[100 100 1200 400]);
for i = 1:length(Ns)
    N = Ns(i);
    n = 0:N-1;
    x = sin(2*pi*5*n*Tp) + 0.5*sin(2*pi*10*n*Tp) + 0.25*sin(2*pi*30*n*Tp);
    X = Tp * fft(x);
    f = (0:N-1)*fs/N;
    
    subplot(1,3,i); 
    stem(f(1:N/2), 2*abs(X(1:N/2))/N,  'LineWidth',1.5);
    title(sprintf('N=%d (%.1f okresów 5Hz)', N, okresy_5Hz(i)));
    xlabel('f [Hz]'); ylabel('|X(f)|');
    xlim([0 50]); grid on;
    legend('Piki: 5,10,30 Hz');
end
sgtitle('Wpływ długości okna N na widmo amplitudowe DFT');

%% pojawiły się dodatkowe prążki widma - świadczy to o "wycieku"