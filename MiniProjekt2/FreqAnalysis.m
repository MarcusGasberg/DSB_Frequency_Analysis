function [X_out,f_out] = FreqAnalysis(S,T)
%FREQANALYSIS 
%S: signal
%T: Amount of time to analyze
%X_out: Fourier transformation X
%f_out: Frequency Axis
[y, fs] = audioread(S);
fprintf('Frequency Analysis for %s\n', S);
T_stop = T;

%Time axis for 0:T_stop
t=0:1/fs:T_stop-1/fs;

%samples to analyze
s1 = y(1:fs*T_stop);

%Plot signal
figure(1);
subplot(2,2,1);
plot(t,s1);
ylabel('Amplitude');
xlabel('Time (s)');
title(S);

%FFT
X = fft(s1,length(s1));
delta_f = fs/length(s1);
f_axis = 0:delta_f:fs-delta_f;



%Plot frequency spectrum for X
figure(1);
subplot(2,2,2);
semilogx(f_axis(1:0.5*end), 20*log10(abs((2/length(X))*X(1:0.5*end))))
xlabel('Frequency [Hz]')
ylabel('[dB]')
title('DFT magnitude')

hold on
%  Smooths the spectrum in X sampled by fs in 1/M octaves
%  The two-element vector stst contains start and stop frequencies,
[F, X_smooth] = oct_smooth(X, fs, 18,[1 20000]);
semilogx(F, 20*log10(abs((2/length(X))*X_smooth)))
hold off

E_low = 1/length(X)*sum(abs(X(1:round(40/delta_f)).^2));
E_high = 1/length(X)*sum(abs(X(round(40/delta_f):end).^2));
fprintf('Energy of low freqs of s1: %f\n', E_low);
fprintf('Energy of high freqs of s1: %f\n', E_high);
fprintf('Energy relationship: %f\n', E_low/E_high);

%Hanning window
s1_hann = s1.*hann(length(s1))';
X_hann = fft(s1_hann,length(s1_hann));

figure(1);
subplot(2,2,3);
plot(t,s1_hann);
ylabel('Amplitude');
xlabel('Time (s)');
title('Hanning window on signal');

%plot frequency spectrum for X_hann
figure(1);
subplot(2,2,4);
semilogx(f_axis(1:0.5*end), 20*log10(abs((2/length(X_hann))*X_hann(1:0.5*end))))
xlabel('Frequency [Hz]')
ylabel('[dB]')
title('DFT magnitude(Hanning)')

hold on
%  Smooths the spectrum in X sampled by fs in 1/M octaves
%  The two-element vector stst contains start and stop frequencies,
[F, X_smooth_hann] = oct_smooth(X_hann, fs, 18,[1 20000]);
semilogx(F, 20*log10(abs((2/length(X_hann))*X_smooth_hann)))
hold off

X_out = X;
f_out = f_axis;
end

