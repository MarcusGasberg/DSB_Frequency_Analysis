%Miniprojekt 2 - Wind turbine
%%
clear;
clc;
[y_turb, fs_turb] = audioread('Wind turbine sound - high quality audio.wav');

T_stop = 10;

%Time axis for first 10 sec
t_turb=0:1/fs_turb:T_stop-1/fs_turb;

%samples for first 10 sec
s1 = y_turb(1:fs_turb*T_stop);

%Plot turbine
plot(t_turb,s1);
ylabel('Amplitude');
xlabel('Time (s)');

%% Frequency analysis for wind turbine

%fft
X_turb = fft(s1,length(s1));
delta_f = fs_turb/length(s1);
f_axis = 0:delta_f:fs_turb-delta_f;


%Plot for signals and frequency spectrum
figure(1);
%Plot s2
subplot(2,2,1);
plot(t_turb,s1);
xlabel('Time [s]')
ylabel('Amplitude')
title('S2')

%Plot frequency spectrum for X
subplot(2,2,2);
semilogx(f_axis(1:0.5*end), 20*log10(abs((2/length(X_turb))*X_turb(1:0.5*end))))
xlabel('Frequency [Hz]')
ylabel('[dB]')
title('DFT magnitude')

%Calculate and print low and high frequency energy
fprintf('\nEnergy of signals for X\n');
E_low = 1/length(X_turb)*sum(abs(X_turb(1:round(40/delta_f)).^2));
E_high = 1/length(X_turb)*sum(abs(X_turb(round(40/delta_f):end).^2));
fprintf('Energy of low freqs of s1: %f\n', E_low);
fprintf('Energy of high freqs of s1: %f\n', E_high);
fprintf('Energy relationship: %f\n', E_low/E_high);

%Hanning window
s1_hann = s1.*hann(length(s1))';
X_turb_hann = fft(s1_hann,length(s1_hann));

%plot s2 hanning
subplot(2,2,3);
plot(t_turb,s1_hann);
xlabel('Time [s]')
ylabel('Amplitude')
title('S2 with hanning')

%plot frequency spectrum for X_hann
subplot(2,2,4);
semilogx(f_axis(1:0.5*end), 20*log10(abs((2/length(X_turb_hann))*X_turb_hann(1:0.5*end))))
xlabel('Frequency [Hz]')
ylabel('[dB]')
title('DFT magnitude(Hanning)')


fprintf('\nEnergy of signal X_hann\n');
E_low = 1/length(X_turb_hann)*sum(abs(X_turb_hann(1:round(40/delta_f)).^2));
E_high = 1/length(X_turb_hann)*sum(abs(X_turb_hann(round(40/delta_f):end).^2));

fprintf('Energy of low freqs of s1: %f\n', E_low);
fprintf('Energy of high freqs of s1: %f\n', E_high);
fprintf('Energy relationship: %f\n', E_low/E_high);

%% Dial up

[y_dial, fs_dial] = audioread('Phone Dialing Sound Effect.wav');
figure(2);
img = imread('DTMF_Freq.PNG');
image(img); %Dial up freqs

T_stop = 5;

%Time axis for first 10 sec
t_dial=0:1/fs_dial:T_stop-1/fs_dial;

%samples for first 10 sec
s2 = y_dial(1:fs_dial*T_stop);

%fft
X_dial = fft(s2,length(s2));
delta_f = fs_dial/length(s2);
f_axis = 0:delta_f:fs_dial-delta_f;

%Plot for signals and frequency spectrum
figure(3);
%Plot s2
subplot(2,1,1);
plot(t_dial,s2);
xlabel('Time [s]')
ylabel('Amplitude')
title('S2')

%Plot frequency spectrum for X
subplot(2,1,2);
semilogx(f_axis(1:0.5*end), 20*log10(abs((2/length(s2))*X_dial(1:0.5*end))))
xlabel('Frequency [Hz]')
ylabel('[dB]')
title('DFT magnitude')

%% Music

%Show case of function
FreqAnalysis('Old RuneScape Soundtrack Sea Shanty2.wav',30);
pause
FreqAnalysis('Wind turbine sound - high quality audio.wav',10);
pause
FreqAnalysis('Phone Dialing Sound Effect.wav',5);


