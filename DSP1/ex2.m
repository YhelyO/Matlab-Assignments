function [] = ex2()
%EX2 Summary of this function goes here
%   Detailed explanation goes here
load('filter_0.25_101.mat');

Fs = 60; % according to 2.1
t = 0:1/Fs:34.12; % 2048 samples
r_n = cos(2*pi*5*t)+cos(2*pi*10*t);
 
R_f = fftshift(fft(r_n));
f = (0:length(R_f)-1)*Fs/length(R_f);
h1 = [h zeros(1,1946)];
H_f = fftshift(fft(h1));
 
 
S_f = R_f .* H_f;
 
figure3 = figure ; 
 
subplot(3,1,1)
plot(abs(R_f));
title('R_f');
xlabel('f(Hz)');
ylabel('R_f');
axis tight;
 
subplot(3,1,2)
plot(abs(H_f));
title('H_f');
xlabel('f(Hz)');
ylabel('H_f');
axis tight;
 
subplot(3,1,3)
plot(f,abs(S_f));
title('S_f');
xlabel('f(Hz)');
ylabel('S_f');
axis tight;
 end

