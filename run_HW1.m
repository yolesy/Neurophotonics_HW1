clear; clc; close all;

SDS = 3;
tissueType = 'adult_head';
plotChannelIdx = [1 2];

file1 = 'FN_031_V2_Postdose2_Nback.mat';
file2 = 'FN_032_V1_Postdose1_Nback.mat';

extFile = 'ExtinctionCoefficientsData.csv';
dpfFile = 'DPFperTissue.txt';
relFile = 'RelativeDPFCoefficients.csv';

[dHbR1, dHbO1, fig1] = CalcNIRS(file1, SDS, tissueType, plotChannelIdx, ...
    extFile, dpfFile, relFile);

[dHbR2, dHbO2, fig2] = CalcNIRS(file2, SDS, tissueType, plotChannelIdx, ...
    extFile, dpfFile, relFile);


% -------- Fourier + SNR (first file, channel 1) --------

load('FN_031_V2_Postdose2_Nback.mat')

Fs = 1 / mean(diff(t));

x = dHbO1(:,1);   % you can also try dHbR1(:,1)
x = x - mean(x);

N = length(x);
Y = abs(fft(x));
f = (0:N-1) * Fs/N;

% take half spectrum
half = 1:floor(N/2);
f = f(half);
Y = Y(half);

figure;
plot(f, Y);
xlabel('Frequency [Hz]');
ylabel('Magnitude');
title('Fourier - Channel 1');
grid on;

% -------- SNR --------

% heart rate ~0.8–2 Hz
heartIdx = f >= 0.8 & f <= 2;
[signalAmp, idx] = max(Y(heartIdx));
heartFreqs = f(heartIdx);
heartFreq = heartFreqs(idx);

noiseAmp = mean(Y(f > 2.5));

SNR = signalAmp / noiseAmp;

fprintf('Heart frequency = %.3f Hz\n', heartFreq);
fprintf('SNR = %.3f\n', SNR);