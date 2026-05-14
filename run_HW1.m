clear; clc; close all;

SDS = 3;
tissueType = 'adult_head';
plotChannelIdx = [1 2];

file1 = 'FN_031_V2_Postdose2_Nback.mat';
file2 = 'FN_032_V1_Postdose1_Nback.mat';

extFile = 'ExtinctionCoefficientsData.csv';
dpfFile = 'DPFperTissue.txt';
relFile = 'RelativeDPFCoefficients.csv';

%% -------- Run CalcNIRS on first file --------
[dHbR1, dHbO1, fig1] = CalcNIRS(file1, SDS, tissueType, plotChannelIdx, ...
    extFile, dpfFile, relFile);

exportgraphics(fig1, 'Hb_channels_file1.tif', 'Resolution', 600);

%% -------- Run CalcNIRS on second file --------
[dHbR2, dHbO2, fig2] = CalcNIRS(file2, SDS, tissueType, plotChannelIdx, ...
    extFile, dpfFile, relFile);

exportgraphics(fig2, 'Hb_channels_file2.tif', 'Resolution', 600);

%% -------- Fourier + SNR, first file, channel 1 --------
load(file1)

Fs = 1 / mean(diff(t));

x = dHbO1(:,1);   % Channel 1, oxygenated hemoglobin
x = x - mean(x);  % Remove DC component

N = length(x);
Y = abs(fft(x));
f = (0:N-1) * Fs/N;

% Take half spectrum
half = 1:floor(N/2);
f = f(half);
Y = Y(half);

%% -------- SNR calculation --------
% Heart-rate frequency range
heartIdx = f >= 0.8 & f <= 2;

[signalAmp, idx] = max(Y(heartIdx));
heartFreqs = f(heartIdx);
heartFreq = heartFreqs(idx);

% Noise is defined as the average magnitude above 2.5 Hz
noiseAmp = mean(Y(f > 2.5));

SNR = signalAmp / noiseAmp;

%% -------- Plot Fourier with peak annotation --------
figFFT = figure;

plot(f, Y, 'LineWidth', 1.2);
hold on;

plot(heartFreq, signalAmp, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);

text(heartFreq, signalAmp, sprintf('  Peak = %.2f Hz', heartFreq), ...
    'VerticalAlignment', 'bottom', ...
    'HorizontalAlignment', 'left');

xlabel('Frequency [Hz]');
ylabel('Magnitude [a.u.]');
title('Fourier - Channel 1');
grid on;
xlim([0 4]);

exportgraphics(figFFT, 'FFT_Channel1.tif', 'Resolution', 600);

%% -------- Print results --------
fprintf('Heart frequency = %.3f Hz\n', heartFreq);
fprintf('Heart rate = %.2f BPM\n', heartFreq * 60);
fprintf('SNR = %.3f\n', SNR);
