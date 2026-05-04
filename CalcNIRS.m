function [dHbR, dHbO, fig] = CalcNIRS(dataFile, SDS, tissueType, plotChannelIdx, ...
    extinctionCoefficientsFile, DPFperTissueFile, relDPFfile)

%% ---------- input checks ----------
if ~(ischar(dataFile) || isstring(dataFile))
    error('dataFile must be a filename.');
end
if ~isfile(dataFile)
    error('dataFile does not exist.');
end

if ~isnumeric(SDS) || ~isscalar(SDS) || SDS <= 0
    error('SDS must be a positive scalar.');
end

if ~isnumeric(plotChannelIdx) || any(plotChannelIdx < 1) || any(plotChannelIdx > 20)
    error('plotChannelIdx must contain channel numbers between 1 and 20.');
end

%% ---------- load data ----------
S = load(dataFile);

if ~isfield(S,'SD')
    error('Missing field SD.');
end
if ~isfield(S,'t')
    error('Missing field t.');
end
if ~isfield(S,'d')
    error('Missing field d.');
end

lambda = S.SD.Lambda;
t = S.t(:);
d = S.d;

if numel(lambda) ~= 2
    error('SD.Lambda must contain exactly two wavelengths.');
end
if size(d,2) ~= 40
    error('d must have 40 columns: 20 channels per wavelength.');
end
if size(d,1) ~= length(t)
    error('Length of t must match number of rows in d.');
end

%% ---------- load extinction coefficients ----------
E = readtable(extinctionCoefficientsFile);

% assumes columns: lambda, HbO, HbR
epsHbO = interp1(E.wavelength, E.HbO2, lambda);
epsHbR = interp1(E.wavelength, E.HHb, lambda);

epsMat = [epsHbR(:), epsHbO(:)];   % 2 x 2

%% ---------- DPF ----------
DPFtable = readtable(DPFperTissueFile);

idx = strcmpi(DPFtable.Tissue, tissueType);
if ~any(idx)
    error('tissueType not found in DPF file.');
end

DPF = DPFtable.DPF(idx);
DPF = DPF(1);

if isscalar(DPF)
    DPF = [DPF; DPF];
end

%% ---------- modified Beer-Lambert ----------
I1 = d(:,1:20);      % first wavelength
I2 = d(:,21:40);     % second wavelength

I01 = mean(I1(1:10,:),1);
I02 = mean(I2(1:10,:),1);

OD1 = -log10(I1 ./ I01);
OD2 = -log10(I2 ./ I02);

dHbR = zeros(size(I1));
dHbO = zeros(size(I1));

L = SDS; % cm

for ch = 1:20
    OD = [OD1(:,ch), OD2(:,ch)]';   % 2 x time
    
    A = diag(DPF * L) * epsMat;     % 2 x 2
    conc = A \ OD;                  % 2 x time
    
    dHbR(:,ch) = conc(1,:)';
    dHbO(:,ch) = conc(2,:)';
end

%% ---------- plot selected channels ----------
fig = figure;

for k = 1:length(plotChannelIdx)
    ch = plotChannelIdx(k);
    
    subplot(length(plotChannelIdx),1,k)
    plot(t, dHbR(:,ch), 'b', 'LineWidth', 1.2); hold on;
    plot(t, dHbO(:,ch), 'r', 'LineWidth', 1.2);
    xlabel('Time [sec]');
    ylabel('\Delta Hb');
    title(['Channel ', num2str(ch)]);
    legend('\DeltaHbR','\DeltaHbO');
    grid on;
end

end
