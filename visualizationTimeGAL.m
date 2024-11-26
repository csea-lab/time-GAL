function visualizationTimeGAL(timeGALoutput, varargin)
% Visualization of results from TimeGALdecoding. This function creates 3
%   Figures:
%
%       1. Matrices of GAL and Correlation, i.e. the backward and forward model,
%       respectively.
%
%       2. Circular Connectivity plots for direct and inverse GAL and brain
%       topographical reconstruction of decoding rate.
%       
%       3. Time-GAL and topography across time, describing the temporal
%       dynamics of the 
%
%   Parameters:
%       - 'Title' = '', title for each figure.
%       - 'TimeSteps' = 10, indicate how many time-GAL windows visualize.
%       - 'AlphaGAL' = [], indicate a new alphaGAL different from the one contained at the timeGALoutput. 
%       - 'AlphaCorrelation' = [], indicate a new alphaGAL different from the one contained at the timeGALoutput. 
%       - 'FileName' = '', name for the file saving the visualization figures.
%       - 'Resolution' = [], number of dots per inch (DPI).
%       - 'Figures' = [1 2 3], indicates which figures to show.
%       - 'TopoPercentile'= 20, indicate how much % of the topography is shown in the scale.
%
%  Use:
%       visualization(timeGALoutput)
%

%% Arguments
arg = inputParser;
addParameter(arg, 'Title', '', @ischar); %
addParameter(arg, 'TimeSteps', 10, @numeric); %
addParameter(arg, 'AlphaGAL', [], @numeric); %
addParameter(arg, 'AlphaCorrelation', [], @numeric); %
addParameter(arg, 'Resolution', 150, @numeric); %
addParameter(arg, 'FileName', '', @ischar); %
addParameter(arg, 'TopoPercentile', 20, @numeric); %
addParameter(arg, 'Figures', [1 2 3], @numeric); %



% Access to arguments
parse(arg, varargin{:});
titleText = arg.Results.Title;
steps = arg.Results.TimeSteps;
filename = arg.Results.FileName;
alphaGAL = arg.Results.AlphaGAL;
alphaCorr = arg.Results.AlphaCorrelation;
topoPrctile = arg.Results.TopoPercentile;
figures2plot = arg.Results.Figures;
picResolution = arg.Results.Resolution;

% Parameters
chanceLevel = 0.5;
channels = timeGALoutput.Parameters.Channels;
screenSize = get(0, 'ScreenSize'); % [x y width height]


% Take data from timeGALoutput
GAL = timeGALoutput.GeneralizationMatrix.GAL;

if isempty(alphaGAL)
    % Use masks in timeGALoutput
    maskpos = timeGALoutput.GeneralizationMatrix.GALmaskPos;
    maskneg = timeGALoutput.GeneralizationMatrix.GALmaskNeg;
    maskposneg = timeGALoutput.GeneralizationMatrix.GALmask;
else
    % T-Test contrast with new alpha values
    alpha_corrected = alphaGAL/125;
    [~,p] = ttest(GAL(:,channels,channels), chanceLevel, 'tail', 'right');
    maskpos = squeeze(p < alpha_corrected);
    [~,p] = ttest(GAL(:,channels,channels), chanceLevel, 'tail', 'left');
    maskneg = (squeeze(p < alpha_corrected) );
    maskposneg = maskpos + (maskneg .* -1);
end
if isempty(alphaCorr)
    alphaCorr = timeGALoutput.Parameters.alphaCorrelation;
end
CorrelationR = timeGALoutput.CorrelationMatrix.CorrelationR;
CorrelationP = timeGALoutput.CorrelationMatrix.CorrelationP;




%% FIG 1. Matrices

if ismember(1,figures2plot)

proportion = [0.1, 0.2, 0.8, 0.5]; % [x, y, width, height]
figurePosition = [screenSize(1) + proportion(1) * screenSize(3), ...
                  screenSize(2) + proportion(2) * screenSize(4), ...
                  proportion(3) * screenSize(3), ...
                  proportion(4) * screenSize(4)];
Fig1 = figure('Position', figurePosition);

subplot(1,3,1)
imagesc(squeeze(mean(GAL(:,channels,channels)))), %colorbar
title('GAL matrix (Backward model)')
h_ax = gca;
c = colorbar;
c.Position = [h_ax.Position(1) + h_ax.Position(3) + 0.01, ...
              h_ax.Position(2), 0.017, h_ax.Position(4)];
h_ax.XLabel.String = 'Channels';
h_ax.YLabel.String = 'Channels';
cb = [  [linspace(192, 255, 100)' linspace(100, 255, 100)' linspace(100, 255, 100)'] / 255;  ones(1,3)  ;  [linspace(255, 1, 100)' linspace(255, 1, 100)' linspace(255, 1, 100)'] / 255;  ];
h_ax.Colormap = cb;
minmaxs = abs([min(min(mean(GAL(:,channels,channels)))), max(max(mean(GAL(:,channels,channels))))]);
clim([0.5 - (max(minmaxs) - 0.5) , 0.5 + (max(minmaxs) - 0.5) ])
h_ax_c = axes('position', get(h_ax, 'position'), 'Color', 'none');
[M, c] = contour(flipud(maskpos), [0.5 0.5]);%, 'ShowText', 'on');
c.LineWidth = 2;
h_ax_c.Color = 'none';
h_ax_c.Colormap = [0 0 0];
h_ax_c.XTick = [];
h_ax_c.YTick = [];
h_ax_c2 = axes('position', get(h_ax, 'position'), 'Color', 'none');
[M2, c2] = contour(flipud(maskneg), [0.5 0.5]);%, 'ShowText', 'on');
c2.LineWidth = 2;
h_ax_c2.Color = 'none';
h_ax_c2.Colormap = [1 0 0];
h_ax_c2.XTick = [];
h_ax_c2.YTick = [];
h_ax.Box = 'off';
hold on; plot(fliplr(1:0.1:125), 1:0.1:125, 'Color',[0.5 0.5 0.5], 'LineWidth', 1)
h_ax = gca; h_ax.Box = 'off';

% Correlation matrix
subplot(1,3,2:3)
imagesc(CorrelationR(channels,:)), colorbar
title('Correlation matrix (Forward model)')
h_ax = gca;
h_ax.XLabel.String = 'Time (s)';
h_ax.YTick = [];
h_ax.Colormap = cb;
minmaxs= abs([ min(min(CorrelationR(channels,:))), max(max(CorrelationR(channels,:)))  ] );
clim([max(minmaxs) * -1 , max(minmaxs)])
h_ax.Box = 'off';
hold on, 
[M3, c3] = contour(CorrelationP(channels,:) < alphaCorr, 'k');
c3.LineWidth = 0.3;
h_ax_c3.Color = 'none';
h_ax_c3.Colormap = [0.5 0 0];

if ~isempty(titleText)
    sgtitle(titleText);
end
if ~isempty(filename)
    exportgraphics(gcf, [filename '.Matrices.png'], 'Resolution', picResolution);
end

end




%% FIG 2. Connectivity Circles and Brain

if ismember(2,figures2plot)

proportion = [0.0136, 0.3067, 1.5789, 0.4200];
figurePosition = [proportion(1) * screenSize(3), ...
                  proportion(2) * screenSize(4), ...
                  proportion(3) * screenSize(3), ...
                  proportion(4) * screenSize(4)];

Fig2 = figure('Position', figurePosition);

subplot(1,3,1)
circularConnectivity125static(squeeze(mean(GAL(:,channels,channels))), maskpos(:,:), 1, 'Gradient')
subplot(1,3,2)
circularConnectivity125static(squeeze(mean(GAL(:,channels,channels))), maskneg(:,:) .* -1, 1, 'Gradient')
subplot(1,3,3)
topo = diag(squeeze(mean(GAL))); topo(125:128) = 0.5; %chance level
plot_brain(topo, prctile(topo, 100-topoPrctile),-90,90); % prctile ensure that only 20% higher topography is displayed
h_ax = gca();
c = colorbar;
c.Position = [h_ax.Position(1) + h_ax.Position(3) + 0.01, ...
              h_ax.Position(2), 0.017, h_ax.Position(4)];

if ~isempty(titleText)
    sgtitle(titleText);
end
if ~isempty(filename)
    exportgraphics(gcf, [filename '.Connectivity.png'], 'Resolution', picResolution);
end

end



%% FIG 3. Time GAL

if ismember(3,figures2plot)

proportion = [0.0, 0.2393, 1.8040, 0.4874];
figurePosition = [proportion(1) * screenSize(3), ...
                  proportion(2) * screenSize(4), ...
                  proportion(3) * screenSize(3), ...
                  proportion(4) * screenSize(4)];
Fig3 = figure('Position', figurePosition);

maskGAL = abs(maskposneg);
maskTime = CorrelationP(channels,:) < alphaCorr;
% look for those channels with any statistically significant GAL connection
indcs = find(sum(maskGAL)>=1);
[row, col] = find(maskGAL);
TimeGAL = zeros(125,125,size(maskTime,2));
% looping over channels, we mark those indicated by time correlation matrix
for i = 1:length(row)
    TimeGAL(row(i),col(i), maskTime(row(i),:) == 1) = 1;
end
% and apply the value (positive/negative)
TimeGAL = TimeGAL.*maskposneg; %%%(channels,channels);

TopographyAcrossTime = zeros([129, size(CorrelationR, 2)]);
for t = 1:size(CorrelationR,2)
    TopographyAcrossTime(channels,t) = abs(CorrelationR(channels,t)) .* diag(squeeze(mean(GAL(:,channels,channels))));
end

vec = 1:size(TimeGAL,3); % Tu vector original
numSegments = steps;
segmentSize = floor(length(vec) / numSegments); % Tamaño base de cada segmento
remainder = mod(length(vec), numSegments); % Elementos sobrantes

segments = cell(1, numSegments);
startIndex = 1;

for i = 1:numSegments
    if i <= remainder
        % Añadir un elemento extra a los primeros 'remainder' segmentos
        endIndex = startIndex + segmentSize;
    else
        endIndex = startIndex + segmentSize - 1;
    end
    segments{i} = vec(startIndex:endIndex);
    startIndex = endIndex + 1;
end

for st = 1:steps
    timeGAL2plot = sum(TimeGAL(:,:,segments{st}),3);
    timeGAL2plotPos = timeGAL2plot; timeGAL2plotPos(timeGAL2plot < 0) = 0;
    timeGAL2plotNeg = timeGAL2plot; timeGAL2plotNeg(timeGAL2plot > 0) = 0;
    subplot(4,steps,st)
    circularConnectivity125times(timeGAL2plotPos, length(segments{st}), 0, 0, 'Gradient')
    
    subplot(4,steps,st + steps)
    circularConnectivity125times(timeGAL2plotNeg, length(segments{st}), 0, 0, 'Gradient')

    Topographies2plot(:,st) = mean(TopographyAcrossTime(:,segments{st}), 2);
    subplot(4,steps,st + steps * 2)
    plot_brain(Topographies2plot(:,st), [prctile(TopographyAcrossTime(:), 100-topoPrctile) prctile(TopographyAcrossTime(:), 100)], 180, 0);
    
    subplot(4,steps,st + steps * 3)
    plot_brain(Topographies2plot(:,st), [prctile(TopographyAcrossTime(:), 100-topoPrctile) prctile(TopographyAcrossTime(:), 100)], 0, 0);
end

if ~isempty(titleText)
    sgtitle(titleText);
end
if ~isempty(filename)
    exportgraphics(Fig3, [filename '.timeGAL.png'], 'Resolution', picResolution);
end

end






