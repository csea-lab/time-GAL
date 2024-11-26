function[h2]=plot_brain(data,threshold,view1,view2)

% First load the kernel

% This one is ImagingKernel from Minimum Norm Estiamtion (MNE) with eye
% noise cov mat.
load('kernelPlotBrain.mat', 'ImagingKernel')
kernel =  ImagingKernel'; 

% Alternatively, the leadfield could also be used to plot topography
% head = load('headmodel_surf_openmeeg.mat'); 
% kernel = head.Gain


% Constrain the Unconstrained kernel model
c = 1;
for i = 1:3:9004
    tmp = sqrt(kernel(:,i).^2 + kernel(:,i+1).^2 + kernel(:,i+2).^2);
    W(:,c) = tmp;
    c = c+1;
end

% Normalize and multiply Weigths matrix (W) by data
W_norm = W ./ sum(W); 
projection = W_norm' * data;% topography;
vector = projection;


% Define colors
tmpcmp1 = [
    1 1 1
    1 1 0
    1 0 0
    0 0 0
];
cmap_personalize = flipud([linspace(tmpcmp1(1,1), tmpcmp1(2,1), 85)' linspace(tmpcmp1(1,2), tmpcmp1(2,2), 85)' linspace(tmpcmp1(1,3), tmpcmp1(2,3), 85)' ;
    linspace(tmpcmp1(2,1), tmpcmp1(3,1), 85)' linspace(tmpcmp1(2,2), tmpcmp1(3,2), 85)' linspace(tmpcmp1(2,3), tmpcmp1(3,3), 85)';
    linspace(tmpcmp1(3,1), tmpcmp1(4,1), 86)' linspace(tmpcmp1(3,2), tmpcmp1(4,2), 86)' linspace(tmpcmp1(3,3), tmpcmp1(4,3), 86)'; ...
]);

% In case that there is a threshold, we add one final colkor gray to
% colarate as brain surface
if exist('threshold') 
    if length(threshold) > 0
    cmap_personalize = flipud([linspace(tmpcmp1(1,1), tmpcmp1(2,1), 85)' linspace(tmpcmp1(1,2), tmpcmp1(2,2), 85)' linspace(tmpcmp1(1,3), tmpcmp1(2,3), 85)' ;
        linspace(tmpcmp1(2,1), tmpcmp1(3,1), 85)' linspace(tmpcmp1(2,2), tmpcmp1(3,2), 85)' linspace(tmpcmp1(2,3), tmpcmp1(3,3), 85)';
        linspace(tmpcmp1(3,1), tmpcmp1(4,1), 85)' linspace(tmpcmp1(3,2), tmpcmp1(4,2), 85)' linspace(tmpcmp1(3,3), tmpcmp1(4,3), 85)';
        0.5 0.5 0.5]);
    end
end


% Now we normalize the vector of the data to caculate the scales and colors
% for each point in the brain.
vector_norm = (vector - min(vector)) / (max(vector) - min(vector));

% and get the vector of brain points re-scaled again to the incoming data
% so we have  same values as original data.
rescaled = ( vector_norm * ( max(data)-min(data) ) ) + min(data);

% set limits for colorbar, or later modify if threshold differences
cblims = [min(rescaled), max(rescaled)];

% if there is threshold(s), set new color of each brain point based on
% them. Distinctly if they are only one threshold [treshold1] or two
% [treshold1 threshold2]
if exist('threshold')
    if length(threshold) == 2 % [treshold1 threshold2]
        rescaledTH = (rescaled - threshold(1)) / (threshold(2) - threshold(1)); % here we use threshold2 as upper lim
        rescaledTH(rescaledTH < 0) = 0; rescaledTH(rescaledTH > 1) = 1;
        rescaled(rescaled < threshold(1)) = [];
        rescaled(rescaled > threshold(2)) = [];
        cblims = [threshold(1), threshold(2)];
        vector_norm = rescaledTH;
    elseif length(threshold) == 1 % [treshold1] Only down threshold
        rescaledTH = (rescaled - threshold(1)) / (max(rescaled) - threshold(1)); % here we use max cause there is no upper lim
        rescaledTH(rescaledTH < 0) = 0; rescaledTH(rescaledTH > 1) = 1;
        rescaled(rescaled < threshold(1)) = [];
        cblims = [threshold(1), max(rescaled)];
        vector_norm = rescaledTH;
    end
end


% Now get the vector_norm (after thresholding if apply) and select a color 
indices = round(vector_norm * 255) + 1;
color_values = cmap_personalize(indices, :);
%disp(color_values);

% Now we load brain model with faces and vertices
brain = load('brain_cortex_3000V.mat');

%%%% [BRAINSTORM NEEDED HERE]
%%%% IS POSSIBLE TO SMOOTH BRAIN. FUTURE VERSIONS OF THIS SCRIPT MIGHT PRE-CALCULATE
%%%% A FEW SMOOTHED BRAINS AND INCLUDE AS POSSIBLE FILES TO USE.
%%%% smooth_value = 5
%%%% [Vertices_sm, A] = tess_smooth(brain.Vertices, 1, smooth_value, brain.VertConn);
%%%% brain.Vertices = Vertices_sm;
%%%% SMOOOTH_VALUE OF 0 OR NO USING IMPLY NO SMOOTH, 10 IS HIGH SMOOTH.


% Draw brain with the vertices, faces and colors calculated before
h2=patch('vertices',brain.Vertices,'faces',brain.Faces,'FaceVertexCData',color_values,'EdgeColor','none','FaceColor','interp','FaceLighting','phong','AmbientStrength',0.3,'SpecularColorReflectance',0.2);

% Set colormap using the previous created here and set lims based on real data 
colormap(gca, cmap_personalize);
caxis(cblims);
%%%c = colorbar; % Uncomment this to automatically create colorbar on
% function execution

% If there are values for view variables, rotate brain
if exist('view1') && exist('view2')
    view(view1,view2);
end

% Remove grid and axis
grid off;
axis off;

% Set material properties
material([ 0.0 0.9 0.00 1.0 0.0 ]);

% Set lights properties for brian illumination
lighting phong; %gouraud;
hl = [];
hl(1) = camlight(0,40,'local');
hl(2) = camlight(180,40,'local');
hl(3) = camlight(0,-90,'local');
hl(4) = camlight(90,0,'local');
hl(4) = camlight(-90,0,'local');
cb.LineWidth = 0.001;
cb.TickLength = 0;