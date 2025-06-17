%% Operations_Exp2.m
% This program is used to generate plots for papers/presentations.
clc; clear; close all;

one_px = 12.446e-6; % meters

dataDir = 'E:/SFF/Exp2_MeasuringTracks_DiffStartsStops/1_ConstantLaserPower_40W_0o5mms_fp_7mm/';
load(fullfile(dataDir, 'XZ_profile_visual_allFrames.mat'));
load(fullfile(dataDir, 'XZ_profile_confocal_allFrames.mat'));

X_conf = X_2d;
Z_conf = Z_2d;

clear Z_2d X_2d;

% Create figure for raw points and processed curves
fig1 = figure('Color', 'w', 'Units', 'inches', 'Position', [1, 1, 12, 5]);
hold on;

% Plot a representative gray marker for legend (before loop)
hGray = plot(nan, nan, 'o', 'MarkerEdgeColor', [0.6 0.6 0.6], ...
    'MarkerSize', 5, 'DisplayName', 'Raw data');

numFrames = height(results);
initial_shift = 49; % mm
increment = 0.5;          % mm per frame

allX = [];
allZ = [];

for i = 1:numFrames
    X_vis = results.X_vis{i};
    Z_vis = results.Z_vis{i};

    targetMaxX = initial_shift - (i - 1) * increment;
    currentMaxX = max(X_vis);
    X_shifted = X_vis - (currentMaxX - targetMaxX);

    allX = [allX; X_shifted(:)];
    allZ = [allZ; Z_vis(:)];

    % Raw scatter plot (light gray, small markers)
    plot(X_shifted, Z_vis, 'o', ...
        'MarkerEdgeColor', [0.6 0.6 0.6], ...
        'MarkerSize', 5, ...
        'HandleVisibility', 'off'); 
end

% Median filtering
xRounded = round(allX, 1); % round to nearest 0.1mm for binning
uniqueX = unique(xRounded);
medianZ = zeros(size(uniqueX));

for j = 1:length(uniqueX)
    mask = xRounded == uniqueX(j);
    medianZ(j) = median(allZ(mask));
end

% Plot confocal data (solid blue)
plot(X_conf, Z_conf, '-', 'Color', [0 0.45 0.74], 'LineWidth', 1.5, 'DisplayName', 'Confocal data');

% Plot visual data (dashed dark red)
plot(uniqueX, medianZ, '--', 'Color', [0.65 0 0], 'LineWidth', 1.5, 'DisplayName', 'Visual data');

% Formatting for academic presentation
xlabel('X (mm)');
ylabel('Z (mm)');
ylim([0 2]);
yticks(0:0.25:2);
grid on;
box on;
set(gca, 'FontSize', 16);
legend('Location', 'best', 'FontSize', 14);

% Save first figure
set(fig1, 'PaperUnits', 'inches', 'PaperPosition', [0 0 12 5]);
exportgraphics(fig1, fullfile(dataDir, 'XZ_ConfocalVisualHeight.png'), 'Resolution', 300);
% Uncomment for PDF:
% exportgraphics(fig1, fullfile(dataDir, 'Figure1_ProfileComparison.pdf'), 'ContentType', 'vector');

%% Height Variation Between Visual and Confocal Profiles
% Interpolate confocal Z values at visual X locations
Z_conf_interp = interp1(X_conf, Z_conf, uniqueX, 'linear', NaN);

% Calculate the difference (error) between visual and confocal
Z_diff = medianZ - Z_conf_interp;

% Remove NaNs (from out-of-range interpolation)
valid = ~isnan(Z_conf_interp);
mae = mean(abs(Z_diff(valid)));
rmse = sqrt(mean(Z_diff(valid).^2));
max_error = max(abs(Z_diff(valid)));

fig2 = figure('Color', 'w', 'Units', 'inches', 'Position', [1, 1, 12, 5]);
plot(uniqueX(valid), Z_diff(valid), '-', 'LineWidth', 1.5);
yline(0, '--k', 'LineWidth', 1.5);
xlabel('X (mm)');
ylabel('\DeltaZ (mm)');
xlim([0 50]);
grid on;
set(gca, 'FontSize', 16);

% Add annotation in bottom-right corner
statsText = sprintf(['Mean Abs Error: %.3f mm\n' ...
                     'RMSE: %.3f mm\n' ...
                     'Max Error: %.3f mm'], mae, rmse, max_error);
text(0.98, 0.02, statsText, ...
    'Units', 'normalized', ...
    'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize', 14, ...
    'Margin', 6, ...
    'Interpreter', 'none');

% Save second figure
set(fig2, 'PaperUnits', 'inches', 'PaperPosition', [0 0 12 5]);
exportgraphics(fig2, fullfile(dataDir, 'XZ_ErrorPlot.png'), 'Resolution', 300);
% Uncomment for PDF:
% exportgraphics(fig2, fullfile(dataDir, 'Figure2_ZDifferencePlot.pdf'), 'ContentType', 'vector');
