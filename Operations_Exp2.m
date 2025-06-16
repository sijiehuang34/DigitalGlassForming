%% Operations_Exp2.m
% This program is used to generate plots for papers/presentations.
clc; clear; close all;

one_px = 12.446e-6; % meters

dataDir = 'E:/SFF/Exp2_MeasuringTracks_DiffStartsStops/3_ConstantLaserPower_50W_0o5mms_fp_7mm/';
load(fullfile(dataDir, 'XZ_profile_visual_allFrames.mat'));

% Create figure for raw points and processed curves
figure('Color', 'w'); % good size for papers
hold on;

% Plot a representative gray marker for legend (before loop)
hGray = plot(nan, nan, 'o', 'MarkerEdgeColor', [0.6 0.6 0.6], ...
    'MarkerSize', 5, 'DisplayName', 'Raw data');

numFrames = height(results);
initial_shift = 15; % mm
increment = 0.5;     % mm per frame

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

% Plot median curve (bold, academic style)
plot(uniqueX, medianZ, 'r-', 'LineWidth', 2, 'DisplayName', 'Visual method');

% Formatting for academic presentation
xlabel('X (mm)');
ylabel('Z (mm)');
xlim([-20 20]);
ylim([0 1.5]);
xticks(-20:5:20);
yticks(0:0.25:1.5);
grid on;
box on;
set(gca, 'FontSize', 16);

legend('Location', 'northeast', 'FontSize', 14);


% % Shift X_2d so that its minimum is at 0
% X_2d = X_2d - min(X_2d);

% 
% % Plot confocal data
% plot(X_2d, Z_2d, 'o', 'MarkerSize', 6, 'LineWidth', 1.5, ...
%     'DisplayName', 'Confocal Profile');
% 
