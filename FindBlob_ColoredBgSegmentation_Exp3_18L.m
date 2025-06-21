%% FindBlob_ColoredBgSegmentation_Exp3_18L.m
% This program is used for testing methods on individual frame.
clc; clear; close all;

%% Height Profile %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% User Input
% Load image
expName = 'E:\SFF\Exp3_JoiningTracks\';
expCase = '18L_0.5mmsPRINT_40W_20mm_4mmAboveFocus_BlobTip_0mmOverlap_1mmsRETRACT\';
expImg = 'images\Acquisition620.jpeg';
saveDir = 'E:\SFF\Exp3_JoiningTracks\18L_0.5mmsPRINT_40W_20mm_4mmAboveFocus_BlobTip_0mmOverlap_1mmsRETRACT\';
fullImagePath = [expName, expCase, expImg];
image = imread(fullImagePath);

% Define substrate region
substrate_start_row = 1235;
w_hypodermic = 2e-3; % m
w_filament = 1.5e-3; % m
known_hypodermic_depth = 190; % px
known_filament_start_depth = known_hypodermic_depth + 100; % px
one_px = 12.446e-6; % m

%% Color-based segmentation
[image_height, image_width, ~] = size(image);
substrate_mask = false(image_height, image_width);
substrate_mask(substrate_start_row:end, :) = true;

%% Binary filter based on background color (excluding substrate)
region_height = substrate_start_row - 1;
region_width = 0.5 * (image_width - 200);

% Convert to LAB color space
lab_image = rgb2lab(image);
a_channel = lab_image(:,:,2);
b_channel = lab_image(:,:,3);

% Define background regions of distinct colors
% Region 1 (lower part)
region1_a = a_channel(700:region_height, 1:region_width);
region1_b = b_channel(700:region_height, 1:region_width);
region1_a2 = a_channel(800:region_height-100, image_width-region_width:image_width);
region1_b2 = b_channel(800:region_height-100, image_width-region_width:image_width);

% Region 2 (upper part)
region2_a = a_channel(1:500, 1:region_width);
region2_b = b_channel(1:500, 1:region_width);
region2_a2 = a_channel(1:500, image_width-region_width:image_width);
region2_b2 = b_channel(1:500, image_width-region_width:image_width);

% Mean 'a' and 'b' values for each region
% Region 1 (lower part)
bg1_a = mean([mean2(region1_a), mean2(region1_a2)]); 
bg1_b = mean([mean2(region1_b), mean2(region1_b2)]);
% Region 2 (upper part)
bg2_a = mean([mean2(region2_a), mean2(region2_a2)]); 
bg2_b = mean([mean2(region2_b), mean2(region2_b2)]);

% Compute Euclidean distance from both background colors
distance_bg1 = sqrt((a_channel - bg1_a).^2 + (b_channel - bg1_b).^2);
distance_bg2 = sqrt((a_channel - bg2_a).^2 + (b_channel - bg2_b).^2);

% Define threshold
threshold = 12;

% Create masks
bg_mask_1 = distance_bg1 < threshold; 
bg_mask_2 = distance_bg2 < threshold; 

% Remove substrate region from background
bg_mask_1(substrate_mask) = false;
bg_mask_2(substrate_mask) = false;

% Object mask is whatever is not in either background or substrate
object_mask = ~bg_mask_1 & ~bg_mask_2 & ~substrate_mask;

% Create segmented images
bg_image_1 = image;
bg_image_1(repmat(~bg_mask_1, [1 1 3])) = 0;

bg_image_2 = image;
bg_image_2(repmat(~bg_mask_2, [1 1 3])) = 0;

object_image = image;
object_image(repmat(~object_mask, [1 1 3])) = 0;

% Create binary image (white object, black background)
bw = uint8(object_mask) * 255;
minSize = 2000;
bw = bwareaopen(bw, minSize);
bw = imfill(bw, 'holes'); % Fill holes

% Smooth the bw edges slightly to disconnect glass tip & substrate
se = strel('disk', 4); % Keep the radius small
bw_smooth = imopen(bw, se);  % Removes small protrusions
bw_smooth = imclose(bw_smooth, se); % Fills small gaps
bw_smooth = bwareaopen(bw_smooth, minSize);

% Loop through the first 100 rows and fill between first white pixels in each row
for row = 1:min(100, size(bw_smooth, 1))
    thisRow = bw_smooth(row, :);
    % Find indices of white pixels (logical 1s)
    whiteIdx = find(thisRow);
    if ~isempty(whiteIdx)
        left = whiteIdx(1);               % First white pixel from left
        right = whiteIdx(end);            % First white pixel from right
        % Set all pixels between left and right to white
        bw_smooth(row, left:right) = 1;
    end
end

boundary = [];  % Initialize an empty array to hold boundary coordinates

for col = 1:image_width
    col_data = bw_smooth(:, col);
    row_idx = find(col_data, 1, 'first');  % First white pixel from the top
    if ~isempty(row_idx)
        boundary = [boundary; col, row_idx];  % Store as (x, y)
    end
end

% Sort boundary by descending x (column)
[sortedX, sortIdx] = sort(boundary(:,1), 'descend');
sortedY = boundary(sortIdx, 2);

% Define incontinuity threshold
incontThreshold = 300; % pixels

% Compute differences in Y
dy = diff(sortedY);

% Find all jump indices where the Y-jump exceeds the threshold
jumpIdxAll = find(abs(dy) > incontThreshold);

if length(jumpIdxAll) >= 2
    % Get the first two incontinuity points
    idx1 = jumpIdxAll(1);
    idx2 = jumpIdxAll(2);

    % Get corresponding X values (sorted from right to left)
    x1 = sortedX(idx1);
    x2 = sortedX(idx2);

    % Get the boundary indices in between those Xs
    x_min_jump = min(x1, x2);
    x_max_jump = max(x1, x2);
else
    x_min_jump = -Inf;
    x_max_jump = -Inf;
end

% Additional exclusion range (in pixels)
x_min_custom = 400;
x_max_custom = 650;

% Combine both exclusion ranges
excludeMask = (boundary(:,1) > x_min_jump & boundary(:,1) < x_max_jump) | ...
              (boundary(:,1) > x_min_custom & boundary(:,1) < x_max_custom);
% excludeMask = (boundary(:,1) < 250) | ...
%               (boundary(:,1) > x_min_custom & boundary(:,1) < x_max_custom);

% Keep only the points outside both exclusion zones
boundary = boundary(~excludeMask, :);

% Height and length profile
Z_vis = (substrate_start_row - boundary(:,2)); % px
X_vis = boundary(:,1); % px

% Pixel to metric conversion
Z_vis = Z_vis * one_px * 1000; % mm
X_vis = X_vis * one_px * 1000; % mm

%% Apply DBSCAN to remove noise (assumes 2 clusters of interest)
dataPoints = [X_vis, Z_vis];
epsilon = 0.30;      % Distance threshold in mm (tune as needed)
minpts = 5;          % Minimum number of points per cluster
labels = dbscan(dataPoints, epsilon, minpts);

% Keep only the two largest clusters
validClusters = labels > 0;
uniqueClusters = unique(labels(validClusters));
clusterSizes = arrayfun(@(c) sum(labels == c), uniqueClusters);

[~, sortedIdx] = sort(clusterSizes, 'descend');
topTwoLabels = uniqueClusters(sortedIdx(1:min(2, numel(sortedIdx))));

% Create mask for points in the two largest clusters
clusterMask = ismember(labels, topTwoLabels);

% Filtered coordinates
X_vis = X_vis(clusterMask);
Z_vis = Z_vis(clusterMask);

%% Figure
figure;
subplot(2,4,1);
imshow(image), title('OG');
subplot(2,4,2);
imshow(bg_image_1), title('BG 1');
subplot(2,4,3);
imshow(bg_image_2), title('BG 2');
subplot(2,4,4);
imshow(object_image), title('Object');
subplot(2,4,5);
imshow(bw), title('Binary');
subplot(2,4,6);
imshow(bw_smooth), title('Denoised');
subplot(2,4,7);
imshow(bw_smooth), title('Height Profile');
hold on;
plot(boundary(:,1), boundary(:,2), 'r.', 'MarkerSize', 5);
hold off;

figure;
hold on;
plot(X_vis, Z_vis, 'o');
xlabel('X (mm)'); xlim([0 14]);
ylabel('Z (mm)'); ylim([0 2]);
grid on;
hold off;

% Save the data
save(fullfile(saveDir, 'XZ_profile_visual.mat'), 'X_vis', 'Z_vis');

clear region_height region_width lab_image region1_a region1_b region2_a ...
      region2_b background_a background_b distance threshold background_mask...
      minSize a_channel b_channel se row right left

%% Work zone data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find outline via bwboundaries

% Extract object boundary points
[Boundary, ~, ~] = bwboundaries(bw_smooth, 4, 'noholes');

% Filter boundaries to only keep those with at least one Y > 1100
Boundary = Boundary(cellfun(@(b) any(b(:,1) > 1100), Boundary));

% Sort in descending order of number of rows
[~, sortIdx] = sort(cellfun(@(x) size(x, 1), Boundary), 'descend');
Boundary = Boundary(sortIdx);

boundaryPoints = Boundary{1};

% %temp
% figure;
% imshow(image);
% hold on;
% for i = 1:length(Boundary)
%     boundaryPoints = Boundary{i};
%     plot(boundaryPoints(:,2), boundaryPoints(:,1), 'r-', 'LineWidth', 1.5);
% end
% hold off;

clear sortIdx perimeterLengths sortedIndices topIndices i

%% Identify workzone 
% Sort boundaryPoints by ascending order of y value
[sortedY, sortIdx] = sort(boundaryPoints(:,1), 'ascend');
sortedX = boundaryPoints(sortIdx, 2);
boundaryPoints = [sortedY, sortedX];
x_all = boundaryPoints(:,2);
y_all = boundaryPoints(:,1);
% Find unique y values and count them
uniqueY = unique(boundaryPoints(:,1), 'sorted'); % Unique y-values in ascending order
n = length(uniqueY); % Number of unique y-values

% Initialize workzonePoints
workzonePoints = [];

% Loop through each unique y value
for i = 1:length(uniqueY)
    thisY = uniqueY(i);
    x_values = x_all(y_all == thisY);  % Find x-values for current y
    if isempty(x_values)
        continue;
    end
    x_min = min(x_values);
    x_max = max(x_values);

    % Ensure x range is within bounds
    x_min = max(1, round(x_min));
    x_max = min(image_width, round(x_max));

    % Loop through pixels in this row between x_min and x_max
    for x = x_min:x_max
        pixel = squeeze(object_image(thisY, x, :));  % RGB pixel values
        if all(pixel >= 240)  % Check if pixel is close to white
            workzonePoints = [workzonePoints; x, thisY];  % Store (x, y)
        end
    end
end

% Remove outlier points from workzonePoints using DBSCAN
if ~isempty(workzonePoints)
    % Normalize coordinates for more consistent clustering
    coords = double(workzonePoints); % (x, y)
    
    % Apply DBSCAN clustering
    epsilon = 5;  % max distance between points in the same cluster (in pixels)
    minpts = 50;   % minimum number of points to form a cluster
    clusterLabels = dbscan(coords, epsilon, minpts);

    % Find the largest cluster (excluding label -1 for noise)
    validLabels = unique(clusterLabels(clusterLabels > 0));
    if ~isempty(validLabels)
        clusterSizes = arrayfun(@(x) sum(clusterLabels == x), validLabels);
        [~, largestClusterIdx] = max(clusterSizes);
        keepLabel = validLabels(largestClusterIdx);
        workzonePoints = workzonePoints(clusterLabels == keepLabel, :);
    else
        warning('No dense clusters found. All points may be outliers.');
        workzonePoints = [];
    end
end

% Create binary mask from workzonePoints for outlining
workzoneMask = false(size(image,1), size(image,2)); % same size as image
for k = 1:size(workzonePoints,1)
    y = workzonePoints(k,2);  % row
    x = workzonePoints(k,1);  % column
    if x > 0 && y > 0 && y <= size(workzoneMask,1) && x <= size(workzoneMask,2)
        workzoneMask(y, x) = true;
    end
end

wz_boundaries = bwboundaries(workzoneMask, 'noholes');

% Calculate workzone centroid
if ~isempty(workzonePoints)
    centroidX = mean(workzonePoints(:,1));
    centroidY = mean(workzonePoints(:,2));
    workzoneCentroid = [centroidX, centroidY];  % [x, y] format
else
    workzoneCentroid = [NaN, NaN];
    warning('Workzone centroid could not be computed: no valid points found.');
end

clear clusterLabels clusterSizes coords i k keepLabel thisY x_values x_max...
      x_min x y workzoneMask whiteIdx thisRow validLabels pixel minpts

%% Identify printed part - plot width vs. depth
% Initialize y_vs_width matrix (n x 2)
wid_vs_dep = zeros(n, 2);
wid_vs_dep(:,1) = uniqueY; % First column: unique y-values

% Compute width for each y value
for i = 1:n
    yVal = uniqueY(i); % Current y-value
    xVals = boundaryPoints(boundaryPoints(:,1) == yVal, 2); % Corresponding x-values
    width = abs(max(xVals) - min(xVals)); % Compute absolute width
    wid_vs_dep(i,2) = width; % Store in second column
end

clear sortedY sortIdx sortedX uniqueY n yVal xVals i

depth = wid_vs_dep(:,1); % px
width = wid_vs_dep(:,2); % px

% %temp
% figure;
% plot(depth, width, 'k.', 'MarkerSize', 10); 

% Step 1: Fit the horizontal lines (between x = 0 and known_hypodermic_depth: constant)
idx1 = (depth >= 0) & (depth <= known_hypodermic_depth);
y1_mean = mean(width(idx1));   % First horizontal line y-value

% Convert known widths from mm to px
w_hypodermic = w_hypodermic / one_px; %px
w_filament = w_filament / one_px; %px
w_diff = w_hypodermic - w_filament;
y2_mean = y1_mean - w_diff;    % Second horizontal line y-value

% % Step 2: Fit the curve for the blob
% 
% % Extract depth of the 1st point (after width = 600) that goes over y1_mean
% valid_idx = find(depth > known_filament_start_depth & width > y1_mean, 1, 'first');
% if ~isempty(valid_idx)
%     rough_blob_start = depth(valid_idx);
% else
%     rough_blob_start = NaN; 
% end
% 
% idx_curve = (depth >= rough_blob_start) & (depth <= max(depth));
% x_curve = depth(idx_curve);
% y_curve = width(idx_curve);
% 
% % Fit a polynomial (2nd or 4th order)
% p = polyfit(x_curve, y_curve, 2);
% % p = polyfit(x_curve, y_curve, 4);
% 
% % Extend the parabola until it hits the x-axis (y = 0)
% x_extended = linspace(min(depth), max(depth)*1.2, 1000); % Extend slightly beyond max(x) ~20%
% y_fit = polyval(p, x_extended);
% 
% % Step 3: Find intersections with 2nd horizontal lines
% tp = roots([p(1), p(2), p(3) - y2_mean]);
% % tp = roots([p(1), p(2), p(3), p(4), p(5) - y2_mean]);
% 
% tp = tp(imag(tp) == 0);  % Only real solutions
% tp_ref = min(tp);        % Choose leftmost intersection
% 
% clear idx1 idx_curve valid_idx
% 
% %% Filter workzone points
% % Keep only the boundary points where y > blob_tp
% printedPart_outline = boundaryPoints(boundaryPoints(:,1) > tp_ref, :);
% 
% x_printedPart = printedPart_outline(:,2);
% y_printedPart = printedPart_outline(:,1);
% 
% y_opening = min(y_printedPart);
% x_opening = x_printedPart(y_printedPart == y_opening);
% x_opening_start = min(x_opening);
% x_opening_end = max(x_opening);
% 
% x_range = x_opening_start:x_opening_end;
% y_range = y_opening * ones(size(x_range));
% 
% clear printedPart_outline y_opening x_opening x_opening_start x_opening_end
% 
%% Display results
width = width * one_px * 1000; %mm
depth = depth * one_px * 1000; %mm
y1_mean = y1_mean * one_px * 1000; %mm
y2_mean = y2_mean * one_px * 1000; %mm
% x_extended = x_extended * one_px * 1000; %mm
% y_fit = y_fit * one_px * 1000; %mm
% tp_ref = tp_ref * one_px * 1000; %mm

figure;
plot(depth, width, 'k.', 'MarkerSize', 10, 'DisplayName', 'width-depth data'); 
hold on;
yline(y1_mean, 'b:', 'LineWidth', 2, 'DisplayName', 'hypodermic tube');  % First horizontal line
yline(y2_mean, 'b--', 'LineWidth', 2, 'DisplayName', 'cold filament'); % Second horizontal line
% plot(x_extended, y_fit, 'g-', 'LineWidth', 2, 'DisplayName', 'work zone'); % Ellipse arc
% plot(tp_ref, y2_mean, 'ro', 'MarkerSize', 15, 'LineWidth', 3, 'DisplayName', 'work zone start point');
legend("Location", "best");

% xlabel('Depth (mm)', 'FontSize', 40);
% ylabel('Width (mm)', 'FontSize', 40);
xlabel('Depth (mm)');
ylabel('Width (mm)');

grid on;
% xlim([0 15]);
ylim([0 8]);

% ax = gca;
% % ax.FontSize = 40;       % Axis tick label size
% ax.LineWidth = 1.5;       % Thicken the outer box
% box on;                 % Make sure the box is visible
set(gca, "FontSize", 14);
hold off;

% exportgraphics(gcf, 'myFigure.png', 'Resolution', 300);
% % Save with DPI control
% dataDir = 'D:/ND Lab/'; % New dataset location
% print(fig, fullfile(dataDir, 'Plot'), '-djpeg', '-r300'); % 300 DPI JPEG
% savefig(fig, fullfile(dataDir, 'Plot.fig'));

figure;
subplot(1,3,1);
imshow(image);
hold on;
plot(x_all, y_all, 'r.', "MarkerSize", 3);
hold off;
title("OI ONLY");

subplot(1,3,2);
imshow(image);
hold on;
for i = 1:length(wz_boundaries)
    wz_boundary = wz_boundaries{i};
    plot(wz_boundary(:,2), wz_boundary(:,1), 'r.', "MarkerSize", 3);
end
if ~any(isnan(workzoneCentroid))
    plot(workzoneCentroid(1), workzoneCentroid(2), 'rx', 'MarkerSize', 5);
end
hold off;
title("Workzone with Centroid");

% subplot(1,3,3);
% imshow(image);
% hold on;
% plot(x_printedPart, y_printedPart, 'r.', "MarkerSize", 3);
% plot(x_range, y_range, 'r.', 'MarkerSize', 3);
% hold off;
% title("Printed Part");

clear i x_range y_range

%% Final overlay on original image
figure;
imshow(image); hold on;

% Convert X_vis and Z_vis back to pixel coordinates
X_px = X_vis / (one_px * 1000);
Y_px = substrate_start_row - (Z_vis / (one_px * 1000));

% Get cluster labels for filtered points
finalLabels = labels(clusterMask);  % Get labels for plotted points
colors = ['b', 'r'];                % Track 1 = blue, Track 2 = red
styles = {'--', '-.'};              % Line styles: dashed and dotted
clusterPoints = [X_px, Y_px];

% Plot height profile by cluster
for i = 1:length(topTwoLabels)
    clusterID = topTwoLabels(i);
    mask = finalLabels == clusterID;
    plot(clusterPoints(mask,1), clusterPoints(mask,2), ...
         'LineStyle', styles{i}, 'Color', colors(i), ...
         'LineWidth', 1.5);
end

% Plot centroid
if ~any(isnan(workzoneCentroid))
    plot(workzoneCentroid(1), workzoneCentroid(2), 'kx', ...
         'MarkerSize', 10, 'LineWidth', 2);
end

% Plot work zone outline (label only once)
labelPlotted = false;
for i = 1:length(wz_boundaries)
    wz_boundary = wz_boundaries{i};
    if ~labelPlotted
        plot(wz_boundary(:,2), wz_boundary(:,1), 'k-', 'LineWidth', 1.5);
        labelPlotted = true;
    else
        plot(wz_boundary(:,2), wz_boundary(:,1), 'k-', 'LineWidth', 1.5);
    end
end

legend('First track', 'Second track', 'Centroid', 'Work zone', 'Location', 'northwest');
hold off;
set(gca, 'FontSize', 15);
