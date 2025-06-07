%% FindBlob_ColoredBgSegmentation_Exp2.m
clc; clear; close all;

%% User Input
% Load image
expName = 'E:\SFF\Exp2_MeasuringTracks_DiffStartsStops\';
expCase = '1_ConstantLaserPower_40W_0o5mms_fp_7mm\';
expImg = 'images\Acquisition1000.jpeg';
fullImagePath = [expName, expCase, expImg];
image = imread(fullImagePath);

% %temp
% imshow(image);

% Define substrate region
substrate_start_row = 1350;
% Known widths of hypodermic tube & filament
w_hypodermic = 2e-3; %m        % w_hypodermic = 1.7e-3; %m
w_filament = 1.35e-3; %m        % w_filament = 1.4e-3; %m
known_hypodermic_depth = 400; %px (no need to be exact, but need to be sure, i.e. less is ok)
known_filament_start_depth = known_hypodermic_depth + 100; %px
one_px = 12.446e-6; % m

%% Color-based segmentation
[image_height, image_width, ~] = size(image);
substrate_mask = false(image_height, image_width);
substrate_mask(substrate_start_row:end, :) = true;
substrate_image = image;
substrate_image(repmat(~substrate_mask, [1, 1, 3])) = 0;
% %temp
% imshow(substrate_image);

%% Binary filter based on background color (excluding substrate)
% Define region size
region_height = substrate_start_row - 1; % Only consider above the substrate
region_width = 0.5 * (image_width - 200); % 1/2 of (width - thickness of tube)

% Convert to LAB color space
lab_image = rgb2lab(image);

% Extract 'a' and 'b' channels
a_channel = lab_image(:,:,2);
b_channel = lab_image(:,:,3);

% Define two background regions:
% Region 1: Top-left corner
region1_a = a_channel(1:region_height, 1:region_width);
region1_b = b_channel(1:region_height, 1:region_width);
% %temp
% region1_a = a_channel(1:1200, 1:300);
% region1_b = b_channel(1:1200, 1:300);

% Region 2: Top-right corner
region2_a = a_channel(1:region_height, image_width-region_width:image_width);
region2_b = b_channel(1:region_height, image_width-region_width:image_width);
% %temp
% region2_a = a_channel(1:1200, 700:1000);
% region2_b = b_channel(1:1200, 700:1000);

% Compute mean 'a' and 'b' values for both background regions
background_a = mean([mean2(region1_a), mean2(region2_a)]);
background_b = mean([mean2(region1_b), mean2(region2_b)]);

% Compute Euclidean distance from background color
distance = sqrt((a_channel - background_a).^2 + (b_channel - background_b).^2);

% Threshold to separate background and objects
threshold = 10;  % Adjust this value as needed
background_mask = distance < threshold;

% Exclude substrate region from background and object masks
background_mask(substrate_mask) = false;

% Create segmented images
background_image = image;
background_image(repmat(~background_mask, [1 1 3])) = 0;

object_mask = ~background_mask & ~substrate_mask;

object_image = image;
object_image(repmat(~object_mask, [1 1 3])) = 0;

% Create binary image (white object, black background)
bw = uint8(object_mask) * 255;
minSize = 1000;
bw = bwareaopen(bw, minSize);
bw = imfill(bw, 'holes'); % Fill holes

% Smooth the bw edges slightly to disconnect glass tip & substrate
se = strel('disk', 2); % Keep the radius small
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

% %temp
% subplot(1,4,1);
% imshow(image);
% subplot(1,4,2);
% imshow(bw);
% subplot(1,4,3);
% imshow(bw_smooth);

clear region_height region_width lab_image region1_a region1_b region2_a ...
      region2_b background_a background_b distance threshold background_mask...
      minSize a_channel b_channel se row right left

%% Find outline via bwboundaries

% Extract object boundary points
[Boundary, ~, ~] = bwboundaries(bw_smooth, 4, 'noholes');

% Sort in descending order of rows
[~, sortIdx] = sort(cellfun(@(x) size(x, 1), Boundary), 'descend');
Boundary = Boundary(sortIdx);

% %temp (disable if removed OI)

% Remove small noise: boundaries with fewer than 2000 points
Boundary = Boundary(cellfun(@(x) size(x,1) >= 2000, Boundary));
% Remove left edge noise
Boundary = Boundary(~cellfun(@(x) any(x(:,2) <= 250), Boundary));

boundaryPoints = Boundary{1};

% %temp
% subplot(1,4,4);
% % figure;
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

% Step 2: Fit the curve for the blob

% Extract depth of the 1st point (after width = 600) that goes over y1_mean
valid_idx = find(depth > known_filament_start_depth & width > y1_mean, 1, 'first');
if ~isempty(valid_idx)
    rough_blob_start = depth(valid_idx);
else
    rough_blob_start = NaN; 
end

idx_curve = (depth >= rough_blob_start) & (depth <= max(depth));
x_curve = depth(idx_curve);
y_curve = width(idx_curve);

% Fit a polynomial (2nd or 4th order)
p = polyfit(x_curve, y_curve, 2);
% p = polyfit(x_curve, y_curve, 4);

% Extend the parabola until it hits the x-axis (y = 0)
x_extended = linspace(min(depth), max(depth)*1.2, 1000); % Extend slightly beyond max(x) ~20%
y_fit = polyval(p, x_extended);

% Step 3: Find intersections with 2nd horizontal lines
tp = roots([p(1), p(2), p(3) - y2_mean]);
% tp = roots([p(1), p(2), p(3), p(4), p(5) - y2_mean]);

tp = tp(imag(tp) == 0);  % Only real solutions
tp_ref = min(tp);        % Choose leftmost intersection

clear idx1 idx_curve valid_idx

%% Filter workzone points
% Keep only the boundary points where y > blob_tp
printedPart_outline = boundaryPoints(boundaryPoints(:,1) > tp_ref, :);

x_printedPart = printedPart_outline(:,2);
y_printedPart = printedPart_outline(:,1);

y_opening = min(y_printedPart);
x_opening = x_printedPart(y_printedPart == y_opening);
x_opening_start = min(x_opening);
x_opening_end = max(x_opening);

x_range = x_opening_start:x_opening_end;
y_range = y_opening * ones(size(x_range));

clear printedPart_outline y_opening x_opening x_opening_start x_opening_end

%% Display results
width = width * one_px * 1000; %mm
depth = depth * one_px * 1000; %mm
y1_mean = y1_mean * one_px * 1000; %mm
y2_mean = y2_mean * one_px * 1000; %mm
x_extended = x_extended * one_px * 1000; %mm
y_fit = y_fit * one_px * 1000; %mm
tp_ref = tp_ref * one_px * 1000; %mm

figure;
plot(depth, width, 'k.', 'MarkerSize', 10, 'DisplayName', 'width-depth data'); 
hold on;
yline(y1_mean, 'b:', 'LineWidth', 2, 'DisplayName', 'hypodermic tube');  % First horizontal line
yline(y2_mean, 'b--', 'LineWidth', 2, 'DisplayName', 'cold filament'); % Second horizontal line
plot(x_extended, y_fit, 'g-', 'LineWidth', 2, 'DisplayName', 'work zone'); % Ellipse arc
plot(tp_ref, y2_mean, 'ro', 'MarkerSize', 15, 'LineWidth', 3, 'DisplayName', 'work zone start point');
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
subplot(2,4,1);
imshow(image);
title("Original");

subplot(2,4,2);
imshow(background_image);
title("Background");

subplot(2,4,3);
imshow(object_image);
title("OI Noisy");

subplot(2,4,4);
imshow(bw);
title("Binarize & Denoise");

subplot(2,4,5);
imshow(image);
hold on;
plot(x_all, y_all, 'r.', "MarkerSize", 3);
hold off;
title("OI ONLY");

subplot(2,4,6);
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

subplot(2,4,7);
imshow(image);
hold on;
plot(x_printedPart, y_printedPart, 'r.', "MarkerSize", 3);
plot(x_range, y_range, 'r.', 'MarkerSize', 3);
hold off;
title("Printed Part");

clear i x_range y_range