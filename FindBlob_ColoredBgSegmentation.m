%% FindBlob_ColoredBgSegmentation.m
clc; clear; close all;

%% Color-based segmentation
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Load image (working frames: 1760-3756, except those connecting to the left)

image = imread('E:\3_PowerRamp_1W\3_PowerRamp_1W\images\Acquisition1760.jpeg');
[image_height, image_width, ~] = size(image);

% Define substrate region (bottom 1/6 of the image)
finetune_val = 35;
substrate_start_row = floor(image_height * (5/6)) + 1 + finetune_val;
substrate_mask = false(image_height, image_width);
substrate_mask(substrate_start_row:end, :) = true;

substrate_image = image;
substrate_image(repmat(~substrate_mask, [1, 1, 3])) = 0;

% Binary filter based on background color (excluding substrate)
% Define region size
region_height = substrate_start_row - 1; % Only consider above the substrate
region_width = 400; % 1/2 of (width - thickness of tube)

% Convert to LAB color space
lab_image = rgb2lab(image);

% Extract 'a' and 'b' channels
a_channel = lab_image(:,:,2);
b_channel = lab_image(:,:,3);

% Define two background regions:
% Region 1: Top-left corner
region1_a = a_channel(1:region_height, 1:region_width);
region1_b = b_channel(1:region_height, 1:region_width);

% Region 2: Top-right corner
region2_a = a_channel(1:region_height, end-region_width:end);
region2_b = b_channel(1:region_height, end-region_width:end);

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
se = strel('disk', 3); % Keep the radius small
bw_smooth = imopen(bw, se);  % Removes small protrusions
bw_smooth = imclose(bw_smooth, se); % Fills small gaps
bw_smooth = bwareaopen(bw_smooth, minSize);

%temp
subplot(1,4,1);
imshow(image);
subplot(1,4,2);
imshow(bw);
subplot(1,4,3);
imshow(bw_smooth);

clear region_height region_width lab_image region1_a region1_b region2_a ...
      region2_b background_a background_b distance threshold background_mask...
      minSize a_channel b_channel se

%% Find outline via bwboundaries

% Extract object boundary points
[Boundary, Label, NumObj] = bwboundaries(bw_smooth, 4, 'noholes');

% Sort in descending order of rows
[~, sortIdx] = sort(cellfun(@(x) size(x, 1), Boundary), 'descend');
Boundary = Boundary(sortIdx);

% Remove small noise: boundaries with fewer than 2000 points
Boundary = Boundary(cellfun(@(x) size(x,1) >= 2000, Boundary));

% Remove left edge noise: boundaries that contain x = 1
Boundary = Boundary(~cellfun(@(x) any(x(:,2) == 1), Boundary));
% Remove right edge noise: boundaries that contain x = image_width
Boundary = Boundary(~cellfun(@(x) any(x(:,2) == image_width), Boundary));

% boundaryPoints = Boundary{1};

%temp
subplot(1,4,4);
% figure;
imshow(image);
hold on;
for i = 1:length(Boundary)
    boundaryPoints = Boundary{i};
    plot(boundaryPoints(:,2), boundaryPoints(:,1), 'r-', 'LineWidth', 1.5);
end
hold off;

clear sortIdx perimeterLengths sortedIndices topIndices i

%% Plot width vs. depth
% Sort boundaryPoints by ascending order of y value
[sortedY, sortIdx] = sort(boundaryPoints(:,1), 'ascend'); % Sort by y-values
sortedX = boundaryPoints(sortIdx, 2); % Reorder x-values accordingly
boundaryPoints = [sortedY, sortedX]; % Reconstruct sorted boundaryPoints
x_all = boundaryPoints(:,2);
y_all = boundaryPoints(:,1);

% Find unique y values and count them
uniqueY = unique(boundaryPoints(:,1), 'sorted'); % Unique y-values in ascending order
n = length(uniqueY); % Number of unique y-values

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

% Find blob diameter
blob_D = max(wid_vs_dep(:,2));

clear sortedY sortIdx sortedX uniqueY n yVal xVals i

%% Fit curves & find intersections (2nd order - parabola)
depth = wid_vs_dep(:,1);
width = wid_vs_dep(:,2);

% Step 1: Fit the horizontal lines (between x = 0 and 500 <-- constant)
idx1 = (depth >= 0) & (depth <= 500);
y1_mean = mean(width(idx1));   % First horizontal line y-value
% Known widths of hypodermic tube & filament
w_hypodermic = 1.7e-3; %m
w_filament = 1.4e-3; %m
% Convert to px
one_px = 12.446e-6; % um
w_hypodermic = w_hypodermic / one_px; %px
w_filament = w_filament / one_px; %px
w_diff = w_hypodermic - w_filament;
y2_mean = y1_mean - w_diff;    % Second horizontal line y-value

% Step 2: Fit the curve for the blob

% Extract depth of the 1st point (after width = 600) that goes over y1_mean
valid_idx = find(depth > 600 & width > y1_mean, 1, 'first');
if ~isempty(valid_idx)
    rough_blob_start = depth(valid_idx);
else
    rough_blob_start = NaN; 
end

idx_curve = (depth >= rough_blob_start) & (depth <= max(depth));
x_curve = depth(idx_curve);
y_curve = width(idx_curve);

% Fit a polynomial (4th order)
p = polyfit(x_curve, y_curve, 2);
% p = polyfit(x_curve, y_curve, 4);
% p = polyfit(x_curve, y_curve, 6);
% p = polyfit(x_curve, y_curve, 8);

% Extend the parabola until it hits the x-axis (y = 0)
x_extended = linspace(min(depth), max(depth)*1.2, 1000); % Extend slightly beyond max(x) ~20%
y_fit = polyval(p, x_extended);

% Step 3: Find intersections with 2nd horizontal lines
tp = roots([p(1), p(2), p(3) - y2_mean]);
% tp = roots([p(1), p(2), p(3), p(4), p(5) - y2_mean]);
% tp = roots([p(1), p(2), p(3), p(4), p(5), p(6), p(7) - y2_mean]);
% tp = roots([p(1), p(2), p(3), p(4), p(5), p(6), p(7), p(8), p(9) - y2_mean]);

tp = tp(imag(tp) == 0);  % Only real solutions
tp_ref = min(tp);        % Choose leftmost intersection

clear idx1 idx_curve x_curve y_curve tp p wid_vs_dep valid_idx

%% Find height of blob above substrate
blob_substrate_dist = substrate_start_row - max(depth); %px
blob_substrate_dist = blob_substrate_dist * one_px * 1000; %mm

%% Filter blob points
% Keep only the boundary points where y > blob_tp
blob_outline = boundaryPoints(boundaryPoints(:,1) > tp_ref, :);

x_blob = blob_outline(:,2);
y_blob = blob_outline(:,1);

y_opening = min(y_blob);
x_opening = x_blob(y_blob == y_opening);
x_opening_start = min(x_opening);
x_opening_end = max(x_opening);

x_range = x_opening_start:x_opening_end;
y_range = y_opening * ones(size(x_range));

clear blob_outline y_opening x_opening x_opening_start x_opening_end

%% Display results
width = width * one_px * 1000; %mm
depth = depth * one_px * 1000; %mm
y1_mean = y1_mean * one_px * 1000; %mm
y2_mean = y2_mean * one_px * 1000; %mm
x_extended = x_extended * one_px * 1000; %mm
y_fit = y_fit * one_px * 1000; %mm
tp_ref = tp_ref * one_px * 1000; %mm

figure;
plot(depth, width, 'k.', 'MarkerSize', 5); 
hold on;
plot(x_extended, y_fit, 'g-', 'LineWidth', 2); % Ellipse arc
yline(y1_mean, 'b-', 'LineWidth', 1.5);       % First horizontal line
yline(y2_mean, 'b-', 'LineWidth', 1.5);       % Second horizontal line
plot(tp_ref, y2_mean, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
xlabel('Depth (mm)');
ylabel('Width (mm)');
grid on;
xlim([0 20]);
ylim([0 3]);
hold off;

figure;
subplot(2,3,1);
imshow(image);
title("Original Image");

subplot(2,3,2);
imshow(background_image);
title("Background Only");

subplot(2,3,3);
imshow(object_image);
title("Glass Only");

subplot(2,3,4);
imshow(bw);
title("Denoised Glass Only");

subplot(2,3,5);
imshow(image);
hold on;
plot(x_all, y_all, 'r.', "MarkerSize", 3);
hold off;
title("Entire object outlined");

subplot(2,3,6);
imshow(image);
hold on;
plot(x_blob, y_blob, 'r.', "MarkerSize", 3);
plot(x_range, y_range, 'r.', 'MarkerSize', 3);
hold off;
title("Blob outlined");
