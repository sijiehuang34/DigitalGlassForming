%% FindBlob_ColoredBgSegmentation.m
clc; clear; close all;

%% User Input
% Load image
image = imread('E:\SFF\Exp1_ConstantPowerExp\3_30W_2mmAboveFocus_FlatCut\images\Acquisition200.jpeg');
% Define substrate region
substrate_start_row = 1750;
% Known widths of hypodermic tube & filament
w_hypodermic = 2e-3; %m        % w_hypodermic = 1.7e-3; %m
w_filament = 1.4e-3; %m        % w_filament = 1.4e-3; %m
known_hypodermic_depth = 250; %px (no need to be exact, but need to be sure, i.e. less is ok)
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
se = strel('disk', 3); % Keep the radius small
bw_smooth = imopen(bw, se);  % Removes small protrusions
bw_smooth = imclose(bw_smooth, se); % Fills small gaps
bw_smooth = bwareaopen(bw_smooth, minSize);

% Loop through the first 50 rows and fill between first white pixels in each row
for row = 1:min(50, size(bw_smooth, 1))
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

% %temp (disabled)

% Remove small noise: boundaries with fewer than 2000 points
Boundary = Boundary(cellfun(@(x) size(x,1) >= 2000, Boundary));
% Remove left edge noise: boundaries that contain x = 1
Boundary = Boundary(~cellfun(@(x) any(x(:,2) == 1), Boundary));
% % Remove right edge noise: boundaries that contain x = image_width
% Boundary = Boundary(~cellfun(@(x) any(x(:,2) == image_width), Boundary));

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

clear sortedY sortIdx sortedX uniqueY n yVal xVals i

%% Fit curves & find intersections (2nd order - parabola)
depth = wid_vs_dep(:,1); % px
width = wid_vs_dep(:,2); % px

%temp
figure;
plot(depth, width, 'k.', 'MarkerSize', 10); 

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

% Step 4: Calculate blob area (area under the parabola and above y = 0)

% Find real roots of the parabola (where it intersects y = 0)
x_roots = roots(p);
x_roots = x_roots(imag(x_roots) == 0); % only real
if length(x_roots) < 2
    error('Parabola does not intersect y=0 twice — cannot compute blob area.');
end

% Sort roots to define integration limits
x1 = min(x_roots);
x2 = max(x_roots);

% Define fine x range for numerical integration
x_blob = linspace(x1, x2, 1000);
y_blob = polyval(p, x_blob);

% Only keep area where y > 0 (in case part of the parabola dips below)
y_blob(y_blob < 0) = 0;

% Numerically integrate (area in px^2)
area_px2 = trapz(x_blob, y_blob);

% Convert to mm²
area_mm2 = area_px2 * one_px^2 * 1e6; % (12.446e-6 m)^2 → mm²

fprintf('Work zone area: %.4f mm²\n', area_mm2);

clear idx1 idx_curve x_curve y_curve tp wid_vs_dep valid_idx x_roots x1 x2 x_blob y_blob

%% Find height of blob above substrate
blob_substrate_dist = substrate_start_row - max(depth); %px
blob_substrate_dist = blob_substrate_dist * one_px * 1000; %mm
fprintf('Work zone height: %.4f mm\n', blob_substrate_dist);

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

% fig = figure;
% % Set figure size in inches
% width_inch = 12; % Width in inches
% height_inch = 9; % Height in inches
% set(fig, 'Units', 'inches', 'Position', [1, 1, width_inch, height_inch]);

figure;
plot(depth, width, 'k.', 'MarkerSize', 10); 
hold on;
plot(x_extended, y_fit, 'g-', 'LineWidth', 2); % Ellipse arc
yline(y1_mean, 'b:', 'LineWidth', 2);         % First horizontal line
yline(y2_mean, 'b--', 'LineWidth', 2);        % Second horizontal line
plot(tp_ref, y2_mean, 'ro', 'MarkerSize', 15, 'LineWidth', 3);

% xlabel('Depth (mm)', 'FontSize', 40);
% ylabel('Width (mm)', 'FontSize', 40);
xlabel('Depth (mm)');
ylabel('Width (mm)');

grid on;
% xlim([0 15]);
ylim([0 5]);

ax = gca;
% ax.FontSize = 40;       % Axis tick label size
ax.LineWidth = 1.5;       % Thicken the outer box
box on;                 % Make sure the box is visible

% legend("width-depth data", "work zone", "hypodermic tube", "cold filament", "work zone start point","FontSize", 20)
legend("width-depth data", "work zone", "hypodermic tube", "cold filament", "work zone start point", "Location", "best")

hold off;
% exportgraphics(gcf, 'myFigure.png', 'Resolution', 300);
% % Save with DPI control
% dataDir = 'D:/ND Lab/'; % New dataset location
% print(fig, fullfile(dataDir, 'Plot'), '-djpeg', '-r300'); % 300 DPI JPEG
% savefig(fig, fullfile(dataDir, 'Plot.fig'));

figure;
subplot(2,3,1);
imshow(image);
title("Original");

subplot(2,3,2);
imshow(background_image);
title("Background");

subplot(2,3,3);
imshow(object_image);
title("Object of Interest (OI)");

subplot(2,3,4);
imshow(bw);
title("Binarize & Denoise");

subplot(2,3,5);
imshow(image);
hold on;
plot(x_all, y_all, 'r.', "MarkerSize", 3);
hold off;
title("OI Outlined");

subplot(2,3,6);
imshow(image);
hold on;
plot(x_blob, y_blob, 'r.', "MarkerSize", 3);
plot(x_range, y_range, 'r.', 'MarkerSize', 3);
hold off;
title("Work Zone Outlined");

%% Other plot arrangements
% figure;
% plot(depth, width, 'k.', 'MarkerSize', 5); 
% hold on;
% plot(x_extended, y_fit, 'g-', 'LineWidth', 1.25); % Ellipse arc
% yline(y1_mean, 'b--', 'LineWidth', 1.25);       % First horizontal line
% yline(y2_mean, 'b-.', 'LineWidth', 1.25);       % Second horizontal line
% plot(tp_ref, y2_mean, 'ro', 'MarkerSize', 8, 'LineWidth', 1.25);
% xlabel('Depth (mm)');
% ylabel('Width (mm)');
% grid on;
% xlim([0 20]);
% ylim([0 3]);
% legend("width-depth data", "work zone", "hypodermic tube", "cold filament", "work zone starting point")
% set(gca, "FontSize", 12);
% hold off;

% figure;
% plot(depth, width, 'k.', 'MarkerSize', 10); 
% hold on;
% plot(x_extended, y_fit, 'g-', 'LineWidth', 2); % Ellipse arc
% yline(y1_mean, 'b--', 'LineWidth', 2);       % First horizontal line
% yline(y2_mean, 'b-.', 'LineWidth', 2);       % Second horizontal line
% plot(tp_ref, y2_mean, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
% xlabel('Depth (mm)');
% ylabel('Width (mm)');
% grid on;
% xlim([0 20]);
% ylim([0 3]);
% legend("width-depth data", "work zone", "hypodermic tube", "filament", "work zone start point")
% set(gca, "FontSize", 20);
% hold off;


% figure;
% subplot(1,3,1);
% imshow(image);
% title("Original");
% set(gca, "FontSize", 20);
% 
% subplot(1,3,2);
% imshow(background_image);
% title("Background");
% set(gca, "FontSize", 20);
% 
% subplot(1,3,3);
% imshow(object_image);
% title("Object of Interest (OI)");
% set(gca, "FontSize", 20);
% 
% figure;
% subplot(1,3,1);
% imshow(bw);
% title("Binarize & Denoise");
% set(gca, "FontSize", 20);
% 
% subplot(1,3,2);
% imshow(image);
% hold on;
% plot(x_all, y_all, 'r.', "MarkerSize", 5);
% hold off;
% title("OI Outlined");
% set(gca, "FontSize", 20);
% 
% subplot(1,3,3);
% imshow(image);
% hold on;
% plot(x_blob, y_blob, 'r.', "MarkerSize", 5);
% plot(x_range, y_range, 'r.', 'MarkerSize', 5);
% hold off;
% title("Work Zone Outlined");
% set(gca, "FontSize", 20);

% figure;
% subplot(1,6,1);
% imshow(image);
% title("Original");
% set(gca, "FontSize", 16);
% 
% subplot(1,6,2);
% imshow(background_image);
% title("Background");
% set(gca, "FontSize", 16);
% 
% subplot(1,6,3);
% imshow(object_image);
% title("Object of Interest (OI)");
% set(gca, "FontSize", 16);
% 
% subplot(1,6,4);
% imshow(bw);
% title("Binarize & Denoise");
% set(gca, "FontSize", 16);
% 
% subplot(1,6,5);
% imshow(image);
% hold on;
% plot(x_all, y_all, 'r.', "MarkerSize", 5);
% hold off;
% title("OI Outlined");
% set(gca, "FontSize", 16);
% 
% subplot(1,6,6);
% imshow(image);
% hold on;
% plot(x_blob, y_blob, 'r.', "MarkerSize", 5);
% plot(x_range, y_range, 'r.', 'MarkerSize', 5);
% hold off;
% title("Work Zone Outlined");
% set(gca, "FontSize", 16);