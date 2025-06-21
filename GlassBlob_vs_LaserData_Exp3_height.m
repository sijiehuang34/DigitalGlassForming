%% GlassBlob_vs_LaserData_Exp3_height_Spatiotemporal.m
% This script processes a series of image frames to extract spatiotemporal
% height profile data and visualize it as a 3D waterfall plot.
clc; clear; close all;

%% Specify dataset location
baseDir = 'E:/SFF/Exp3_JoiningTracks/18L_0.5mmsPRINT_40W_20mm_4mmAboveFocus_BlobTip_0mmOverlap_1mmsRETRACT/';
dataDir = baseDir;
imagesFolder = fullfile(dataDir, 'images');

if ~isfolder(dataDir) || ~isfolder(imagesFolder)
    error('Check that data and image folders exist.');
end

%% Load CPP log
cpp_file = fullfile(dataDir, 'CPPlog.csv');
data = readtable(cpp_file);
data = data(1:end-1,:);

startTime = mergeTimeData(data.StartTime_sec_, data.StartTime_nanosec_);
endTime = mergeTimeData(data.EndTimeLoop_sec_, data.EndTimeLoop_nanosec_);
data.timeVector = startTime - startTime(1);
data.IterationNum_unit_ = data.IterationNum_unit_ + 1;

%% Map image paths
imageFiles = dir(fullfile(imagesFolder, 'Acquisition*.jpeg'));
imageNames = {imageFiles.name};
imageNumbers = cellfun(@(x) sscanf(x, 'Acquisition%d.jpeg'), imageNames);
[~, sortIdx] = sort(imageNumbers);
imageFiles = imageFiles(sortIdx);
imageMap = containers.Map(imageNumbers, fullfile(imagesFolder, imageNames));

data.image_path = repmat({NaN}, height(data), 1);
for i = 1:height(data)
    iterNum = data.IterationNum_unit_(i);
    if isKey(imageMap, iterNum)
        data.image_path{i} = imageMap(iterNum);
    end
end

%% Constants and parameters
one_px = 12.446e-6;
substrate_start_row = 1235;
sampleInterval = 10;
se = strel('disk', 4);
threshold = 12;
minSize = 2000;
height_x_exclude_start = 400;
height_x_exclude_end = 650;
epsilon = 0.30;    % Distance threshold in mm (adjust as needed)

firstValid = 480;
lastValid = 610;

%% Initialize output
Z_profiles = {}; X_profiles = {}; T_values = [];

%% Frame processing loop
for idx = firstValid:sampleInterval:lastValid
    rowIdx = find(data.IterationNum_unit_ == idx);
    if isempty(rowIdx) || isnumeric(data.image_path{rowIdx})
        continue;
    end

    image = imread(data.image_path{rowIdx});
    [h, w, ~] = size(image);
    lab = rgb2lab(image);
    a = lab(:,:,2); b = lab(:,:,3);

    regW = 0.5 * (w - 200);
    bg1_a = mean([mean2(a(1100:substrate_start_row-1, 1:regW)), mean2(a(800:substrate_start_row-101, end-regW:end))]);
    bg1_b = mean([mean2(b(1100:substrate_start_row-1, 1:regW)), mean2(b(800:substrate_start_row-101, end-regW:end))]);
    bg2_a = mean([mean2(a(1:500, 1:regW)), mean2(a(1:500, end-regW:end))]);
    bg2_b = mean([mean2(b(1:500, 1:regW)), mean2(b(1:500, end-regW:end))]);

    d1 = sqrt((a - bg1_a).^2 + (b - bg1_b).^2);
    d2 = sqrt((a - bg2_a).^2 + (b - bg2_b).^2);
    mask1 = d1 < threshold; mask2 = d2 < threshold;
    submask = false(h, w); submask(substrate_start_row:end,:) = true;
    mask1(submask) = false; mask2(submask) = false;
    objMask = ~mask1 & ~mask2 & ~submask;

    bw = uint8(objMask) * 255;
    bw = bwareaopen(bw, minSize);
    bw = imfill(bw, 'holes');
    bw = imopen(bw, se);
    bw = imclose(bw, se);
    bw = bwareaopen(bw, minSize);

    for row = 1:min(100, size(bw,1))
        whiteIdx = find(bw(row,:));
        if ~isempty(whiteIdx)
            bw(row, whiteIdx(1):whiteIdx(end)) = 1;
        end
    end

    boundary = [];
    for col = 1:w
        row_idx = find(bw(:,col), 1, 'first');
        if ~isempty(row_idx)
            boundary = [boundary; col, row_idx];
        end
    end

    [sortedX, sortIdx] = sort(boundary(:,1), 'descend');
    sortedY = boundary(sortIdx, 2);
    dy = diff(sortedY);
    jumps = find(abs(dy) > 200);

    if length(jumps) >= 2
        x1 = sortedX(jumps(1)); x2 = sortedX(jumps(2));
        x_ex1 = min(x1, x2); x_ex2 = max(x1, x2);
    else
        x_ex1 = -Inf; x_ex2 = -Inf;
    end
    mask = (boundary(:,1) > x_ex1 & boundary(:,1) < x_ex2) | ...
           (boundary(:,1) > height_x_exclude_start & boundary(:,1) < height_x_exclude_end);
    boundary = boundary(~mask,:);

    % Convert to physical units
    Z = (substrate_start_row - boundary(:,2)) * one_px * 1000;  % mm
    X = boundary(:,1) * one_px * 1000;                          % mm
    
    % Apply DBSCAN to remove noise (assumes 2 dominant clusters)
    dataPoints = [X, Z];
    minpts = 5;        % Minimum number of points per cluster
    labels = dbscan(dataPoints, epsilon, minpts);
    
    % Keep only the two largest clusters
    validClusters = labels > 0;
    if any(validClusters)
        uniqueClusters = unique(labels(validClusters));
        clusterSizes = arrayfun(@(c) sum(labels == c), uniqueClusters);
    
        [~, sortedIdx] = sort(clusterSizes, 'descend');
        topTwoLabels = uniqueClusters(sortedIdx(1:min(2, numel(sortedIdx))));
        clusterMask = ismember(labels, topTwoLabels);
    
        % Filter coordinates
        X = X(clusterMask);
        Z = Z(clusterMask);
    else
        % If DBSCAN finds no valid clusters, skip frame
        continue;
    end
    
    % Store results
    Z_profiles{end+1} = Z(:)';
    X_profiles{end+1} = X(:)';
    T_values(end+1) = data.timeVector(rowIdx);

end

%% Pad for plotting
maxLen = max(cellfun(@length, Z_profiles));
Z_mat = NaN(length(Z_profiles), maxLen);
X_mat = NaN(length(X_profiles), maxLen);
for i = 1:length(Z_profiles)
    Z_mat(i,1:length(Z_profiles{i})) = Z_profiles{i};
    X_mat(i,1:length(X_profiles{i})) = X_profiles{i};
end

%% Waterfall plot (lines for each cluster, no connecting across gaps)
figure('Color','w','Units','inches','Position',[1 1 11 4]); hold on;
colors = lines(length(Z_profiles));  % One color per frame

for i = 1:length(Z_profiles)
    x_row = X_profiles{i};
    z_row = Z_profiles{i};
    t_row = T_values(i) * ones(size(x_row));

    % Apply DBSCAN again per row to break into clusters (optional but safe)
    if length(x_row) < 2
        continue;
    end

    dataPoints = [x_row(:), z_row(:)];
    labels = dbscan(dataPoints, 0.15, 5);
    clusterIDs = unique(labels(labels > 0));

    for c = 1:numel(clusterIDs)
        cid = clusterIDs(c);
        mask = (labels == cid);
        plot3(x_row(mask), t_row(mask), z_row(mask), '-', ...
              'Color', colors(i,:), 'LineWidth', 1.5);
    end
end

xlabel('X (mm)', 'FontSize', 16);
ylabel('Time (s)', 'FontSize', 16);
zlabel('Z (mm)', 'FontSize', 16);
grid on;
view(20, 50);
set(gca, 'FontSize', 16);

exportgraphics(gcf, fullfile(dataDir, 'Figure_11_HeightProfile_Spatiotemporal_2.png'), 'Resolution', 300);
savefig(gcf, fullfile(dataDir, 'Figure_11_HeightProfile_Spatiotemporal_2.fig'));

%% Helper
function totalSeconds = mergeTimeData(sec, nsec)
    totalSeconds = double(sec) + double(nsec) * 1e-9;
end
