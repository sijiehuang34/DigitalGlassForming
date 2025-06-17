%% GlassBlob_vs_LaserData_Exp2_workzone.m
% This program is used to process multiple frames.
clc; clear; close all;

%% Specify the dataset location here
dataDir = 'E:/SFF/Exp2_MeasuringTracks_DiffStartsStops/3_ConstantLaserPower_50W_0o5mms_fp_7mm/';
if ~isfolder(dataDir)
    error('Directory "%s" does not exist. Please check the path and try again.', dataDir);
end

%% Read the csv file and create time vector here
disp('Now reading CPPlog.csv file...');
cpp_file = [dataDir  '/CPPlog.csv'];
data_nonTrim = readtable(cpp_file);
data_nonTrim = data_nonTrim(1:end-1,:); % Remove last entry which is populated with 0s

% Create time vector and add the column to the dataset
startTimeSec = data_nonTrim.StartTime_sec_;
startTimeNanosec = data_nonTrim.StartTime_nanosec_;
endTimeProcessSec = data_nonTrim.EndTimeProcess_sec_;
endTimeProcessNanosec = data_nonTrim.EndTimeProcess_nanosec_;
endTimeLoopSec = data_nonTrim.EndTimeLoop_sec_;
endTimeLoopNanosec = data_nonTrim.EndTimeLoop_nanosec_;

totalStartTime = mergeTimeData(startTimeSec, startTimeNanosec);
totalEndTimeProcess = mergeTimeData(endTimeProcessSec, endTimeProcessNanosec);
totalEndTimeLoop = mergeTimeData(endTimeLoopSec, endTimeLoopNanosec);

data_nonTrim.process_time = totalEndTimeProcess - totalStartTime;
data_nonTrim.loop_time = totalEndTimeLoop - totalStartTime;
data_nonTrim.IterationNum_unit_ =  data_nonTrim.IterationNum_unit_ + 1;
data_nonTrim.iter_counter = data_nonTrim.IterationNum_unit_(1:end);

timeVector = totalStartTime - totalStartTime(1);
data_nonTrim.timeVector = timeVector;

%% Read and associate image files
disp('Now associating image files...');

imagesFolder = fullfile(dataDir, 'images');
if ~isfolder(imagesFolder)
    error('Folder "images" does not exist in the specified directory.');
end

imageFiles = dir(fullfile(imagesFolder, 'Acquisition*.jpeg'));
if isempty(imageFiles)
    error('No JPEG images found inside the "images" folder.');
end

imageNames = {imageFiles.name};
imageNumbers = cellfun(@(x) sscanf(x, 'Acquisition%d.jpeg'), imageNames);
[~, sortIdx] = sort(imageNumbers);
imageFiles = imageFiles(sortIdx);

imageMap = containers.Map(imageNumbers, fullfile(imagesFolder, imageNames));
data_nonTrim.image_path = repmat({NaN}, height(data_nonTrim), 1);

for i = 1:height(data_nonTrim)
    iterNum = data_nonTrim.IterationNum_unit_(i);
    if isKey(imageMap, iterNum)
        data_nonTrim.image_path{i} = imageMap(iterNum);
    end
end

%% Visual processing parameters
vidSpeed = 10;
sampleInterval = 10;
errorThreshold = 1.15;
one_px = 12.446e-6; % meters
w_hypodermic = 2e-3 / one_px; %px
w_filament = 1.35e-3 / one_px; %px
minSize = 2000;
se = strel('disk', 3);
substrate_start_row = 1350;
bg_oi_threshold = 10;

firstValidAcq = 390;
% lastValidAcq = find(~cellfun(@(x) isnumeric(x) && isnan(x), data_nonTrim.image_path), 1, 'last');
lastValidAcq = 1070;
fprintf('The largest row index without NaN in image_path is: %d\n', lastValidAcq);

data_nonTrim.wz_area = NaN(height(data_nonTrim), 1);
data_nonTrim.wz_centroid = NaN(height(data_nonTrim), 2);

disp('Now looping through frames...');
idx = firstValidAcq;

while idx <= lastValidAcq
    rowIdx = find(data_nonTrim.IterationNum_unit_ == idx);
    if isempty(rowIdx) || (isnumeric(data_nonTrim.image_path{rowIdx}) && isnan(data_nonTrim.image_path{rowIdx}))
        idx = idx + sampleInterval;
        continue;
    end

    image = imread(data_nonTrim.image_path{rowIdx});
    [image_height, image_width, ~] = size(image);

    lab_image = rgb2lab(image);
    a_channel = lab_image(:,:,2);
    b_channel = lab_image(:,:,3);

    region_width = 0.5 * (image_width - 200);
    region1_a = a_channel(1200:substrate_start_row-1, 1:region_width);
    region1_b = b_channel(1200:substrate_start_row-1, 1:region_width);
    region1_a2 = a_channel(800:substrate_start_row-101, image_width-region_width:end);
    region1_b2 = b_channel(800:substrate_start_row-101, image_width-region_width:end);
    region2_a = a_channel(1:500, 1:region_width);
    region2_b = b_channel(1:500, 1:region_width);
    region2_a2 = a_channel(1:500, image_width-region_width:end);
    region2_b2 = b_channel(1:500, image_width-region_width:end);

    bg1_a = mean([mean2(region1_a), mean2(region1_a2)]); 
    bg1_b = mean([mean2(region1_b), mean2(region1_b2)]);
    bg2_a = mean([mean2(region2_a), mean2(region2_a2)]); 
    bg2_b = mean([mean2(region2_b), mean2(region2_b2)]);

    distance_bg1 = sqrt((a_channel - bg1_a).^2 + (b_channel - bg1_b).^2);
    distance_bg2 = sqrt((a_channel - bg2_a).^2 + (b_channel - bg2_b).^2);

    threshold = bg_oi_threshold;
    bg_mask_1 = distance_bg1 < threshold;
    bg_mask_2 = distance_bg2 < threshold;
    substrate_mask = false(image_height, image_width);
    substrate_mask(substrate_start_row:end, :) = true;
    bg_mask_1(substrate_mask) = false;
    bg_mask_2(substrate_mask) = false;

    object_mask = ~bg_mask_1 & ~bg_mask_2 & ~substrate_mask;
    bw = uint8(object_mask) * 255;
    bw = bwareaopen(bw, minSize);
    bw = imfill(bw, 'holes');
    bw_smooth = imopen(bw, se);
    bw_smooth = imclose(bw_smooth, se);
    bw_smooth = bwareaopen(bw_smooth, minSize);

    for row = 1:min(100, size(bw_smooth, 1))
        row_data = bw_smooth(row,:);
        w = find(row_data);
        if ~isempty(w)
            bw_smooth(row, w(1):w(end)) = 1;
        end
    end

    object_image = image;
    object_image(repmat(~logical(bw_smooth), [1 1 3])) = 0;
    [r, c, ~] = find(all(object_image >= 240, 3));
    workzonePoints = [c, r];
    wz_centroid = [NaN, NaN];
    wz_area = NaN;
    if ~isempty(workzonePoints)
        clusterLabels = dbscan(double(workzonePoints), 5, 50);
        validLabels = unique(clusterLabels(clusterLabels > 0));
        if ~isempty(validLabels)
            clusterSizes = arrayfun(@(x) sum(clusterLabels == x), validLabels);
            [~, maxIdx] = max(clusterSizes);
            keepLabel = validLabels(maxIdx);
            validPoints = workzonePoints(clusterLabels == keepLabel, :);
            wz_area = size(validPoints,1) * one_px^2 * 1e6; % mm^2
   
            % --- Compute hypodermic_left_edge from top 10 rows of bw_smooth ---
            whitePixelXs = [];
            
            for row = 1:min(10, size(bw_smooth, 1))
                for col = 1:size(bw_smooth, 2)
                    if bw_smooth(row, col) > 0
                        whitePixelXs(end+1) = col; % store X coordinate
                        break; % stop at first white pixel in this row
                    end
                end
            end
            
            % Take the minimum X found
            if ~isempty(whitePixelXs)
                hypodermic_left_edge = min(whitePixelXs); % px
            else
                hypodermic_left_edge = NaN; % fallback if no white pixel is found
            end
            
            % Compute origin for X and shift both axes before converting to mm
            centroidX_rel_origin = hypodermic_left_edge + (0.5 * w_hypodermic);  % px
            
            wz_centroid_px = mean(validPoints, 1);
            wz_centroid_px(1) = wz_centroid_px(1) - centroidX_rel_origin;         % X shift
            wz_centroid_px(2) = substrate_start_row - wz_centroid_px(2);         % Y shift
            
            wz_centroid = wz_centroid_px * one_px * 1e3; % convert to mm

        end
    end

    data_nonTrim.wz_area(rowIdx) = wz_area;
    data_nonTrim.wz_centroid(rowIdx, :) = wz_centroid;

    fprintf('Frame %d processed.\n', idx);
    idx = idx + sampleInterval;
end

%% Save and plot
save(fullfile(dataDir, 'processed_blob_wz_summary.mat'), 'data_nonTrim');

figure('Units', 'inches', 'Position', [1, 1, 12, 4.5], 'Color', 'w');

% --- Plot 1: Workzone Area vs Time ---
subplot(1,2,1);
plot(data_nonTrim.timeVector, data_nonTrim.wz_area, 'k.', 'MarkerSize', 10);
xlabel('Time (s)');
ylabel('Area (mm²)');
grid on;
ylim([0 10]);
set(gca, 'FontSize', 16);

% --- Plot 2: Centroid X and Y vs Time with dual y-axes ---
subplot(1,2,2);
ax = gca; % get current axis

yyaxis left;
ax.YAxis(1).Color = [0 0.4470 0.7410]; % MATLAB blue
plot(data_nonTrim.timeVector, data_nonTrim.wz_centroid(:,1), 'o', ...
     'Color', ax.YAxis(1).Color, 'LineWidth', 1);
ylabel('X (mm)');
ylim([0 3]);

yyaxis right;
ax.YAxis(2).Color = [0.8500 0.3250 0.0980]; % MATLAB red
plot(data_nonTrim.timeVector, data_nonTrim.wz_centroid(:,2), 's', ...
     'Color', ax.YAxis(2).Color, 'LineWidth', 1);
ylabel('Y (mm)');
ylim([0 3]);
xlabel('Time (s)');
grid on;
set(gca, 'FontSize', 16);
legend('Centroid X', 'Centroid Y', 'FontSize', 14, 'Location', 'best');

% Save figure
exportgraphics(gcf, fullfile(dataDir, 'Figure_10_AreaCentroidTime_DualAxis.png'), 'Resolution', 300);
savefig(gcf, fullfile(dataDir, 'Figure_10_AreaCentroidTime_DualAxis.fig'));

disp('Program done executing.');

%% Helper function
function totalSeconds = mergeTimeData(secondsVec, nanosecondsVec)
    totalSeconds = double(secondsVec) + double(nanosecondsVec) * 1e-9;
end
