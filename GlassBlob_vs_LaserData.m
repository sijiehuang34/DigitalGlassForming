%% GlassBlob_vs_LaserData.m
clc; clear; close all;

%% Specify the dataset location here
dataDir = 'E:/SFF/Exp1_ConstantPowerExp/12_20W_Filament2mmBelowFocus_Substrate6mmBelowFocus_BlobTip/'; % New dataset location
if ~isfolder(dataDir)
    error('Directory "%s" does not exist. Please check the path and try again.', dataDir);
end

%% Read the csv file and create time vector here
disp('Now reading CPPlog.csv file...');
cpp_file = [dataDir  '/CPPlog.csv'];
data_nonTrim = readtable(cpp_file);
clear cpp_file

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

clear startTimeSec startTimeNanosec endTimeProcessSec endTimeProcessNanosec endTimeLoopSec endTimeLoopNanosec totalStartTime totalEndTimeProcess totalEndTimeLoop timeVector

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

data_nonTrim.image_path = repmat({NaN}, height(data_nonTrim), 1);

% Create a map from IterationNum_unit_ to the corresponding image path
imageMap = containers.Map(imageNumbers, fullfile(imagesFolder, imageNames));

% Loop through the rows of the data_nonTrim table
for i = 1:height(data_nonTrim)
    iterNum = data_nonTrim.IterationNum_unit_(i);
    if isKey(imageMap, iterNum)
        data_nonTrim.image_path{i} = imageMap(iterNum);  % Assign the image path
    end
end

% Find the last valid image acquisition
lastValidAcq = NaN;  % Initialize
for i = 1:height(data_nonTrim)
    if ~isnan(data_nonTrim.image_path{i})
        lastValidAcq = i;
    end
end
fprintf('The largest row index without NaN in image_path is: %d\n', lastValidAcq);


clear imageNumbers imageNames imageMap i iterNum

%% Visual processing here
% Video specs
vidSpeed = 10;

% User-defined constants
sampleInterval = 10;
skipFlag = false;
errorThreshold = 1.15; % remove area that's 1.25x bigger than previous one

one_px = 12.446e-6; % meters
w_hypodermic = 2e-3 / one_px;
w_filament = 1.4e-3 / one_px;

minSize = 1000;
se = strel('disk', 3);
substrate_start_row = 2000;
bg_oi_threshold = 10;

known_hypodermic_depth = 190;
known_filament_start_depth = known_hypodermic_depth + 100; %px

firstValidAcq = 1;

% Preallocate new columns
data_nonTrim.depth = NaN(height(data_nonTrim), 1);
data_nonTrim.width = NaN(height(data_nonTrim), 1);
data_nonTrim.blob_area = NaN(height(data_nonTrim), 1);
data_nonTrim.blob_height = NaN(height(data_nonTrim), 1);
data_nonTrim.x_blob = cell(height(data_nonTrim), 1);
data_nonTrim.y_blob = cell(height(data_nonTrim), 1);
data_nonTrim.x_range = cell(height(data_nonTrim), 1);
data_nonTrim.y_range = cell(height(data_nonTrim), 1);
data_nonTrim.x_all = cell(height(data_nonTrim), 1);
data_nonTrim.y_all = cell(height(data_nonTrim), 1);

disp('Now looping through frames...');

% Initialize loop
prevArea = NaN;
idx = firstValidAcq;

while idx <= lastValidAcq

    rowIdx = find(data_nonTrim.IterationNum_unit_ == idx);
    if isempty(rowIdx) || (isnumeric(data_nonTrim.image_path{rowIdx}) && isnan(data_nonTrim.image_path{rowIdx}))
        idx = idx + sampleInterval;
        continue;
    end


    image = imread(data_nonTrim.image_path{rowIdx});
    [image_height, image_width, ~] = size(image);

%     substrate_start_row = floor(image_height * (5/6)) + 1 + finetune_substrate;
    substrate_mask = false(image_height, image_width);
    substrate_mask(substrate_start_row:end, :) = true;

    lab_image = rgb2lab(image);
    a_channel = lab_image(:,:,2);
    b_channel = lab_image(:,:,3);
    region_height = substrate_start_row - 1;
    region_width = 0.5 * (image_width - 200); % 1/2 of (width - thickness of tube)
    region1_a = a_channel(1:region_height, 1:region_width);
    region1_b = b_channel(1:region_height, 1:region_width);
    region2_a = a_channel(1:region_height, end-region_width:end);
    region2_b = b_channel(1:region_height, end-region_width:end);
    background_a = mean([mean2(region1_a), mean2(region2_a)]);
    background_b = mean([mean2(region1_b), mean2(region2_b)]);
    distance = sqrt((a_channel - background_a).^2 + (b_channel - background_b).^2);
    background_mask = distance < bg_oi_threshold;
    background_mask(substrate_mask) = false;
    object_mask = ~background_mask & ~substrate_mask;
    bw = uint8(object_mask) * 255;
    bw = bwareaopen(bw, minSize);
    bw = imfill(bw, 'holes');
    bw_smooth = imopen(bw, se);
    bw_smooth = imclose(bw_smooth, se);
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

    [Boundary, ~, ~] = bwboundaries(bw_smooth, 4, 'noholes');
    [~, sortIdxB] = sort(cellfun(@(x) size(x, 1), Boundary), 'descend');
    Boundary = Boundary(sortIdxB);
    Boundary = Boundary(cellfun(@(x) size(x,1) >= 2000, Boundary));
    Boundary = Boundary(~cellfun(@(x) any(x(:,2) == 1 | x(:,2) == image_width), Boundary));

    if isempty(Boundary)
        idx = idx + sampleInterval;
        continue;
    end

    boundaryPoints = Boundary{1};
    [sortedY, sortIdx] = sort(boundaryPoints(:,1), 'ascend');
    sortedX = boundaryPoints(sortIdx, 2);
    boundaryPoints = [sortedY, sortedX];
    x_all = boundaryPoints(:,2);
    y_all = boundaryPoints(:,1);

    uniqueY = unique(boundaryPoints(:,1), 'sorted');
    n = length(uniqueY);
    wid_vs_dep = zeros(n, 2);
    wid_vs_dep(:,1) = uniqueY;
    for j = 1:n
        yVal = uniqueY(j);
        xVals = boundaryPoints(boundaryPoints(:,1) == yVal, 2);
        wid_vs_dep(j,2) = abs(max(xVals) - min(xVals));
    end

    depth = wid_vs_dep(:,1);
    width = wid_vs_dep(:,2);
    idx1 = (depth >= 0) & (depth <= known_hypodermic_depth);
    y1_mean = mean(width(idx1));
    y2_mean = y1_mean - (w_hypodermic - w_filament);
    valid_idx = find(depth > known_filament_start_depth & width > y1_mean, 1, 'first');
    if isempty(valid_idx), idx = idx + sampleInterval; continue; end
    rough_blob_start = depth(valid_idx);
    idx_curve = (depth >= rough_blob_start);
    x_curve = depth(idx_curve);
    y_curve = width(idx_curve);
    p = polyfit(x_curve, y_curve, 2);
%     p = polyfit(x_curve, y_curve, 4);
    tp = roots([p(1), p(2), p(3) - y2_mean]);
%     tp = roots([p(1), p(2), p(3), p(4), p(5) - y2_mean]);
    tp = tp(imag(tp) == 0);
    if isempty(tp), idx = idx + sampleInterval; continue; end
    tp_ref = min(tp);

    x_roots = roots(p);
    x_roots = x_roots(imag(x_roots) == 0);
    if length(x_roots) < 2, idx = idx + sampleInterval; continue; end
    x1 = min(x_roots);
    x2 = max(x_roots);
    x_blob = linspace(x1, x2, 1000);
    y_blob = polyval(p, x_blob);
    y_blob(y_blob < 0) = 0;
    area_px2 = trapz(x_blob, y_blob);
    area_mm2 = area_px2 * one_px^2 * 1e6;
    blob_substrate_dist = (substrate_start_row - max(depth)) * one_px * 1000;

%     % Error flag for unreliable area detection
%     if ~isnan(prevArea) && area_mm2 > errorThreshold * prevArea
%         skipFlag = true;
%         fprintf('Data from frame %d is unreliable and abandoned.\n', idx);
%         idx = idx + sampleInterval;
%         continue;
%     end

    prevArea = area_mm2;

    blob_outline = boundaryPoints(boundaryPoints(:,1) > tp_ref, :);
    x_blob_pts = blob_outline(:,2);
    y_blob_pts = blob_outline(:,1);
    y_opening = min(y_blob_pts);
    x_opening = x_blob_pts(y_blob_pts == y_opening);
    x_range = min(x_opening):max(x_opening);
    y_range = y_opening * ones(size(x_range));

    data_nonTrim.depth(rowIdx) = max(depth) * one_px * 1000;
    data_nonTrim.width(rowIdx) = max(width) * one_px * 1000;
    data_nonTrim.blob_area(rowIdx) = area_mm2;
    data_nonTrim.blob_height(rowIdx) = blob_substrate_dist;
    data_nonTrim.x_blob{rowIdx} = x_blob_pts;
    data_nonTrim.y_blob{rowIdx} = y_blob_pts;
    data_nonTrim.x_range{rowIdx} = x_range;
    data_nonTrim.y_range{rowIdx} = y_range;
    data_nonTrim.x_all{rowIdx} = x_all;
    data_nonTrim.y_all{rowIdx} = y_all;

    fprintf('Frame %d processed.\n', idx);
    idx = idx + sampleInterval;

    % Reset error flag
    skipFlag = false;
end

clear a_channel b_channel background_a background_b area_mm2 area_px2 ...
bg_oi_threshold blob_substrate_dist depth idx idx1 ...
idx_curve image_height image_width j lab_image minSize n region1_a region1_b...
region2_a region2_b region_height region_width rough_blob_start rowIdx se...
sortedX sortedY sortIdx sortIdxB substrate_start_row valid_idx width ...
wid_vs_dep x1 x2 x_blob_pts y_blob_pts x_range y_range x_all y_all xVals yVal

%% Plot only
figure;
subplot(1,2,1);
plot(data_nonTrim.timeVector, data_nonTrim.blob_area, 'b.-', 'LineWidth', 2, 'MarkerSize', 15);
xlabel('Time (s)');
ylabel('Area (mm²)');
grid on;
% set(gca, "FontSize", 30);
ax = gca;
% ax.FontSize = 40;       % Axis tick label size
ax.LineWidth = 1.5;       % Thicken the outer box
box on;                 % Make sure the box is visible
% xlim([250 400]);
subplot(1,2,2);
plot(data_nonTrim.timeVector, data_nonTrim.blob_height, 'r.-', 'LineWidth', 2, 'MarkerSize', 15);
xlabel('Time (s)');
ylabel('Height (mm)');
grid on;
% set(gca, "FontSize", 30);
ax = gca;
% ax.FontSize = 40;       % Axis tick label size
ax.LineWidth = 1.5;       % Thicken the outer box
box on;                 % Make sure the box is visible
hold off;
% xlim([250 400]);

% Save the figure to the specified folder
exportgraphics(gcf, fullfile(dataDir, 'Figure_10_AreaHeightTime.png'), 'Resolution', 300);
savefig(gcf, fullfile(dataDir, 'Figure_10_AreaHeightTime.fig'));

% %%%%%%%%%%%%%%%%%
% 
% 
% % plot(depth, width, 'k.', 'MarkerSize', 10); 
% % hold on;
% % plot(x_extended, y_fit, 'g-', 'LineWidth', 2); % Ellipse arc
% % yline(y1_mean, 'b:', 'LineWidth', 2);         % First horizontal line
% % yline(y2_mean, 'b--', 'LineWidth', 2);        % Second horizontal line
% % plot(tp_ref, y2_mean, 'ro', 'MarkerSize', 15, 'LineWidth', 2);
% % 
% % xlabel('Depth (mm)', 'FontSize', 40);
% % ylabel('Width (mm)', 'FontSize', 40);
% % 
% % grid on;
% % xlim([0 20]);
% % ylim([0 3]);
% % 
% % legend("width-depth data", "work zone", "hypodermic tube", "cold filament", "work zone start point")
% % 
% % ax = gca;
% % ax.FontSize = 20;       % Axis tick label size
% % ax.LineWidth = 1.5;       % Thicken the outer box
% % box on;                 % Make sure the box is visible
% % 
% % hold off;
% % 
% % % Save with DPI control
% % dataDir = 'D:/ND Lab/'; % New dataset location
% % print(fig, fullfile(dataDir, 'Plotssss'), '-djpeg', '-r300'); % 300 DPI JPEG
% % savefig(fig, fullfile(dataDir, 'Plot.fig'));
% 
% %%%%%%%%%%%%%%%
% 
% % xlim([200 max(timeVecAll)]);
% % ylim([0 max(areaVecAll) * 1.1]);
% 
% %% Create a video of width vs. depth + blob outline
% % disp('Now creating Video 1...');
% % 
% % outputVideo = VideoWriter(fullfile(dataDir, 'BlobWidthVsDepth.avi'));
% % outputVideo.FrameRate = vidSpeed;
% % open(outputVideo);
% % 
% % for idx = 1760:sampleInterval:3756
% %     rowIdx = find(data_nonTrim.IterationNum_unit_ == idx);
% %     if isempty(rowIdx) || ...
% %        isnan(data_nonTrim.depth(rowIdx)) || ...
% %        isempty(data_nonTrim.x_blob{rowIdx}) || ...
% %        isempty(data_nonTrim.x_all{rowIdx})
% %         continue;
% %     end
% % 
% %     % Load data
% %     image = imread(data_nonTrim.image_path{rowIdx});
% %     x_all = data_nonTrim.x_all{rowIdx};
% %     y_all = data_nonTrim.y_all{rowIdx};
% %     x_blob = data_nonTrim.x_blob{rowIdx};
% %     y_blob = data_nonTrim.y_blob{rowIdx};
% %     x_range = data_nonTrim.x_range{rowIdx};
% %     y_range = data_nonTrim.y_range{rowIdx};
% % 
% %     % Recompute depth and width
% %     uniqueY = unique(y_all, 'sorted');
% %     width_vec = zeros(size(uniqueY));
% %     for j = 1:length(uniqueY)
% %         yVal = uniqueY(j);
% %         xVals = x_all(y_all == yVal);
% %         width_vec(j) = abs(max(xVals) - min(xVals));
% %     end
% %     depth_mm = uniqueY * one_px * 1000;
% %     width_mm = width_vec * one_px * 1000;
% % 
% %     % Parabola fit
% %     y1_mean = mean(width_vec(uniqueY >= 0 & uniqueY <= 500));
% %     w_hypodermic = 1.7e-3 / one_px;
% %     w_filament = 1.4e-3 / one_px;
% %     y2_mean = y1_mean - (w_hypodermic - w_filament);
% %     valid_idx = find(uniqueY > 600 & width_vec > y1_mean, 1, 'first');
% %     if isempty(valid_idx), continue; end
% %     rough_blob_start = uniqueY(valid_idx);
% %     idx_curve = (uniqueY >= rough_blob_start);
% %     p = polyfit(uniqueY(idx_curve), width_vec(idx_curve), 2);
% %     tp = roots([p(1), p(2), p(3) - y2_mean]);
% %     tp = tp(imag(tp) == 0);
% %     if isempty(tp), continue; end
% %     tp_ref = min(tp);
% % 
% %     % Prepare figure
% %     figure('Visible', 'off', 'Position', [100 100 1200 500]);
% % 
% %     % Left plot: Width vs. depth
% %     subplot(1,2,1);
% %     x_ext = linspace(0, max(uniqueY)*1.2, 1000);
% %     y_fit = polyval(p, x_ext);
% %     plot(depth_mm, width_mm, 'k.', 'MarkerSize', 5);
% %     hold on;
% %     plot(x_ext * one_px * 1000, y_fit * one_px * 1000, 'g-', 'LineWidth', 2);
% %     yline(y1_mean * one_px * 1000, 'b-', 'LineWidth', 1.5);
% %     yline(y2_mean * one_px * 1000, 'b-', 'LineWidth', 1.5);
% %     plot(tp_ref * one_px * 1000, y2_mean * one_px * 1000, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
% %     xlabel('Depth (mm)');
% %     ylabel('Width (mm)');
% % %     title(sprintf('Width vs. Depth — Frame %d', idx));
% %     grid on;
% %     xlim([0 20]); ylim([0 3]);
% %     hold off;
% % 
% %     % Right plot: Original image + blob outline
% %     subplot(1,2,2);
% %     imshow(image);
% %     hold on;
% %     plot(x_blob, y_blob, 'r.', 'MarkerSize', 3);
% %     plot(x_range, y_range, 'r.', 'MarkerSize', 3);
% %     hold off;
% % %     title('Blob Outlined on Original Image');
% % 
% %     % Get frame
% %     frame = getframe(gcf);
% %     writeVideo(outputVideo, frame);
% %     close(gcf);
% % end
% % 
% % close(outputVideo);
% % fprintf('Video saved to: %s\n', fullfile(dataDir, 'BlobWidthVsDepth.avi'));
% 
% %% Create a video of blob area vs. time (growing plot) + blob outline
% % disp('Now creating Video 2...');
% % 
% % outputVideo2 = VideoWriter(fullfile(dataDir, 'BlobAreaVsTime.avi'));
% % outputVideo2.FrameRate = vidSpeed;
% % open(outputVideo2);
% % 
% % % Get full time/area vectors
% % timeVecAll = data_nonTrim.timeVector;
% % areaVecAll = data_nonTrim.blob_area;
% % iterVecAll = data_nonTrim.IterationNum_unit_;
% % 
% % % Get valid indices to loop through
% % validFrames = 1760:sampleInterval:3756;
% % 
% % for idx = validFrames
% %     rowIdx = find(data_nonTrim.IterationNum_unit_ == idx);
% %     if isempty(rowIdx) || ...
% %        isnan(data_nonTrim.blob_area(rowIdx)) || ...
% %        isempty(data_nonTrim.x_blob{rowIdx})
% %         continue;
% %     end
% % 
% %     image = imread(data_nonTrim.image_path{rowIdx});
% %     x_blob = data_nonTrim.x_blob{rowIdx};
% %     y_blob = data_nonTrim.y_blob{rowIdx};
% %     x_range = data_nonTrim.x_range{rowIdx};
% %     y_range = data_nonTrim.y_range{rowIdx};
% % 
% %     % Get all data up to and including this frame
% %     time_so_far = timeVecAll(iterVecAll <= idx);
% %     area_so_far = areaVecAll(iterVecAll <= idx);
% % 
% %     % Create figure
% %     figure('Visible', 'off', 'Position', [100 100 1200 500]);
% % 
% %     % Left plot: Growing blob area vs. time
% %     subplot(1,2,1);
% %     plot(time_so_far, area_so_far, 'b.-', 'LineWidth', 1.5, 'MarkerSize', 12);
% %     xlabel('Time (s)');
% %     ylabel('Area (mm²)');
% % %     title(sprintf('Blob Area vs. Time — Frame %d', idx));
% %     grid on;
% %     xlim([200 max(timeVecAll)]);
% %     ylim([0 max(areaVecAll)*1.1]);
% % 
% %     % Right plot: Original image with blob outline
% %     subplot(1,2,2);
% %     imshow(image);
% %     hold on;
% %     plot(x_blob, y_blob, 'r.', 'MarkerSize', 3);
% %     plot(x_range, y_range, 'r.', 'MarkerSize', 3);
% %     hold off;
% % %     title('Blob Outlined on Original Image');
% % 
% %     % Write to video
% %     frame = getframe(gcf);
% %     writeVideo(outputVideo2, frame);
% %     close(gcf);
% % end
% % 
% % close(outputVideo2);
% % fprintf('Video saved to: %s\n', fullfile(dataDir, 'BlobAreaVsTime.avi'));
% 
% %% Create a video of blob height vs. time (growing plot) + blob outline
% % disp('Now creating Video 3...');
% % 
% % outputVideo3 = VideoWriter(fullfile(dataDir, 'BlobHeightVsTime.avi'));
% % outputVideo3.FrameRate = vidSpeed;
% % open(outputVideo3);
% % 
% % % Get time and height vectors for all frames
% % timeVecAll = data_nonTrim.timeVector;
% % heightVecAll = data_nonTrim.blob_height;
% % iterVecAll = data_nonTrim.IterationNum_unit_;
% % 
% % validFrames = 1760:sampleInterval:3756;
% % 
% % for idx = validFrames
% %     rowIdx = find(data_nonTrim.IterationNum_unit_ == idx);
% %     if isempty(rowIdx) || ...
% %        isnan(data_nonTrim.blob_height(rowIdx)) || ...
% %        isempty(data_nonTrim.x_blob{rowIdx})
% %         continue;
% %     end
% % 
% %     image = imread(data_nonTrim.image_path{rowIdx});
% %     x_blob = data_nonTrim.x_blob{rowIdx};
% %     y_blob = data_nonTrim.y_blob{rowIdx};
% %     x_range = data_nonTrim.x_range{rowIdx};
% %     y_range = data_nonTrim.y_range{rowIdx};
% % 
% %     % Get all data up to and including this frame
% %     time_so_far = timeVecAll(iterVecAll <= idx);
% %     height_so_far = heightVecAll(iterVecAll <= idx);
% % 
% %     % Create figure
% %     figure('Visible', 'off', 'Position', [100 100 1200 500]);
% % 
% %     % Left plot: Growing blob height vs. time
% %     subplot(1,2,1);
% %     plot(time_so_far, height_so_far, 'k.-', 'LineWidth', 1.5, 'MarkerSize', 12);
% %     xlabel('Time (s)');
% %     ylabel('Height (mm)');
% % %     title(sprintf('Blob Height vs. Time — Frame %d', idx));
% %     grid on;
% %     xlim([200 max(timeVecAll)]);
% %     ylim([0 max(heightVecAll)*1.1]);
% % 
% %     % Right plot: Original image with blob outline
% %     subplot(1,2,2);
% %     imshow(image);
% %     hold on;
% %     plot(x_blob, y_blob, 'r.', 'MarkerSize', 3);
% %     plot(x_range, y_range, 'r.', 'MarkerSize', 3);
% %     hold off;
% % %     title('Blob Outlined on Original Image');
% % 
% %     % Write to video
% %     frame = getframe(gcf);
% %     writeVideo(outputVideo3, frame);
% %     close(gcf);
% % end
% % 
% % close(outputVideo3);
% % fprintf('Video saved to: %s\n', fullfile(dataDir, 'BlobHeightVsTime.avi'));
% 
% % %% Create a video of blob area + height vs. time (yyaxis) + blob outline
% % disp('Now creating Video 4...');
% % 
% % outputVideo4 = VideoWriter(fullfile(dataDir, 'BlobAreaAndHeightVsTime.avi'));
% % outputVideo4.FrameRate = vidSpeed;
% % open(outputVideo4);
% % 
% % % Get time, area, and height vectors for all frames
% % timeVecAll = data_nonTrim.timeVector;
% % areaVecAll = data_nonTrim.blob_area;
% % heightVecAll = data_nonTrim.blob_height;
% % iterVecAll = data_nonTrim.IterationNum_unit_;
% % 
% % validFrames = 1760:sampleInterval:3756;
% % 
% % for idx = validFrames
% %     rowIdx = find(data_nonTrim.IterationNum_unit_ == idx);
% %     if isempty(rowIdx) || ...
% %        isnan(data_nonTrim.blob_area(rowIdx)) || ...
% %        isnan(data_nonTrim.blob_height(rowIdx)) || ...
% %        isempty(data_nonTrim.x_blob{rowIdx})
% %         continue;
% %     end
% % 
% %     image = imread(data_nonTrim.image_path{rowIdx});
% %     x_blob = data_nonTrim.x_blob{rowIdx};
% %     y_blob = data_nonTrim.y_blob{rowIdx};
% %     x_range = data_nonTrim.x_range{rowIdx};
% %     y_range = data_nonTrim.y_range{rowIdx};
% % 
% %     % Get all data up to and including this frame
% %     time_so_far = timeVecAll(iterVecAll <= idx);
% %     area_so_far = areaVecAll(iterVecAll <= idx);
% %     height_so_far = heightVecAll(iterVecAll <= idx);
% % 
% %     % Create figure
% %     figure('Visible', 'off', 'Position', [100 100 1200 500]);
% % 
% %     % Left plot: yyaxis Area and Height vs. Time
% %     subplot(1,2,1);
% %     yyaxis left
% %     plot(time_so_far, area_so_far, 'b.-', 'LineWidth', 1.5, 'MarkerSize', 10);
% %     ylabel('Area (mm²)');
% %     ylim([0, max(areaVecAll)*1.1]);
% % 
% %     yyaxis right
% %     plot(time_so_far, height_so_far, 'r*-');
% %     ylabel('Height (mm)');
% %     ylim([0, max(heightVecAll)*1.1]);
% % 
% %     xlabel('Time (s)');
% %     grid on;
% %     xlim([200, max(timeVecAll)]);
% % 
% %     % Right plot: Original image with blob outline
% %     subplot(1,2,2);
% %     imshow(image);
% %     hold on;
% %     plot(x_blob, y_blob, 'r.', 'MarkerSize', 3);
% %     plot(x_range, y_range, 'r.', 'MarkerSize', 3);
% %     hold off;
% % 
% %     % Write frame to video
% %     frame = getframe(gcf);
% %     writeVideo(outputVideo4, frame);
% %     close(gcf);
% % end
% % 
% % close(outputVideo4);
% % fprintf('Video saved to: %s\n', fullfile(dataDir, 'BlobAreaAndHeightVsTime.avi'));
% 
% %% Create a video with 3 panels: Area vs. time, Height vs. time, Original blob image
% % disp('Now creating Video 5...');
% % 
% % outputVideo5 = VideoWriter(fullfile(dataDir, 'BlobArea_Height_Image.avi'));
% % outputVideo5.FrameRate = vidSpeed;
% % open(outputVideo5);
% % 
% % % Time, area, and height vectors
% % timeVecAll = data_nonTrim.timeVector;
% % areaVecAll = data_nonTrim.blob_area;
% % heightVecAll = data_nonTrim.blob_height;
% % iterVecAll = data_nonTrim.IterationNum_unit_;
% % 
% % validFrames = 1760:sampleInterval:3756;
% % 
% % for idx = validFrames
% %     rowIdx = find(data_nonTrim.IterationNum_unit_ == idx);
% %     if isempty(rowIdx) || ...
% %        isnan(data_nonTrim.blob_area(rowIdx)) || ...
% %        isnan(data_nonTrim.blob_height(rowIdx)) || ...
% %        isempty(data_nonTrim.x_blob{rowIdx})
% %         continue;
% %     end
% % 
% %     image = imread(data_nonTrim.image_path{rowIdx});
% %     x_blob = data_nonTrim.x_blob{rowIdx};
% %     y_blob = data_nonTrim.y_blob{rowIdx};
% %     x_range = data_nonTrim.x_range{rowIdx};
% %     y_range = data_nonTrim.y_range{rowIdx};
% % 
% %     % All data up to and including this frame
% %     time_so_far = timeVecAll(iterVecAll <= idx);
% %     area_so_far = areaVecAll(iterVecAll <= idx);
% %     height_so_far = heightVecAll(iterVecAll <= idx);
% % 
% %     % Create figure with 3 panels
% %     figure('Visible', 'off', 'Position', [100 100 1600 500]);
% % 
% %     % Left panel: Area vs. time
% %     subplot(1,3,1);
% %     plot(time_so_far, area_so_far, 'b.-', 'LineWidth', 1.5, 'MarkerSize', 10);
% %     xlabel('Time (s)');
% %     ylabel('Area (mm²)');
% %     grid on;
% %     xlim([200 max(timeVecAll)]);
% %     ylim([0 max(areaVecAll) * 1.1]);
% % 
% %     % Middle panel: Height vs. time
% %     subplot(1,3,2);
% %     plot(time_so_far, height_so_far, 'b*-');
% %     xlabel('Time (s)');
% %     ylabel('Height (mm)');
% %     grid on;
% %     xlim([200 max(timeVecAll)]);
% %     ylim([0 max(heightVecAll) * 1.1]);
% % 
% %     % Right panel: Original blob image
% %     subplot(1,3,3);
% %     imshow(image);
% %     hold on;
% %     plot(x_blob, y_blob, 'r.', 'MarkerSize', 3);
% %     plot(x_range, y_range, 'r.', 'MarkerSize', 3);
% %     hold off;
% % 
% %     % Write to video
% %     frame = getframe(gcf);
% %     writeVideo(outputVideo5, frame);
% %     close(gcf);
% % end
% % 
% % close(outputVideo5);
% % fprintf('Video saved to: %s\n', fullfile(dataDir, 'Area_Height_BlobOutline.avi'));
% % 



disp('Program done executing.');
%% Helper functions
function totalSeconds = mergeTimeData(secondsVec, nanosecondsVec)
    if any(secondsVec < 0) || any(nanosecondsVec < 0) || any(nanosecondsVec > 999999999)
        error('Input data contains values outside the valid range.');
    end
    totalSeconds = double(secondsVec) + double(nanosecondsVec) * 1e-9;
end
