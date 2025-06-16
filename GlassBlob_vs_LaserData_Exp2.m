%% GlassBlob_vs_LaserData_Exp2.m
% This program is used to process multiple frames.
clc; clear; close all;

%% User Input
expName = 'E:\SFF\Exp2_MeasuringTracks_DiffStartsStops\';
expCase = '3_ConstantLaserPower_50W_0o5mms_fp_7mm\';
imageDir = [expName, expCase, 'images\'];
substrate_start_row = 1350;
known_hypodermic_depth = 400; % px
one_px = 12.446e-6; % m

frameRange = 400:10:1000;
results = table('Size', [0 3], 'VariableTypes', {'double', 'cell', 'cell'}, ...
                'VariableNames', {'Frame', 'Z_vis', 'X_vis'});

for frame = frameRange
    fprintf('Processing frame %d...\n', frame);
    fullImagePath = fullfile(imageDir, sprintf('Acquisition%d.jpeg', frame));
    if ~isfile(fullImagePath)
        warning('Frame %d not found. Skipping...', frame);
        continue;
    end

    image = imread(fullImagePath);
    [image_height, image_width, ~] = size(image);

    substrate_mask = false(image_height, image_width);
    substrate_mask(substrate_start_row:end, :) = true;

    region_height = substrate_start_row - 1;
    region_width = 0.5 * (image_width - 200);

    lab_image = rgb2lab(image);
    a_channel = lab_image(:,:,2);
    b_channel = lab_image(:,:,3);

    region1_a = a_channel(1200:region_height, 1:region_width);
    region1_b = b_channel(1200:region_height, 1:region_width);
    region1_a2 = a_channel(800:region_height-100, image_width-region_width:image_width);
    region1_b2 = b_channel(800:region_height-100, image_width-region_width:image_width);

    region2_a = a_channel(1:500, 1:region_width);
    region2_b = b_channel(1:500, 1:region_width);
    region2_a2 = a_channel(1:500, image_width-region_width:image_width);
    region2_b2 = b_channel(1:500, image_width-region_width:image_width);

    bg1_a = mean([mean2(region1_a), mean2(region1_a2)]); 
    bg1_b = mean([mean2(region1_b), mean2(region1_b2)]);
    bg2_a = mean([mean2(region2_a), mean2(region2_a2)]); 
    bg2_b = mean([mean2(region2_b), mean2(region2_b2)]);

    distance_bg1 = sqrt((a_channel - bg1_a).^2 + (b_channel - bg1_b).^2);
    distance_bg2 = sqrt((a_channel - bg2_a).^2 + (b_channel - bg2_b).^2);
    threshold = 10;

    bg_mask_1 = distance_bg1 < threshold;
    bg_mask_2 = distance_bg2 < threshold;
    bg_mask_1(substrate_mask) = false;
    bg_mask_2(substrate_mask) = false;
    object_mask = ~bg_mask_1 & ~bg_mask_2 & ~substrate_mask;

    bw = uint8(object_mask) * 255;
    minSize = 2000;
    bw = bwareaopen(bw, minSize);
    bw = imfill(bw, 'holes');

    se = strel('disk', 3);
    bw_smooth = imopen(bw, se);
    bw_smooth = imclose(bw_smooth, se);
    bw_smooth = bwareaopen(bw_smooth, minSize);

    for row = 1:min(100, size(bw_smooth, 1))
        thisRow = bw_smooth(row, :);
        whiteIdx = find(thisRow);
        if ~isempty(whiteIdx)
            bw_smooth(row, whiteIdx(1):whiteIdx(end)) = 1;
        end
    end

    boundary = [];
    for col = 1:image_width
        row_idx = find(bw_smooth(:, col), 1, 'first');
        if ~isempty(row_idx)
            boundary = [boundary; col, row_idx];
        end
    end

    % Filter 1: remove ones that aren't continuous with exisiting cluster
    [sortedX, sortIdx] = sort(boundary(:,1), 'descend');
    sortedY = boundary(sortIdx, 2);
    dy = diff(sortedY);
    incontThreshold = 200;
    jumpIdx = find(abs(dy) > incontThreshold, 1, 'first');

    if ~isempty(jumpIdx)
        stopIdx = jumpIdx + (sortedY(jumpIdx + 1) > sortedY(jumpIdx));
        validIdx = sortIdx(1:stopIdx - 1);
        boundary = boundary(validIdx, :);
    end

    Z_vis = (substrate_start_row - boundary(:,2)); % px
    X_vis = boundary(:,1); % px

    % Filter 2: remove those belonging to work zone (unstable yet)
    validIdx = X_vis >= 850;
    X_vis = X_vis(validIdx);
    Z_vis = Z_vis(validIdx);

    % Pixel to metric conversion
    Z_vis = Z_vis * one_px * 1000; % mm
    X_vis = X_vis * one_px * 1000; % mm

    results = [results; {frame, Z_vis, X_vis}];
end

% Save to file
savePath = fullfile(expName, expCase, 'XZ_profile_visual_allFrames.mat');
save(savePath, 'results');

disp('Finished processing all frames.');
