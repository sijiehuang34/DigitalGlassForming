%% VidGen_OverHalfExpVid.m
clc; clear; close all;

%% Video Writers
vidOut1 = VideoWriter('width_vs_depth_blob_outline.avi');
vidOut2 = VideoWriter('height_vs_frame_blob_outline.avi');
vidOut1.FrameRate = 1;
vidOut2.FrameRate = 1;
open(vidOut1);
open(vidOut2);

% Initialize data
height_values = [];
frame_numbers = [];

%% Loop over frames
for frame_idx = 1760:1:3756
    fprintf('Processing frame %d...\n', frame_idx);

    filename = sprintf('E:\\3_PowerRamp_1W\\3_PowerRamp_1W\\images\\Acquisition%d.jpeg', frame_idx);
    if ~isfile(filename)
        fprintf('File not found: %s\n', filename);
        continue;
    end

    image = imread(filename);
    [image_height, image_width, ~] = size(image);

    % Define substrate
    finetune_val = 35;
    substrate_start_row = floor(image_height * (5/6)) + 1 + finetune_val;
    substrate_mask = false(image_height, image_width);
    substrate_mask(substrate_start_row:end, :) = true;
    one_px = 12.446e-6; % px to m

    % Convert to LAB color space
    lab_image = rgb2lab(image);
    a_channel = lab_image(:,:,2);
    b_channel = lab_image(:,:,3);

    region_height = substrate_start_row - 1;
    region_width = 400;
    region1_a = a_channel(1:region_height, 1:region_width);
    region1_b = b_channel(1:region_height, 1:region_width);
    region2_a = a_channel(1:region_height, end-region_width:end);
    region2_b = b_channel(1:region_height, end-region_width:end);

    background_a = mean([mean2(region1_a), mean2(region2_a)]);
    background_b = mean([mean2(region1_b), mean2(region2_b)]);
    distance = sqrt((a_channel - background_a).^2 + (b_channel - background_b).^2);
    threshold = 10;
    background_mask = distance < threshold;
    background_mask(substrate_mask) = false;

    object_mask = ~background_mask & ~substrate_mask;
    bw = uint8(object_mask) * 255;
    bw = bwareaopen(bw, 1000);
    bw = imfill(bw, 'holes');
    se = strel('disk', 3);
    bw_smooth = imopen(bw, se);
    bw_smooth = imclose(bw_smooth, se);
    bw_smooth = bwareaopen(bw_smooth, 1000);

    % Find boundaries
    [Boundary, ~, ~] = bwboundaries(bw_smooth, 4, 'noholes');
    [~, sortIdx] = sort(cellfun(@(x) size(x,1), Boundary), 'descend');
    Boundary = Boundary(sortIdx);
    Boundary = Boundary(cellfun(@(x) size(x,1) >= 2000, Boundary));
    Boundary = Boundary(~cellfun(@(x) any(x(:,2)==1) | any(x(:,2)==image_width), Boundary));
    if isempty(Boundary)
        continue;
    end
    boundaryPoints = Boundary{1};

    % Width vs. Depth
    [sortedY, sortIdx] = sort(boundaryPoints(:,1), 'ascend');
    sortedX = boundaryPoints(sortIdx, 2);
    boundaryPoints = [sortedY, sortedX];
    x_all = boundaryPoints(:,2);
    y_all = boundaryPoints(:,1);
    uniqueY = unique(boundaryPoints(:,1), 'sorted');
    n = length(uniqueY);
    wid_vs_dep = zeros(n, 2);
    wid_vs_dep(:,1) = uniqueY;

    for i = 1:n
        yVal = uniqueY(i);
        xVals = boundaryPoints(boundaryPoints(:,1) == yVal, 2);
        width = abs(max(xVals) - min(xVals));
        wid_vs_dep(i,2) = width;
    end

    depth = wid_vs_dep(:,1);
    width = wid_vs_dep(:,2);
    idx1 = (depth >= 0) & (depth <= 500);
    y1_mean = mean(width(idx1));
    w_hypodermic = 1.7e-3 / one_px;
    w_filament = 1.4e-3 / one_px;
    w_diff = w_hypodermic - w_filament;
    y2_mean = y1_mean - w_diff;

    valid_idx = find(depth > 600 & width > y1_mean, 1, 'first');
    if isempty(valid_idx)
        continue;
    end
    rough_blob_start = depth(valid_idx);
    idx_curve = (depth >= rough_blob_start) & (depth <= max(depth));
    x_curve = depth(idx_curve);
    y_curve = width(idx_curve);
    p = polyfit(x_curve, y_curve, 2);
    x_extended = linspace(min(depth), max(depth)*1.2, 1000);
    y_fit = polyval(p, x_extended);
    tp = roots([p(1), p(2), p(3) - y2_mean]);
    tp = tp(imag(tp)==0);
    if isempty(tp)
        continue;
    end
    tp_ref = min(tp);

    % Convert units
    depth = depth * one_px * 1000;
    width = width * one_px * 1000;
    y1_mean = y1_mean * one_px * 1000;
    y2_mean = y2_mean * one_px * 1000;
    x_extended = x_extended * one_px * 1000;
    y_fit = y_fit * one_px * 1000;
    tp_ref = tp_ref * one_px * 1000;

    % Height above substrate
    blob_substrate_dist = substrate_start_row - max(wid_vs_dep(:,1));
    blob_substrate_dist = blob_substrate_dist * one_px * 1000;
    height_values(end+1) = blob_substrate_dist;
    frame_numbers(end+1) = frame_idx;

    % Filter blob points
    blob_outline = boundaryPoints(boundaryPoints(:,1) > tp_ref / (one_px * 1000), :);
    x_blob = blob_outline(:,2);
    y_blob = blob_outline(:,1);
    y_opening = min(y_blob);
    x_opening = x_blob(y_blob == y_opening);
    x_opening_start = min(x_opening);
    x_opening_end = max(x_opening);
    x_range = x_opening_start:x_opening_end;
    y_range = y_opening * ones(size(x_range));

    %% Video 1 - Width vs. Depth + Blob Outline
    f1 = figure('Visible','off'); tiledlayout(1,2);
    nexttile;
    plot(depth, width, 'k.', 'MarkerSize', 5);
    hold on;
    plot(x_extended, y_fit, 'g-', 'LineWidth', 2);
    yline(y1_mean, 'b-', 'LineWidth', 1.5);
    yline(y2_mean, 'b-', 'LineWidth', 1.5);
    plot(tp_ref, y2_mean, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
    xlabel('Depth (mm)'); ylabel('Width (mm)');
    title(sprintf('Frame %d', frame_idx)); grid on;
    xlim([0 20]); ylim([0 3]);

    nexttile;
    imshow(image); hold on;
    plot(x_blob, y_blob, 'r.', 'MarkerSize', 3);
    plot(x_range, y_range, 'r.', 'MarkerSize', 3);
    title("Blob Outlined");

    writeVideo(vidOut1, getframe(f1));
    close(f1);

    %% Video 2 - Height vs Frame + Blob Outline
    f2 = figure('Visible','off'); tiledlayout(1,2);
    nexttile;
    plot(frame_numbers, height_values, 'b.-');
    xlabel('Frame'); ylabel('Height (mm)');
    title('Blob Height Over Time');
    xlim([1760 3756]); ylim([0 10]); grid on;

    nexttile;
    imshow(image); hold on;
    plot(x_blob, y_blob, 'r.', 'MarkerSize', 3);
    plot(x_range, y_range, 'r.', 'MarkerSize', 3);
    title("Blob Outlined");

    writeVideo(vidOut2, getframe(f2));
    close(f2);
end

% Close videos
close(vidOut1);
close(vidOut2);

fprintf('All frames processed. Videos saved.\n');
