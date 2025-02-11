%% Glass_BlobUp_vs_LaserData.m
clc; clear; close all;

%% Specify the dataset location here

dataDir = 'D:/ND Lab/';
expDir = 'BlobUp_Experiments/';
testCase = 'BlobUp_Exp2_600C_to_1000C';

dataDir = strcat(dataDir,expDir,testCase); % Path file for the dataset location
if ~isfolder(dataDir)
    error('Directory "%s" does not exist. Please check the path and try again.', dataDir);
end

clear expDir testCase

%% Read the csv file and create time vector here
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
data_nonTrim.IterationNum_unit_ =  data_nonTrim.IterationNum_unit_ + 1; % Increment Image counter by 1 to start from 0
data_nonTrim.iter_counter = data_nonTrim.IterationNum_unit_(1:end);

timeVector = totalStartTime - totalStartTime(1);
data_nonTrim.timeVector = timeVector;

% Clear all unneccesary variables to clean the workspace
clear startTimeSec startTimeNanosec endTimeProcessSec endTimeProcessNanosec endTimeLoopSec endTimeLoopNanosec totalStartTime totalEndTimeProcess totalEndTimeLoop timeVector

%% Read the image files and associate them with the time vector here
imagesFolder = fullfile(dataDir, 'images');

% Check if the 'images' folder exists
if ~isfolder(imagesFolder)
    error('Folder "images" does not exist in the specified directory.');
end

% Get a list of all JPEG files in the 'images' folder
imageFiles = dir(fullfile(imagesFolder, 'Acquisition*.jpeg'));

% Check if any JPEG images are found
if isempty(imageFiles)
    error('No JPEG images found inside the "images" folder.');
end

% Extract numbers from file names and sort them in ascending order
    imageNames = {imageFiles.name};
    imageNumbers = cellfun(@(x) sscanf(x, 'Acquisition%d.jpeg'), imageNames);
    [~, sortIdx] = sort(imageNumbers);

% Sort the file structure array by the image number
    imageFiles = imageFiles(sortIdx);
 
clear imagesFolder imageNumbers imageNames

% Associate the images with index in data_nonTrim 

% Initialize image path column with NaNs
data_nonTrim.image_path = repmat({NaN}, height(data_nonTrim), 1);

% Extract iteration numbers from existing images
imageNumbers = cellfun(@(x) sscanf(x, 'Acquisition%d.jpeg'), {imageFiles.name});
imagePaths = fullfile(dataDir, 'images', {imageFiles.name});

% Create a mapping of image numbers to their paths
imageMap = containers.Map(imageNumbers, imagePaths);

% Associate images with data_nonTrim table
for i = 1:height(data_nonTrim)
    iterNum = data_nonTrim.IterationNum_unit_(i);
    if isKey(imageMap, iterNum)
        data_nonTrim.image_path{i} = imageMap(iterNum); % Assign image path if exists
    end
end
 
clear imageNumbers imagePaths imageMap i iterNum

%% Visual processing here

% Define debug flag (set to true to visualize results)
debugFlag = true;

% Initialize vectors with NaN
blobArea = NaN(height(data_nonTrim), 1);
blobEdge = cell(height(data_nonTrim), 1);

for k = 1220:10:1580 % This is to select which iterations we want to test the algorithm on.
    
    imgPath = data_nonTrim.image_path{k};

    % Check if the image path is valid (i.e., not NaN)
    if ischar(imgPath) && exist(imgPath, 'file') % Ensure it's a string and file exists
        img = imread(imgPath); % Read the image 
        
        % Pre-processing
        gray_img = im2gray(img);
        bw = imbinarize(gray_img); 
        bw = imfill(bw, 'holes'); % Fill holes
        minSize = 1000;
        bw = bwareaopen(bw, minSize);

        % Disconnect blob from everything else:
        % Create a structuring element that is shaped like a disk with a radius of 15 px
        se = strel('disk', 15); 
        % Erode the non-circular looking objects
        bwEroded = imerode(bw, se);
        % Restore the eroded parts of the circle
        bwOpened = imdilate(bwEroded, se);

        % Get area and edge (convex hull) coordinates
        stats = regionprops(bwOpened, 'Area', 'ConvexHull');

        % Save blob data
        area = stats.Area;
        blobArea(k) = area;
        edge_coords = stats.ConvexHull;
        blobEdge{k} = edge_coords;

        % Debugging visualization
        if debugFlag
            
                if exist('figHandle', 'var') && isvalid(figHandle)
                clf(figHandle);  % Clear contents but keep the figure open
                else
                figHandle = figure(1);  % Create figure if it doesn’t exist
                hold on;
                end

            subplot(1, 2, 1);
            imshow(img);
            title('Original Image');

            % maybe another subplot for the dilated bw image
            
            subplot(1, 2, 2);
            imshow(img);
            hold on;
            plot(edge_coords(:,1), edge_coords(:,2), 'r-', 'LineWidth', 2);
            title(sprintf('Identified Blob Edge - Iteration %d', k));

            hold off;
        end
    end
end

% Store the results in the table
data_nonTrim.blobArea = blobArea;
data_nonTrim.blobEdge = blobEdge;

clear area bw bwEroded bwOpened edge_coords debugFlag gray_img k minSize se stats


% Plot Blob Area & Temperature vs Time
figure('Name', 'Fig 1: Blob Area & Temperature vs Time');
plot_x = data_nonTrim.timeVector;

yyaxis left
plot_y = data_nonTrim.blobArea;
plot(plot_x, plot_y, '-*', 'LineWidth', 1, 'DisplayName', 'Blob Area');
ylabel('Blob Area (pixel)');

yyaxis right
plot_y1 = data_nonTrim.ReferenceTemp_C_;
plot(plot_x, plot_y1, '--s', 'LineWidth', 1, 'DisplayName', 'Ref. Temp.'); % Dashed line

hold on;
plot_y2 = data_nonTrim.OptrisMaxTemp_C_;
plot(plot_x, plot_y2, '-', 'LineWidth', 1, 'DisplayName', 'Optris Max Temp.'); % Solid line

ylabel('Temperature (C)');  

% Set common properties
xlabel('Time (sec)');
grid on;
set(gca, 'FontName', 'Arial', 'FontSize', 14, 'FontWeight', 'Bold');
legend();
hold off;

% Plot Blob Area & Laser vs Time
figure('Name', 'Fig 2: Blob Area & Laser Power vs Time');
plot_x = data_nonTrim.timeVector;

yyaxis left
plot_y = data_nonTrim.blobArea;
plot(plot_x, plot_y, '-*', 'LineWidth', 1, 'DisplayName', 'Blob Area');
ylabel('Blob Area (pixel)');

yyaxis right
plot_y1 = data_nonTrim.LaserPowerCommanded_NoLim_W_;
plot(plot_x, plot_y1, '--', 'Color', [1, 0, 0], 'LineWidth', 1, 'DisplayName', 'Laser Power Commanded'); % Bright Red (RGB: [1, 0, 0])

hold on;
plot_y2 = data_nonTrim.LaserPowerFbk_W_;
plot(plot_x, plot_y2, '-', 'Color', [0.6, 0, 0], 'LineWidth', 1, 'DisplayName', 'Laser Power Measured'); % Dark Red (RGB: [0.6, 0, 0])

ylabel('Laser Power (W)');  

% Set common properties
xlabel('Time (sec)');
grid on;
set(gca, 'FontName', 'Arial', 'FontSize', 14, 'FontWeight', 'Bold');
legend();
hold off;

disp('Program finished successfully');

%% Helper functions for this program

% Function to merge time values and create time vector
function totalSeconds = mergeTimeData(secondsVec, nanosecondsVec)

    if any(secondsVec < 0) || any(nanosecondsVec < 0) || any(nanosecondsVec > 999999999)
        error('Input data contains values outside the valid range.');
    end

    nanosecondsConverted = double(nanosecondsVec) * 1e-9;
    totalSeconds = double(secondsVec) + nanosecondsConverted;
end