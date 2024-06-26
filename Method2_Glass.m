%% Method1_ThreeCases.m
clc; clear; close all;

% Exp 1: Standard Warm (triangle1_frame_1800)
% Exp 2: Messy (triangle2_frame_1700)
% Exp 3: Curly (triangle3_frame_2600)
% Exp 4: Barely Connecting (triangle4_frame_2600)
% Exp 5: Standard Cool (triangle5_frame_1400)
% Exp 6:
% Exp 7: Standard Cool (triangle7_frame_1300)
% Exp 8: 
img = imread("frame_1800.jpg");

%--------------------------------------------------------------------------

% Turn into grayscale for pixel histogram analysis
gray_img = im2gray(img);
% Binarize (retains top side of the punch only)
bw = imbinarize(gray_img,"adaptive"); 


% Denoise (removes small objects)
minSize = 2100;
bw = bwareaopen(bw,minSize);

% Make substrate region = 0 px
bw(1810:end, 1:end) = 0;
bw(1785:end, 2300:end) = 0;
% Make random blob = 0 px
bw(1000:1250, 2200:end) = 0;
% Make nozzle region = 0 px
bw(1:1175, 1:1224) = 0;

% Denoise again
minSize = 150;
bw = bwareaopen(bw,minSize);

% Overlay on processed image
imshow(bw)
% % Overlay on original image
% imshow(img)

%--------------------------------------------------------------------------

hold on;

% Find highest & lowest points along the glass structure
highest_Points = []; % initialize a matrix to store the coords
lowest_Points = [];

totalDim = size(bw);
totalRow = totalDim(1);
totalCol = totalDim(2);

%--------------------------------------------------------------------------

% Find the highest points
for iCol = 1:totalCol
    for iRow = 1:totalRow
        if bw(iRow,iCol) == 1
            target = [iCol,iRow];
            highest_Points = [highest_Points; target];
            % exit the inner loop immediately
            % disregard the rest of the rows
            break;
        end
    end
end

% Find the lowest points
for iCol = 1:totalCol
    for iRow = totalRow:-1:1
        if bw(iRow,iCol) == 1
            target = [iCol,iRow];
            lowest_Points = [lowest_Points; target];
            % exit the inner loop immediately
            % disregard the rest of the rows
            break;
        end
    end
end

%--------------------------------------------------------------------------

% % Make the background white
% bw(bw == 0) = 1;
% imshow(bw)
% hold on;

%--------------------------------------------------------------------------

% Initialize continuity check
continuityCheck = true;

% Initialize loop index
i = 120;

while i <= length(highest_Points)
    
    % Plot continuous points
    if continuityCheck == true
        plot(highest_Points(i,1), highest_Points(i,2), 'r.', 'MarkerSize', 10);
    end

    % Check continuity for next point
    if abs(highest_Points(i,2) - highest_Points(i+20,2)) < 50
        continuityCheck = true;
    else
        continuityCheck = false;
        % Find the next continuous point
        discont_Pt = i + 20;
        while discont_Pt < length(highest_Points)-20 && abs(highest_Points(discont_Pt,2) - highest_Points(discont_Pt+20,2)) < 50
            discont_Pt = discont_Pt + 20;
        end
        % Redefine i to be the next correct/continuous point and skip the points in between
        nextCorrectPt = discont_Pt+20;
        i = nextCorrectPt;
        continuityCheck = true;
        continue; % Skip the increment of i at the end of the loop
    end
    
    % Increment loop index
    i = i + 20;
end

%--------------------------------------------------------------------------

% Initialize continuity check
continuityCheck = true;

% Initialize loop index
i = 1;

while i <= length(lowest_Points)-20
    
    % Plot continuous points
    if continuityCheck == true
        plot(lowest_Points(i,1), lowest_Points(i,2), 'r.', 'MarkerSize', 10);
    end

    % Check continuity for next point
    if abs(lowest_Points(i,2) - lowest_Points(i+20,2)) < 50
        continuityCheck = true;
    else
        continuityCheck = false;
        % Find the next continuous point
        discont_Pt = i + 20;
        while discont_Pt < length(lowest_Points)-20 && abs(lowest_Points(discont_Pt,2) - lowest_Points(discont_Pt+20,2)) < 50
            discont_Pt = discont_Pt + 20;
        end
        % Redefine i to be the next correct/continuous point and skip the points in between
        nextCorrectPt = discont_Pt+20;
        i = nextCorrectPt;
        continuityCheck = true;
        continue; % Skip the increment of i at the end of the loop
    end
    
    % Increment loop index
    i = i + 20;
end

%--------------------------------------------------------------------------

% Find appropriate points

% Top vertex 
x1 = 1750; % manual
y1 = 1150; % manual
plot(x1, y1, 'gx', 'MarkerSize', 10, "LineWidth", 2);

% Left vertex
pxConv = 90; % 1mm = 80px (Images saved is shrunk by approximately 2x)
triangleDim = 7 * pxConv; 
x2 = x1 - triangleDim;
y2 = y1 + triangleDim;
plot(x2, y2, 'gx', 'MarkerSize', 10, "LineWidth", 2);

% Right vertex
x3 = x1 + triangleDim;
y3 = y1 + triangleDim;
plot(x3, y3, 'gx', 'MarkerSize', 10, "LineWidth", 2);

% Left tail
tailLength = 2 * triangleDim / 7; % according to the 2mm : 7mm ratio
x4 = x2 - tailLength;
y4 = y2;
plot(x4, y4, 'gx', 'MarkerSize', 10, "LineWidth", 2);

% Right tail
x5 = x3 + tailLength;
y5 = y3;
plot(x5, y5, 'gx', 'MarkerSize', 10, "LineWidth", 2);

%--------------------------------------------------------------------------

% % Line up the vertices
% 
% % Left line
% xArray1 = [x1, x2];
% yArray1 = [y1, y2];
% line(xArray1, yArray1, 'Color', 'g', 'LineWidth', 2);
% 
% % Right line
% xArray2 = [x1, x3];
% yArray2 = [y1, y3];
% line(xArray2, yArray2, 'Color', 'g', 'LineWidth', 2);
% 
% % Left tailing line
% xArray3 = [x2, x4];
% yArray3 = [y2, y4];
% line(xArray3, yArray3, 'Color', 'g', 'LineWidth', 2);
% 
% % Right tailing line
% xArray4 = [x3, x5];
% yArray4 = [y3, y5];
% line(xArray4, yArray4, 'Color', 'g', 'LineWidth', 2);

hold off;
