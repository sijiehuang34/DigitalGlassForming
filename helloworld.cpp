#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>
#include <opencv2/opencv.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/core/hal/interface.h>
#include <filesystem>
#include <cmath>

#define CV_32S 4

using namespace cv;

void clearTerminal() {
    system("cls");  // Windows
}

int main(int argc, char** argv) {
    std::string imagePath = "C:/Users/Cindy Huang/projects/fail_visCam_Exp8_frames/frame_30.png";
    std::string outputFolder = "C:/Users/Cindy Huang/projects/fail_visCam_Exp8_denoised/";

    cv::Mat binary = cv::imread(imagePath, cv::IMREAD_GRAYSCALE);
    if (binary.empty()) {
        std::cerr << "Error: Could not load the image at " << imagePath << std::endl;
        return -1;
    }

    // Perform noise removal
    cv::Mat labels;
    int nLabels = cv::connectedComponents(binary, labels);
    cv::Mat denoised = cv::Mat::zeros(binary.size(), CV_8UC1);

    for (int label = 1; label < nLabels; ++label) {
        cv::Mat componentMask = (labels == label);
        int componentArea = cv::countNonZero(componentMask);
        if (componentArea > 3000) {
            denoised.setTo(255, componentMask);
        }
    }

    // Define a vector to store the points of the left edge
    std::vector<cv::Point> leftEdgePoints;

    // Process the top 50% of rows
    int rowsToProcess = denoised.rows / 2;
    int midCol = denoised.cols / 2;

    for (int row = 0; row < rowsToProcess; ++row) {
        for (int col = 0; col < midCol; ++col) { // Start from the left end
            if (denoised.at<uchar>(row, col) == 255) {
                leftEdgePoints.emplace_back(col, row); // Save the point
                break; // Stop searching in this row after the first white point is found
            }
        }
    }

    // Create a color image to mark the points
    cv::Mat colorImage;
    cv::cvtColor(denoised, colorImage, cv::COLOR_GRAY2BGR);

    // Mark the points with red color
    for (const auto& point : leftEdgePoints) {
        cv::circle(colorImage, point, 2, cv::Scalar(0, 0, 255), -1); // Red color
    }

    // Calculate the angle from the horizontal using the left edge points
    double angle = 0.0;
    if (leftEdgePoints.size() > 1) {
        cv::Point firstPoint = leftEdgePoints.front();
        cv::Point lastPoint = leftEdgePoints.back();

        double deltaY = static_cast<double>(lastPoint.y - firstPoint.y);
        double deltaX = static_cast<double>(lastPoint.x - firstPoint.x);

        angle = std::atan2(deltaY, deltaX) * 180.0 / CV_PI; // Angle in degrees
    } else {
        std::cerr << "Not enough points to calculate the angle." << std::endl;
    }

    // Display the angle on the bottom-right corner of the image
    std::string angleText = "Angle: " + std::to_string(angle) + " degrees";
    int fontFace = cv::FONT_HERSHEY_SIMPLEX;
    double fontScale = 0.7;
    int thickness = 2;

    // Calculate the position for the text
    int baseline = 0;
    cv::Size textSize = cv::getTextSize(angleText, fontFace, fontScale, thickness, &baseline);
    cv::Point textOrg(colorImage.cols - textSize.width - 10, colorImage.rows - 10); // Bottom-right corner

    // Put the text on the image
    cv::putText(colorImage, angleText, textOrg, fontFace, fontScale, cv::Scalar(0, 255, 0), thickness);

    // Save the result to the output folder
    std::string outputPath = outputFolder + "processed_image.png";
    if (!cv::imwrite(outputPath, colorImage)) {
        std::cerr << "Error: Could not save the image to " << outputPath << std::endl;
        return -1;
    }

    std::cout << "Processed image saved to: " << outputPath << std::endl;

    return 0;
}



// int main(int argc, char** argv){
//     ///////////////////////////////// Load Video ////////////////////////////////////////////
//     // Define the folder path
//     std::string folderPath = "C:/Users/Cindy Huang/projects/processed_frames/";

//     // Open the video file
//     cv::VideoCapture video("C:/Users/Cindy Huang/projects/fail_visCam_Exp8.mp4");
//     if (!video.isOpened()) {
//         std::cerr << "Error: Could not open the video file!!!!!" << std::endl;
//         return -1;
//     }

//     cv::Mat frame;
//     int frameCount = 0;

//     // Loop through the video frames
//     while (video.read(frame)) {
//         // Check if it's the 10th frame
//         if (frameCount % 10 == 0) {

//             ////////////////////////////////// Image Processing ////////////////////////////
//             // Calculate the cropped region
//             int cropOffset = frame.cols / 6; // Assumes the image is square
//             int newSize = frame.cols - 2 * cropOffset;

//             // Define the region of interest (ROI)
//             cv::Rect roi(cropOffset, cropOffset, newSize, newSize);

//             // Crop the image
//             cv::Mat croppedImage = frame(roi);
            
//             // Create a Mat to store the grayscale image
//             cv::Mat grayImage;

//             // Convert the color image to grayscale
//             cv::cvtColor(croppedImage, grayImage, cv::COLOR_BGR2GRAY);

//             Mat blurred;
//             GaussianBlur(grayImage, blurred, Size(5, 5), 0);

//             Mat binary;
//             adaptiveThreshold(blurred, binary, 255, ADAPTIVE_THRESH_GAUSSIAN_C, THRESH_BINARY, 501, -1); 
//             // Block Size 501: determines the size of the region used for calculating the threshold
//             // ^^^At least 50% of column/row size --> 501 (this shouldn't need to change unless focal distance changes)
//             // Constant: subtracted from the weighted mean (fine-tuning)
//             // ^^^Keep it <= 0 to keep the white area grainy (easier noise reduction by connected area)

//             // //////////////////////////////// Noise reduction ///////////////////////////////////////
//             // Mat labels;
//             // int nLabels = connectedComponents(binary, labels);
//             // Mat denoised = Mat::zeros(binary.size(), CV_8UC1);
//             // for (int label = 1; label < nLabels; ++label) {
//             //     Mat componentMask = (labels == label);
//             //     // Calculate the area (number of pixels) of the component
//             //     int componentArea = countNonZero(componentMask);
//             //     // Remove components that are smaller than a defined threshold (considered noise)
//             //     if (componentArea > 3000) { 
//             //         denoised.setTo(255, componentMask);
//             //     }
//             // }

//             ////////////////////////////////// Image Saving ////////////////////////////////
//             // Construct the filename for each saved frame
//             std::string filename = folderPath + "frame_" + std::to_string(frameCount) + ".png";
//             // Save the frame
//             if (cv::imwrite(filename, binary)) {
//                 std::cout << "Saved frame " << frameCount << " to " << filename << std::endl;
//             } else {
//                 std::cerr << "Error: Could not save frame!!!!! " << frameCount << std::endl;
//             }
//         }
//         // Increment the frame counter
//         frameCount++;
//     }

//     // Release the video capture object
//     video.release();
//     std::cout << "Processing complete!!!!!" << std::endl;
//     /////////////////////////////////////////////////////////////////////////////////////////
//     return 0;
// }
