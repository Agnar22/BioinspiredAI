#include "file.h"

namespace file {
    cv::Mat read_image_to_vector(std::string file_name) {
        // FIXME: For some reason, the third channel is always 2 more than numpy.
        cv::Mat image = cv::imread(file_name, cv::IMREAD_COLOR);
        /*
        for (int row=0; row < 1; ++row) {
            for (int col=0; col < 10; ++col) {
                std::cout << (int)image.at<cv::Vec3b>(row, col).val[0] << " ";
                std::cout << (int)image.at<cv::Vec3b>(row, col).val[1] << " ";
                std::cout << (int)image.at<cv::Vec3b>(row, col).val[2] << ", ";
            }
            std::cout << std::endl;
        }
        */
        return image;
    }
}