#ifndef FILE_H
#define FILE_H

#include <vector>
#include <string>
#include <fstream>
#include <iterator>
#include <tuple>
#include <opencv2/opencv.hpp>

namespace file {
    cv::Mat read_image_to_vector(std::string);
}

#endif