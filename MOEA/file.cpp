#include "file.h"

namespace file {
    cv::Mat read_image_to_vector(std::string file_name) {
        // FIXME: For some reason, the third channel is always 2 more than numpy.
        cv::Mat image = cv::imread(file_name, cv::IMREAD_COLOR);
        return image;
    }
}