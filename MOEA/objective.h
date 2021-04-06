#ifndef OBJECTIVE_H
#define OBJECTIVE_h

#include <vector>
#include <opencv2/opencv.hpp>
#include <map>
#include "dir.h"
#include "individual.h"

namespace obj {
    double edge_value(cv::Mat&, std::vector<Dir>&, std::vector<int>&, int, int);
    double connectivity(std::vector<int>&, int, int);
    std::map<int, cv::Vec3d> find_segment_means(cv::Mat&, std::vector<int> &);
    double overall_deviation(cv::Mat&, std::vector<int>&, int, int);
}

#endif