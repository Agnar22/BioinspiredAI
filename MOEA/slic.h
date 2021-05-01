#ifndef SLIC_H
#define SLIC_H

#include "individual.h"
#include <vector>
#include <opencv2/opencv.hpp>

namespace slic {
    struct RGBCluster {
        private:
            std::vector<double> sum_rgb;
            double sum_x, sum_y;
            int num_pxls;

        public: 
            std::vector<double> rgb;
            double x, y;

            RGBCluster(): rgb{0,0,0}, sum_rgb{0, 0, 0}, x(0), sum_x(0), y(0), sum_y(0), num_pxls(0) {};
            RGBCluster(double x, double y): rgb{0,0,0}, sum_rgb{0, 0, 0}, x(x), sum_x(0), y(y), sum_y(0), num_pxls(0) {};
            RGBCluster(double r, double g, double b, double x, double y): rgb{r,g,b}, sum_rgb{0, 0, 0}, x(x), sum_x(0), y(y), sum_y(0), num_pxls(0) {};

            //void add_pxl(int inp_x, int inp_y, cv::Vec3b &clr) {
            void add_pxl(int inp_x, int inp_y, cv::Point3_<uchar> &clr) {
                ++num_pxls;
                sum_x += inp_x;
                x = sum_x/(double)(num_pxls);
                sum_y += inp_y;
                y = sum_y/(double)(num_pxls);
                /*
                for (int clr_pos=0; clr_pos<3; ++clr_pos) {
                    sum_rgb[clr_pos] += (double)(clr[clr_pos]);
                    rgb[clr_pos] = sum_rgb[clr_pos]/(double)(num_pxls);
                }
                */
                std::cout << (int)(clr.x) << std::endl;
                sum_rgb[0] += (double)(clr.x);
                sum_rgb[1] += (double)(clr.y);
                sum_rgb[2] += (double)(clr.z);
                rgb[0] = sum_rgb[0]/(double)(num_pxls);
                rgb[1] = sum_rgb[1]/(double)(num_pxls);
                rgb[2] = sum_rgb[2]/(double)(num_pxls);
            }
    };

    std::vector<int> slic(cv::Mat&, int);
    std::vector<int> get_neighbours(int, int, int);
    std::vector<int> enforce_connectivity(std::vector<int>, std::vector<RGBCluster>, int, int);
}

#endif