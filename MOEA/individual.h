#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

#include <vector>
#include <stdlib.h>
#include <opencv2/opencv.hpp>
#include "prim.h"

static int euc_dist(cv::Vec3b from, cv::Vec3b to) {
    return std::sqrt(
        std::pow(from[0]-to[0], 2) +
        std::pow(from[1]-to[1], 2) +
        std::pow(from[2]-to[2], 2)
    );
}

class Individual {
    public:
        Individual(cv::Mat);
        std::vector<Dir> genes;
        std::vector<int> root;
        int width, height;

    private:
        void initialize_genes(cv::Mat);
        void find_roots();
        int find_root(int);
};

#endif
