#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

#include <vector>
#include <stdlib.h>
#include <opencv2/opencv.hpp>
#include "prim.h"
#include "objective.h"

static double euc_dist(cv::Vec3b from, cv::Vec3b to) {
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
        double edge_value, connectivity, overall_deviation;

        void calculate_objectives(cv::Mat&);

    private:
        void initialize_genes(cv::Mat);
        void find_roots();
        int find_root(int);
};

inline bool operator<(const Individual &l, const Individual &r) {
    // Edge value should be maximized, on the other hand connectivity and overall deviation should be minimized.
    if (l.edge_value>r.edge_value &&
        l.connectivity<r.connectivity &&
        l.overall_deviation<r.overall_deviation)
        return true;
    return false;
}

#endif