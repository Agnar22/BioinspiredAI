#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

#include <vector>
#include <stdlib.h>
#include <opencv2/opencv.hpp>
#include "prim.h"
#include "kruskal.h"
#include "objective.h"

static double euc_dist(cv::Vec3b from, cv::Vec3b to) {
    return std::sqrt(
        std::pow(from[0]-to[0], 2) +
        std::pow(from[1]-to[1], 2) +
        std::pow(from[2]-to[2], 2)
    );
}

static double euc_dist(cv::Vec3b from, std::vector<double> to) {
    return std::sqrt(
        std::pow((double)(from[0])-to[0], 2) +
        std::pow((double)(from[1])-to[1], 2) +
        std::pow((double)(from[2])-to[2], 2)
    );
}

static double euc_dist(cv::Point3_<uchar> from, std::vector<double> to) {
    return std::sqrt(
        std::pow((double)(from.x)-to[0], 2) +
        std::pow((double)(from.y)-to[1], 2) +
        std::pow((double)(from.z)-to[2], 2)
    );
}

static double euc_dist(double x1, double y1, double x2, double y2) {
    return std::sqrt(
        std::pow(x1-x2, 2) +
        std::pow(y1-y2, 2)
    );
}

class Individual {
    public:
        Individual(cv::Mat, int);
        Individual(Individual&, Individual&, std::vector<int>&, int, int);
        Individual(Individual&, Individual&, int, cv::Mat&);
        std::vector<Dir> genes;
        std::vector<int> root;
        int width, height;
        double edge_value, connectivity, overall_deviation;

        void calculate_objectives(cv::Mat&);
        void mutate(int);
        void set_descendants(int, int);
        std::vector<int> find_children(int);
        int find_root(int);
        void find_roots();

    private:
        std::vector<bool> visited;
        void initialize_genes(cv::Mat, int);
        int root_search(int, int);
};

inline bool operator<(const Individual &l, const Individual &r) {
    // Edge value should be maximized, on the other hand connectivity and overall deviation should be minimized.
    if (l.edge_value>r.edge_value &&
        l.connectivity<r.connectivity &&
        l.overall_deviation<r.overall_deviation)
        return true;
    return false;
}

inline bool operator==(const Individual &l, const Individual &r) {
    // Edge value should be maximized, on the other hand connectivity and overall deviation should be minimized.
    if (l.edge_value==r.edge_value &&
        l.connectivity==r.connectivity &&
        l.overall_deviation==r.overall_deviation &&
        l.genes == r.genes &&
        l.root == r.root)
        return true;
    return false;
}

#endif