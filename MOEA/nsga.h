#ifndef NSGA_H
#define NSGA_H

#include <vector>
#include <opencv2/opencv.hpp>
#include <algorithm>
#include "individual.h"

namespace nsga {
    struct DominationIndividual {
        int ind;
        int times_dominated;
        std::vector<int> dominating;

        DominationIndividual(int n): ind(n), times_dominated(0) {};
        DominationIndividual(): ind(0), times_dominated(0) {};
    };

    struct CrowdingIndividual {
        int ind;
        //double edge_distance, connectivity_distance, overall_deviation_distance;
        std::vector<double> values;
        std::vector<double> distances;
        CrowdingIndividual(int n, std::vector<double> values): ind(n), values{values}, distances(3, 0.0) {};
        CrowdingIndividual(): ind(0), values(3, 0.0), distances(3, 0.0) {};
    };

    std::vector<std::vector<Individual>> sort_and_limit(std::vector<Individual>&, cv::Mat&, int);
    std::vector<std::vector<Individual>> fast_nondominated_sort(std::vector<Individual>&, cv::Mat&, bool recalculate=true);
    std::vector<Individual> crowding_sort(std::vector<Individual>&);
}

#endif