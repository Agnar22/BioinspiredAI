#ifndef NSGA_H
#define NSGA_H

#include <vector>
#include <opencv2/opencv.hpp>
#include "individual.h"

namespace nsga {
    struct DominationIndividual {
        int ind;
        int times_dominated;
        std::vector<int> dominating;

        DominationIndividual(int n): ind(n), times_dominated(0) {};
        DominationIndividual(): ind(0), times_dominated(0) {};
    };

    std::vector<std::vector<Individual>> fast_nondominated_sort(std::vector<Individual>&, cv::Mat&);
}

#endif