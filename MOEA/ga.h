#ifndef GA_H
#define GA_H

#include <vector>
#include <opencv2/opencv.hpp>
#include "individual.h"
#include "nsga.h"
#include "config.h"


class GA {
    public:
        std::vector<Individual> population;
        bool nsga_ii;
        cv::Mat image;


        GA(int, bool, cv::Mat);
        void simulate();

    private:
        std::pair<std::pair<int, int>, std::pair<int, int>> binary_tournament_selection(std::vector<std::vector<Individual>>&);
        std::pair<int, int> select_parent_pos(std::vector<std::vector<Individual>>&);
        std::pair<Individual, Individual> crossover(Individual&, Individual&);
        void mutate(Individual&);

};

template<typename T>
std::vector<T> flatten(const std::vector<std::vector<T>> &orig) {   
    std::vector<T> ret;
    for(const auto &v: orig)
        ret.insert(ret.end(), v.begin(), v.end());
    return ret;
}

#endif
