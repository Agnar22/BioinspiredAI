#include "nsga.h"

/*
nsga::DominationIndividual::DominationIndividual() {
    DominationIndividual(0);
};
nsga::DominationIndividual::DominationIndividual(int ind) {
    ind = ind;
    times_dominated = 0;
};
*/

std::vector<std::vector<Individual>> nsga::sort_and_limit(std::vector<Individual> &pop, cv::Mat &img, int limit) {
    auto nondominated_sorted = fast_nondominated_sort(pop, img);
    int num_ind = 0;
    for (int i=0; i<nondominated_sorted.size(); ++i) {
        num_ind+=nondominated_sorted[i].size();
        if (num_ind>=limit) {
            nondominated_sorted = {nondominated_sorted.begin(), nondominated_sorted.begin()+i+1};
            nondominated_sorted[i] = crowding_sort(nondominated_sorted[i]);
            int ind_to_remove = num_ind-limit;
            nondominated_sorted[i] = {nondominated_sorted[i].begin(), nondominated_sorted[i].end()-ind_to_remove};
            break;
        }
    }
    return nondominated_sorted;
}

std::vector<std::vector<Individual>> nsga::fast_nondominated_sort(std::vector<Individual> &pop, cv::Mat &img, bool recalculate) {
    std::vector<DominationIndividual> dominations(pop.size());
    std::vector<std::vector<DominationIndividual>> fronts(1);
    for (int ind=0; ind<pop.size(); ++ind) {
        dominations[ind] = DominationIndividual(ind);
        if (recalculate)
            pop[ind].calculate_objectives(img);
    }
    
    for (DominationIndividual &p:dominations) {
        for (DominationIndividual &q:dominations) {
            if (pop[p.ind]<pop[q.ind]) {
                p.dominating.push_back(q.ind);
            } else if (pop[q.ind]<pop[p.ind]) {
                ++p.times_dominated;
            }
        }
        if (p.times_dominated==0)
            fronts[0].push_back(p);
    }

    for (int front=0; front<fronts.size(); ++front) {
        std::vector<DominationIndividual> next_front;
        for (DominationIndividual &p:fronts[front]) {
            for (int ind_q:p.dominating) {
                --dominations[ind_q].times_dominated;
                if (dominations[ind_q].times_dominated==0)
                    next_front.push_back(dominations[ind_q]);
            }
        }
        if (next_front.size()!=0)
            fronts.push_back(next_front);
    }

    std::vector<std::vector<Individual>> sorted;
    for (auto &front:fronts) {
        std::vector<Individual> ind_front;
        for (DominationIndividual &p:front) {
            ind_front.push_back(pop[p.ind]);
        }
        sorted.push_back(ind_front);
    }
    return sorted;
}

std::vector<Individual> nsga::crowding_sort(std::vector<Individual> &pop) {
    std::vector<CrowdingIndividual> crowding_pop(pop.size());
    for (int ind=0; ind<pop.size(); ++ind)
        crowding_pop[ind] = CrowdingIndividual(ind, std::vector<double>{pop[ind].edge_value, pop[ind].connectivity, pop[ind].overall_deviation});

    std::cout << std::endl;
    for (int sorting_pos=0; sorting_pos<3; ++sorting_pos) {
        auto sorting_function = [&](const CrowdingIndividual &l, const CrowdingIndividual &r){
            return l.values[sorting_pos] < r.values[sorting_pos];
        };
        std::sort(crowding_pop.begin(), crowding_pop.end(), sorting_function);
        double value_length = crowding_pop[crowding_pop.size()-1].values[sorting_pos] - crowding_pop[0].values[sorting_pos];
        for (int ind=0; ind<crowding_pop.size(); ++ind) {
            if (ind==0 || ind == pop.size()-1) {
                crowding_pop[ind].distances[sorting_pos] = INF;
                continue;
            }
            double neighbour_space = std::abs(crowding_pop[ind+1].values[sorting_pos] - crowding_pop[ind-1].values[sorting_pos]);
            crowding_pop[ind].distances[sorting_pos] = neighbour_space/value_length;
        }
    }

    std::sort(
        crowding_pop.begin(),
        crowding_pop.end(),
        [](const CrowdingIndividual &l, const CrowdingIndividual &r) {
            return std::accumulate(l.distances.begin(), l.distances.end(), 0.0) > std::accumulate(r.distances.begin(), r.distances.end(), 0.0);
        }
    );
    std::vector<Individual> crowding_sorted_individuals;
    for (CrowdingIndividual ind:crowding_pop) {
        crowding_sorted_individuals.push_back(pop[ind.ind]);
    }
    return crowding_sorted_individuals;
}
