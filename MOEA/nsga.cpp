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

std::vector<std::vector<Individual>> nsga::fast_nondominated_sort(std::vector<Individual> &pop, cv::Mat &img) {
    std::vector<DominationIndividual> dominations(pop.size());
    std::vector<std::vector<DominationIndividual>> fronts(1);
    for (int ind=0; ind<pop.size(); ++ind) {
        dominations[ind] = DominationIndividual(ind);
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

