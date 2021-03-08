#ifndef GA_H
#define GA_H

#include <vector>
#include <queue>
#include "individual.h"
#include "problem.h"

class GA {
    private:
        std::vector<Individual> population;
        Problem pr;

    public:
        GA(Problem&, int);
        void initialize_population(int);
        void simulate(int tourname_size, double, double, double, double, int);
        static std::vector<Individual> get_top_n(std::vector<Individual>&, int);
        static std::pair<Individual, Individual> tournament_selection(std::vector<Individual>, int, double);
        static std::pair<Individual, Individual> best_cost_route_crossover(std::pair<Individual, Individual>&);
        void test_fitness();
        static void intra_depot_mutation(Individual&, double, double, double);
        static void inter_depot_mutation(Individual&);
        static void mutate(Individual&, double, double, double, bool);
};

#endif