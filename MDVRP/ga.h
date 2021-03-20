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
        GA() {};
        GA(Problem&, int);
        void initialize_population(int);
        Individual get_individual(int);
        void simulate(int tourname_size, double, double, double, double, int);
        static std::vector<Individual> get_top_n(std::vector<Individual>&, int);
        static std::pair<Individual, Individual> tournament_selection(std::vector<Individual>, int, double);
        static void print_population_statistics(std::vector<Individual>&, std::string);
        static std::pair<Individual, Individual> best_cost_route_crossover(std::pair<Individual, Individual>&, Problem&);
        static void intra_depot_mutation(Individual&, double, double, double, Problem&);
        static void inter_depot_mutation(Individual&);
        static void mutate(Individual&, double, double, double, bool, Problem&);
};

#endif