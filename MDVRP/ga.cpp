#include "ga.h"


GA::GA(Problem &pr, int population_size): pr{pr} {
    initialize_population(population_size);
}

void GA::initialize_population(int population_size){
    for (int i=0; i<population_size; ++i)
        population.push_back(Individual (pr));
}

void GA::simulate(int tourname_size, double stoch_tournament_prob, double prob_rev_mut, double prob_re_routing, double prob_swapping, int inter_depot_swapping) {
    for (int gen=0; gen<500; ++gen) {
        std::vector<Individual> child_gen;
        int num_elites = (int)((double)(population.size())/100.0);
        while (child_gen.size()<population.size()-num_elites) {
            // Selection.
            std::pair<Individual, Individual> parents = GA::tournament_selection(population, tourname_size, stoch_tournament_prob);

            // Recombination.
            std::pair<Individual, Individual> children = GA::best_cost_route_crossover(parents);

            // Mutation.
            GA::mutate(children.first, prob_rev_mut, prob_re_routing, prob_swapping, gen%inter_depot_swapping==0);
            GA::mutate(children.second, prob_rev_mut, prob_re_routing, prob_swapping, gen%inter_depot_swapping==0);

            // Acceptance.
            child_gen.push_back(children.first);
            child_gen.push_back(children.second);
        }

        // Elitism.
        std::vector<Individual> best_parents = GA::get_top_n(population, num_elites);
        child_gen.insert(child_gen.end(), best_parents.begin(), best_parents.end());

        // Next generation.
        population=child_gen;
    }
}

std::vector<Individual> GA::get_top_n(std::vector<Individual> &pop, int n) {
    std::priority_queue<Individual> pq;
    std::vector<Individual> inds;
    for (Individual ind:pop)
        pq.push(ind);
    for (int pos=0; pos<n; ++pos) {
        inds.push_back(pq.top());
        pq.pop();
    }
    return inds;
}

std::pair<Individual, Individual> GA::tournament_selection(std::vector<Individual> population, int tournament_size, double stoch){
    /**
     * Tournament selection with tournament_size number of candidates. Selects the two most fit individuals with probability 1-stoch, else it randomly chooses.
     * Remark: The parents might be equal.
     */
    std::vector<int> tournament_set(tournament_size);
    std::generate(tournament_set.begin(), tournament_set.end(), rand);
    if (stoch < (double)(rand()) / (double)(RAND_MAX)) {
        int parent1=tournament_set[rand()%tournament_size];
        int parent2=tournament_set[rand()%tournament_size];
        return std::make_pair(population[parent1], population[parent2]);
    } else {
        std::vector<Individual> parents = GA::get_top_n(population, 2);
        return std::make_pair(parents[0], parents[1]);
    }
}

std::pair<Individual, Individual> GA::best_cost_route_crossover(std::pair<Individual, Individual> &parents) {

}

void GA::test_fitness() {

}

void GA::intra_depot_mutation(Individual &ind, double prob_rev_mut, double prob_re_routing, double prob_swapping) {
    // Reversal mutation.
    // Single customer re-routing.
    // Swapping (use marginal cost)
}

void GA::inter_depot_mutation(Individual &ind) {}

void GA::mutate(Individual &ind, double prob_rev_mut, double prob_re_routing, double prob_swapping, bool inter_depot_mut) {
    GA::intra_depot_mutation(ind, prob_rev_mut, prob_re_routing, prob_swapping);
    if (inter_depot_mut)
        GA::inter_depot_mutation(ind);
}