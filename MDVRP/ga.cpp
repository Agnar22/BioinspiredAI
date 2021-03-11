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
            std::pair<Individual, Individual> children = GA::best_cost_route_crossover(parents, pr);

            // Mutation.
            GA::mutate(children.first, prob_rev_mut, prob_re_routing, prob_swapping, gen%inter_depot_swapping==0, pr);
            GA::mutate(children.second, prob_rev_mut, prob_re_routing, prob_swapping, gen%inter_depot_swapping==0, pr);

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

std::pair<Individual, Individual> GA::best_cost_route_crossover(std::pair<Individual, Individual> &parents, Problem &pr) {
    int num_depots = parents.first.chromosome_trips.size();
    int depot = rand()%num_depots;
    int pf_route = rand()%parents.first.chromosome_trips[depot].size();
    int ps_route = rand()%parents.second.chromosome_trips[depot].size();

    Individual f_child = parents.first;
    Individual s_child = parents.second;

    // Remove customers in ps_route from pf_route and vice versa.
    std::vector<int> rmd_cust_f = parents.second.chromosome_trips[depot][ps_route];
    std::vector<int> rmd_cust_s = parents.first.chromosome_trips[depot][pf_route];
    f_child.remove_customers(rmd_cust_f, pr);
    s_child.remove_customers(rmd_cust_s, pr);

    // Stochastically add customers, that were previously removed, to depot.
    std::for_each(rmd_cust_f.begin(), rmd_cust_f.end(), [&] (int cust) {f_child.insert_stochastically(cust, 0.8, depot, pr);});
    std::for_each(rmd_cust_s.begin(), rmd_cust_s.end(), [&] (int cust) {s_child.insert_stochastically(cust, 0.8, depot, pr);});

    return std::make_pair(f_child, s_child);
}

void GA::test_fitness() {

}

void GA::mutate(Individual &ind, double prob_rev_mut, double prob_re_routing, double prob_swapping, bool inter_depot_mut, Problem &pr) {
    if (!inter_depot_mut)
        GA::intra_depot_mutation(ind, prob_rev_mut, prob_re_routing, prob_swapping, pr);
    else
        GA::inter_depot_mutation(ind);
}

void GA::intra_depot_mutation(Individual &ind, double prob_rev_mut, double prob_re_routing, double prob_swapping, Problem &pr) {
    // Reversal mutation.
    if (prob_rev_mut<(double)(rand())/(double)(RAND_MAX));
        ind.reversal_mutation(rand()%pr.get_num_depots(), pr);
    // Single customer re-routing.
    if (prob_re_routing<(double)(rand())/(double)(RAND_MAX));
        ind.re_routing_mutation(rand()%pr.get_num_depots(), pr);
    // Swapping (use marginal cost)
    if (prob_swapping<(double)(rand())/(double)(RAND_MAX));
        ind.swapping_mutation(rand()%pr.get_num_depots(), pr);
}

void GA::inter_depot_mutation(Individual &ind) {}
