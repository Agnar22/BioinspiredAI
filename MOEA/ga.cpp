#include "ga.h"

typedef std::pair<int, int> ii;

GA::GA(int pop_size, bool nsga_ii, cv::Mat img): nsga_ii(nsga_ii), image(img) {
    population.clear();
    for (int i=0; i<pop_size; ++i) {
        Individual ind(image);
        population.push_back(ind);
    }
}

void GA::simulate() {
    int pop_size = population.size();
    for (int gen=0; gen<100; ++gen) {
        std::vector<Individual> children;
        std::cout << "Gen: " << gen << std::endl;
        if (nsga_ii) {
            std::cout << "sorting" << population.size() << " " << image.rows << " " << image.cols << " " << pop_size<< std::endl;
            auto sorted_pop = nsga::sort_and_limit(population, image, pop_size);
            std::cout << "flatten sorted pop" << std::endl;
            population = flatten(sorted_pop);
            std::cout << "size: " << population.size() << std::endl;
            // fast_nondomination_sort
            while (children.size() < population.size()) {
                std::cout << "Creating children" << children.size() << std::endl;
                // Binary tournament selection.
                std::pair<ii, ii> parents_pos = binary_tournament_selection(sorted_pop);
                std::cout << "Binary tournament selection" << std::endl;
                Individual p1 = sorted_pop[parents_pos.first.first][parents_pos.first.second];
                Individual p2 = sorted_pop[parents_pos.second.first][parents_pos.second.second];
                // Recombination / crossover.
                std::cout << "Crossover" << std::endl;
                std::pair<Individual, Individual> curr_children = crossover(p1, p2);
                // Mutation.
                std::cout << "Mutate" << std::endl;
                mutate(curr_children.first);
                mutate(curr_children.second);
                std::cout << "Push back" << std::endl;
                children.push_back(curr_children.first);
                children.push_back(curr_children.second);
                std::cout << "Done with curr_children" << std::endl;
            }
            std::cout << "Insert into population" << std::endl;
            population.insert(population.end(), children.begin(), children.end());
            std::cout << "Done with gen " << gen << " population size: " << population.size() << std::endl;
        } else {
            // TODO: Implement SGA.
        }
    }
}

std::pair<ii, ii> GA::binary_tournament_selection(std::vector<std::vector<Individual>> &sorted_parents) {
    return std::make_pair(select_parent_pos(sorted_parents), select_parent_pos(sorted_parents));
}

ii GA::select_parent_pos(std::vector<std::vector<Individual>> &sorted_parents) {
    ii p1 = std::make_pair(rand()%sorted_parents.size(), -1);
    p1.second = rand()%sorted_parents[p1.first].size();
    ii p2 = std::make_pair(rand()%sorted_parents.size(), -1);
    p2.second = rand()%sorted_parents[p2.first].size();
    return p1<p2? p1:p2;
}

std::pair<Individual, Individual> GA::crossover(Individual &p1, Individual &p2) {
    // Single point crossover.
    int crossover_pos = rand()%(image.rows * image.cols);
    std::cout << "crossover_pos " << crossover_pos << std::endl;
    std::cout << "p1.genes.size() " << p1.genes.size() << " p2.genes.size() " << p2.genes.size() << std::endl;
    std::cout << "p1.root.size() " << p1.root.size() << " p2.root.size() " << p2.root.size() << std::endl;
    Individual c1(p1, p2, crossover_pos, image);
    std::cout << "created c1" << std::endl;
    Individual c2(p1, p2, crossover_pos, image);
    std::cout << "created c2" << std::endl;
    return std::make_pair(c1, c2);
}


void GA::mutate(Individual &ind) {
    for (int pos=0; pos<ind.genes.size(); ++pos)
        if ((double)(rand())/(RAND_MAX)<0.1)
            ind.mutate(pos);
}