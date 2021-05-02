#include "ga.h"

typedef std::pair<int, int> ii;

GA::GA(int pop_size, bool nsga_ii, cv::Mat img): nsga_ii(nsga_ii), image(img) {
    population.clear();
    for (int i=0; i<pop_size; ++i) {
        int treshold = MIN_SEGMENTS+rand()%(MAX_SEGMENTS-MIN_SEGMENTS);
        Individual ind(image, treshold);
        std::cout << i << " " << treshold << std::endl;
        population.push_back(ind);
    }
}

void GA::simulate() {
    int pop_size = population.size();
    std::vector<std::vector<Individual>> sorted_pop;
    for (int gen=0; gen<GENERATIONS; ++gen) {
        std::vector<Individual> children;
        std::cout << "Gen: " << gen << std::endl;
        if (nsga_ii) {
            //std::cout << "sorting" << population.size() << " " << image.rows << " " << image.cols << " " << pop_size<< std::endl;
            if (gen==0) {
                sorted_pop = nsga::sort_and_limit(population, image, pop_size);
                //std::cout << "flatten sorted pop" << std::endl;
                population = flatten(sorted_pop);
            }
            std::cout << "size: " << population.size() << std::endl;
            // fast_nondomination_sort
            while (children.size() < population.size()) {
                //std::cout << "Creating children" << children.size() << std::endl;
                // Binary tournament selection.
                std::pair<ii, ii> parents_pos = binary_tournament_selection(sorted_pop);
                Individual p1 = sorted_pop[parents_pos.first.first][parents_pos.first.second];
                Individual p2 = sorted_pop[parents_pos.second.first][parents_pos.second.second];
                // Recombination / crossover.
                //std::cout << "Crossover" << std::endl;
                std::pair<Individual, Individual> curr_children = crossover(p1, p2);
                // Mutation.
                //std::cout << "Mutate" << std::endl;
                mutate(curr_children.first);
                mutate(curr_children.second);
                children.push_back(curr_children.first);
                children.push_back(curr_children.second);
            }
            population.insert(population.end(), children.begin(), children.end());
            std::cout << "Done with gen " << gen << " population size: " << population.size() << std::endl;
            sorted_pop = nsga::sort_and_limit(population, image, pop_size);
            //std::cout << "flatten sorted pop" << std::endl;
            population = flatten(sorted_pop);
        } else {
            std::vector<double> fitnesses = calculate_fitness(population);
            std::vector<Individual> children;
            // TODO: Implement SGA.
            while (children.size() < population.size()) {
                auto parents = tournament_selection(fitnesses);
                Individual p1 = population[parents.first];
                Individual p2 = population[parents.second];
                std::pair<Individual, Individual> curr_children = crossover(p1, p2);
                mutate(curr_children.first);
                mutate(curr_children.second);
                children.push_back(curr_children.first);
                children.push_back(curr_children.second);
            }
            population = children;
            if (gen==GENERATIONS-1)
                calculate_fitness(population);
        }
    }
}

std::vector<double> GA::calculate_fitness(std::vector<Individual> &inds) {
    std::vector<double> min_values(3, 1e12);
    std::vector<double> max_values(3, 0);
    for (Individual &ind:inds) {
        ind.calculate_objectives(image);
        min_values[0] = std::min(ind.edge_value, min_values[0]);
        min_values[1] = std::min(ind.connectivity, min_values[1]);
        min_values[2] = std::min(ind.overall_deviation, min_values[2]);
        max_values[0] = std::max(ind.edge_value, max_values[0]);
        max_values[1] = std::max(ind.connectivity, max_values[1]);
        max_values[2] = std::max(ind.overall_deviation, max_values[2]);
    }
    std::vector<double> fitness;
    for (Individual &ind:inds) {
        fitness.push_back(
            (1-ind.edge_value/(max_values[0]-min_values[0])) +
            ind.connectivity/(max_values[1]-min_values[1]) +
            ind.overall_deviation/(max_values[2]-min_values[2])
        );
    }
    return fitness;
}

std::pair<int, int> GA::tournament_selection(std::vector<double> inds) {
    std::priority_queue<std::pair<double, int>> selected_inds;
    for (int i=0; i<NUM_IN_TOURNAMENT; ++i) {
        int selected_ind = rand()%inds.size();
        selected_inds.push(std::make_pair(-inds[selected_ind], selected_ind));
    }
    int p1 = selected_inds.top().second;
    selected_inds.pop();
    int p2 = selected_inds.top().second;
    return std::make_pair(p1, p2);
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
    //int crossover_pos = rand()%(image.rows * image.cols);
    // kruskal_prim_mst
    UnionFind uf(image.rows*image.cols);
    auto seg = uf.kruskal_prim_mst(2, image).first;
    int seg_1 = seg[0];
    for (int pos=0; pos<seg.size(); ++pos) {
        if (seg[pos]==seg_1)
            seg[pos] = 0;
        else
            seg[pos] = 1;
    }
    //std::cout << "crossover_pos " << crossover_pos << std::endl;
    //std::cout << "p1.genes.size() " << p1.genes.size() << " p2.genes.size() " << p2.genes.size() << std::endl;
    //std::cout << "p1.root.size() " << p1.root.size() << " p2.root.size() " << p2.root.size() << std::endl;
    //Individual c1(p1, p2, crossover_pos, image);
    Individual c1(p1, p2, seg, image.rows, image.cols);
    std::cout << "created c1" << std::endl;
    //Individual c2(p2, p1, crossover_pos, image);
    Individual c2(p2, p1, seg, image.rows, image.cols);
    std::cout << "created c2" << std::endl;
    return std::make_pair(c1, c2);
}


void GA::mutate(Individual &ind) {
    // Split: Find the largest edge and make it point to itself (stochastically).
    if ((double)(rand())/(RAND_MAX)<SPLIT_MUTATION) {
        int max_pos = 0;
        double max_dist = 0;
        for (int pos=0; pos<ind.root.size(); ++pos) {
            Dir actual_dir = get_actual_dir(ind.genes[pos], pos, image.cols, image.rows);
            if (actual_dir==Dir::s)
                continue;
            int parent = find_pos(pos, actual_dir, image.cols, image.rows);
            int cur_x = pos % image.cols;
            int cur_y = pos/ image.cols;
            int parent_x = parent % image.cols;
            int parent_y = parent / image.cols;
            double cur_dist = euc_dist(image.at<cv::Vec3b>(cur_y, cur_x), image.at<cv::Vec3b>(parent_y, parent_x));
            if (cur_dist>max_dist) {
                max_pos = pos;
                max_dist = cur_dist;
            }
        }
        ind.genes[max_pos] = Dir::s;
        ind.find_roots();
    }
    // Connect: Find two neighbouring groups, prim on relevant edges.
    /*
    // find two pixels from different groups.
    // for all pixels 
    //  if pixel in one of the groups
    //      for dir in directions (see main.cpp)
    //          ...
    //          add_edge(...)
    //          ...
    // prim_mst
    // get root of cluster (from one of the pxls)
    // for all pixels
    //  if pixel in one of the groups
    //      change dir to prim_mst
    // find_roots()
    */
    for (int pos=0; pos<ind.genes.size(); ++pos) {
        if ((double)(rand())/(RAND_MAX)<RANDOM_MUTATION) {
            ind.mutate(pos);
            std::cout << "mutation" << std::endl;
        }
    }
    ind.find_roots();
}