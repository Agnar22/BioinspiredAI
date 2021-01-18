#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include <fstream>
#include <bitset>
#include <math.h>
#include <cassert>
#include <numeric>
#include <utility>
#include <stdlib.h>
#include <stdexcept>
#include <algorithm>

#define assertm(exp, msg) assert(((void)msg, exp))
#define PI 3.1415926535


constexpr int BITSTRING_SIZE = 15;
constexpr double MUTATION_PROB = 0.02;

namespace fitness {
    double sin(std::bitset<BITSTRING_SIZE> val){
        // TODO: handle exception when BITSTRING_SIZE > 64.
        assertm(BITSTRING_SIZE<=sizeof(unsigned long long)*8, "The bitstring can be cast to a uint_64.");
        double sin_max = 128;
        unsigned long long bitstring_max = std::bitset<BITSTRING_SIZE>().set().to_ullong();
        return std::sin(((double) val.to_ullong()/(double) bitstring_max) * sin_max) + 1;
    }
}

std::vector<std::vector<double>> read_csv(std::string file_name) {
    std::ifstream fs(file_name);
    std::vector<std::vector<double>> file;
    std::string line;
    while (std::getline(fs, line)) {
        std::istringstream s(line);
        std::string val;
        std::vector<double> row;
        while (std::getline(s, val, ','))
            row.push_back(std::stod(val));
        file.push_back(row);
    }
    fs.close();
    return file;
}

struct Individual {
    std::bitset<BITSTRING_SIZE> genes;
    double fitness;

    Individual () {
        // TODO: random initialization.
        genes = std::bitset<BITSTRING_SIZE>(rand()%(int)std::pow(2, BITSTRING_SIZE));
    }

    Individual (std::string gene_sequence) {
        genes = std::bitset<BITSTRING_SIZE>(std::string(gene_sequence));
    }
};

class SimpleGA {
    public:

        SimpleGA(int population_size, double (*fitness_function)(std::bitset<BITSTRING_SIZE>)) {
            this->fitness_function = fitness_function;
            for (int x=0; x<population_size; ++x)
                population.push_back(Individual());
        }

        void simulate() {
            int x = 300;
            while (x) {
                test_fitness();
                auto selected_pairs = select();
                auto children = crossover(selected_pairs);
                population = mutate(children);
                x--;
            }
        }

        void print_population() {
            for (int y=0; y<population.size(); ++y)
                std::cout << population[y].genes << " ";
            std::cout << std::endl;
        }

        int get_population_size(){
            return population.size();
        }

    private:
        double (*fitness_function)(std::bitset<BITSTRING_SIZE>);
        std::vector<Individual> population;

        void test_fitness() {
            double max_fitness = 0.0;
            for (auto individual:population) {
                individual.fitness = fitness_function(individual.genes);
                max_fitness = std::max(individual.fitness, max_fitness);
            }
        }

        std::vector<std::pair<int, int>> select() {
            double tot_fitness = 0;
            for (auto individual:population)
                tot_fitness += individual.fitness;
            int population_size = population.size();
            std::vector<std::pair<int, int>> selected_pairs;

            for (int x=0; x<population_size/2; ++x) {
                int parent_one = get_pos_from_fitness(get_random_number(tot_fitness));
                int parent_two = get_pos_from_fitness(get_random_number(tot_fitness));
                while (parent_one == parent_two) {
                    parent_one = get_pos_from_fitness(get_random_number(tot_fitness));
                    parent_two = get_pos_from_fitness(get_random_number(tot_fitness));
                }
                selected_pairs.push_back(std::make_pair(parent_one, parent_two));
            }
            return selected_pairs;
        }

        double get_random_number(double scaling) {
            return ((double) rand() / (double)(RAND_MAX)) * scaling;
        }

        int get_pos_from_fitness(double fitness) {
            double cumulative_sum = 0.0;
            for (int x=0; x<population.size(); ++x) {
                cumulative_sum += population[x].fitness;
                if (cumulative_sum >= fitness)
                    return x;
            }
            throw std::invalid_argument("Fitness exceeded cumulative fitness.");
        }

        std::vector<Individual> crossover(std::vector<std::pair<int, int>> selected_pairs) {
            std::vector<Individual> children;
            for (std::pair<int, int> parents:selected_pairs) {
                int idx = rand() % BITSTRING_SIZE;

                for (int x=0; x<2; ++x) {
                    std::string parent_one_genes = population[parents.first].genes.to_string().substr(0,idx);
                    std::string parent_two_genes = population[parents.second].genes.to_string().substr(idx);
                    Individual child(parent_one_genes + parent_two_genes);
                    children.push_back(child);

                    int first = parents.first;
                    parents.first = parents.second;
                    parents.second = first;
                }
            }
            return children;
        }

        std::vector<Individual> mutate(std::vector<Individual> children) {
            for (Individual &indiv:children) {
                for (int pos=0; pos<BITSTRING_SIZE; ++pos) {
                    if (get_random_number(1.0) < MUTATION_PROB) {
                        indiv.genes.flip(pos);
                    }
                }
            }
            return children;
        }
};

int main() {
    auto data = read_csv("data.csv");
    std::cout << "read file " << data.size() << " " << data[0].size() << std::endl;

    SimpleGA SGA(100, &fitness::sin);
    std::cout << "population size " << SGA.get_population_size() << std::endl;
    std::cout << fitness::sin(std::bitset<BITSTRING_SIZE>().set()) << std::endl;
    std::cout << fitness::sin(std::bitset<BITSTRING_SIZE>()) << std::endl;
    std::cout << fitness::sin(std::bitset<BITSTRING_SIZE>("100000")) << std::endl;

    SGA.simulate();
}
