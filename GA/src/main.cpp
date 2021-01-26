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
#include <thread>

#define assertm(exp, msg) assert(((void)msg, exp))
#define PI 3.1415926535


constexpr int BITSTRING_SIZE = 10;
constexpr int POPULATION_SIZE = 100;
constexpr int EPOCHS = 100000;
constexpr double CROSSOVER_PROB = 0.4;
constexpr double MUTATION_PROB = 0.06;
constexpr double KILL_LOWEST = 0.3;


namespace fitness {
    double sin(std::bitset<BITSTRING_SIZE> val){
        // TODO: handle exception when BITSTRING_SIZE > 64.
        assertm(BITSTRING_SIZE<=sizeof(unsigned long long)*8, "The bitstring can be cast to a uint_64.");
        double sin_max = 128;
        unsigned long long bitstring_max = std::bitset<BITSTRING_SIZE>().set().to_ullong();
        //std::cout << bitstring_max << " ";
        //std::cout << val.to_ullong() << " ";
        //std::cout << (double)val.to_ullong() << " ";
        //std::cout << (double)bitstring_max << " ";
        //std::cout << (double) val.to_ullong()/(double) bitstring_max * sin_max << " ";
        return std::sin(((double) val.to_ullong()/(double) bitstring_max) * sin_max) + 1;
    }

    double linear(std::bitset<BITSTRING_SIZE> val){
        // TODO: handle exception when BITSTRING_SIZE > 64.
        assertm(BITSTRING_SIZE<=sizeof(unsigned long long)*8, "The bitstring can be cast to a uint_64.");
        unsigned long long bitstring_max = std::bitset<BITSTRING_SIZE>().set().to_ullong();
        return (double) val.to_ullong()/(double) bitstring_max;
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
        genes = std::bitset<BITSTRING_SIZE>(rand()%(int)std::pow(2, BITSTRING_SIZE));
        //std::cout << "genes " << genes << std::endl;
    }

    Individual (std::string gene_sequence) {
        genes = std::bitset<BITSTRING_SIZE>(std::string(gene_sequence));
        //std::cout << "genes " << gene_sequence << " " <<  genes << std::endl;
    }

    bool operator < (const Individual& other) {
        return fitness < other.fitness;
    }
};

class SimpleGA {
    public:

        SimpleGA(int population_size, double (*fitness_function)(std::bitset<BITSTRING_SIZE>)) {
            this->fitness_function = fitness_function;
            this->population_size = population_size;
            for (int x=0; x<population_size; ++x)
                population.push_back(Individual());
            //std::cout << "created sga." << std::endl;
        }

        void simulate() {
            test_fitness();
            for (int x=0; x<EPOCHS; ++x) {
                auto selected_pairs = select();
                //std::cout << "selected_pairs " << selected_pairs.size() << std::endl;
                auto children = crossover(selected_pairs);
                //std::cout << "children " << children.size() << std::endl;
                population = mutate(children);
                //std::cout << "population " << population.size() << std::endl;
                kill_weakest();
                test_fitness();
            }
        }

        void print_population() {
            for (int y=0; y<population.size(); ++y)
                std::cout << population[y].genes << " " << fitness_function(population[y].genes) << std::endl;
            std::cout << std::endl;
        }

        int get_population_size(){
            return population.size();
        }

        void store_population(std::string file_name) {
            std::ofstream file(file_name);
            for (Individual &indiv:population)
                file << indiv.genes.to_string() << std::endl;
            file.close();
        }

    private:

        double (*fitness_function)(std::bitset<BITSTRING_SIZE>);
        int population_size;
        std::vector<Individual> population;

        void test_fitness() {
            double tot_fitness = 0.0;
            double max_fitness = 0.0;
            for (auto &individual:population) {
                individual.fitness = fitness_function(individual.genes);
                tot_fitness += individual.fitness;
                max_fitness = std::max(individual.fitness, max_fitness);
            }
            std::cout << "Max fitness: " << max_fitness << " avg. fitness: " << tot_fitness/population.size()<< std::endl;
        }

        std::vector<std::pair<int, int>> select() {
            //double tot_fitness = std::accumulate(population.begin(), population.end(), 0.0);
            double tot_fitness = 0;
            for (auto &individual:population)
                tot_fitness += individual.fitness;
            std::vector<std::pair<int, int>> selected_pairs;

            //std::cout << "population_size " << population_size << std::endl;

            for (int x=0; x<population_size/2; ++x) {
                int parent_one = get_pos_from_fitness(get_random_number(tot_fitness));
                int parent_two = get_pos_from_fitness(get_random_number(tot_fitness));
                while (parent_one == parent_two) {
                    parent_one = get_pos_from_fitness(get_random_number(tot_fitness));
                    parent_two = get_pos_from_fitness(get_random_number(tot_fitness));
                }
                //std::cout << "parent " << parent_one << " " << parent_two << std::endl;
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
            for (std::pair<int, int> &parents:selected_pairs) {
                if (get_random_number(1.0) > CROSSOVER_PROB) {
                    children.push_back(population[parents.first]);
                    children.push_back(population[parents.second]);
                } else {
                    int idx = rand() % BITSTRING_SIZE;
                    for (int x=0; x<2; ++x) {
                        //std::cout << "parl " << parents.first << " parr " << parents.second << std::endl;
                        std::string parent_one_genes = population[parents.first].genes.to_string().substr(0,idx);
                        std::string parent_two_genes = population[parents.second].genes.to_string().substr(idx);
                        Individual child(parent_one_genes + parent_two_genes);
                        //std::cout << "parl " << parent_one_genes << " parr " << parent_two_genes << " child " << child.genes << " fitness " << fitness_function(child.genes) <<  std::endl;
                        children.push_back(child);

                        int first = parents.first;
                        parents.first = parents.second;
                        parents.second = first;
                    }
                }
            }
            return children;
        }

        std::vector<Individual> mutate(std::vector<Individual> children) {
            for (Individual &indiv:children) {
                for (int pos=0; pos<BITSTRING_SIZE; ++pos) {
                    if (get_random_number(1.0) < MUTATION_PROB) {
                        //std::cout << pos << "flipped " << indiv.genes;
                        //std::cout << " " << fitness_function(indiv.genes);
                        indiv.genes.flip(pos);
                        //std::cout << " " << indiv.genes;
                        //std::cout << " " << fitness_function(indiv.genes) <<  std::endl;
                    }
                }
            }
            return children;
        }

        void kill_weakest() {
            // TODO: kill the n weakest individuals.

            for (auto &individual:population)
                individual.fitness = fitness_function(individual.genes);
            std::sort(population.begin(), population.end());
            //for (Individual &child:population)
            //  std::cout << child.genes.to_string() << " " << child.fitness << "  ";
            std::cout << std::endl;
            population.erase(population.begin(), population.begin()+(int)(population.size()*KILL_LOWEST));
            std::cout << "size " << population.size() << std::endl;
        }
};

void start_python_fitness() {
    system("./run_python.sh");
}

int main() {
    auto data = read_csv("data.csv");
    std::cout << "read file " << data.size() << " " << data[0].size() << std::endl;


    //SimpleGA SGA(100, &fitness::sin);
    //SimpleGA SGA(POPULATION_SIZE, &fitness::sin);
    //SGA.print_population();
    //std::cout << "population size " << SGA.get_population_size() << std::endl;
    //std::cout << fitness::sin(std::bitset<BITSTRING_SIZE>().set()) << std::endl;
    //std::cout << fitness::sin(std::bitset<BITSTRING_SIZE>()) << std::endl;
    //std::cout << fitness::sin(std::bitset<BITSTRING_SIZE>("100000")) << std::endl;

    //SGA.simulate();
    //SGA.print_population();
    //SGA.store_population("sin.txt");
    std::thread python_fitess(start_python_fitness);
    //system("./run_python.sh");
    std::cout << "test" << std::endl;
    std::chrono::milliseconds timespan(111605);
    std::this_thread::sleep_for(timespan);
    std::cout << "test" << std::endl;
}
