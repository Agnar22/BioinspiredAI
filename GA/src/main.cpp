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

#define assertm(exp, msg) assert(((void)msg, exp))
#define PI 3.1415926535


constexpr int BITSTRING_SIZE = 6;

namespace fitness {
    double sin(std::bitset<BITSTRING_SIZE> val){
        // TODO: handle exception when BITSTRING_SIZE > 64.
        assertm(BITSTRING_SIZE<=sizeof(unsigned long long), "The bitstring can be cast to a uint_64.");
        double sin_max = 128;
        unsigned long long bitstring_max = std::bitset<BITSTRING_SIZE>().set().to_ullong();
        return std::sin((double) val.to_ullong()/(double) bitstring_max * sin_max * PI / 180.0) + 1;
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
        genes = std::bitset<BITSTRING_SIZE>();
    }

    Individual (std::string gene_sequence) {
        genes = std::bitset<BITSTRING_SIZE>(std::string(gene_sequence));
    }
};

class SimpleGA {
    public:

        SimpleGA(int population_size, double (*fitness_function)(std::bitset<BITSTRING_SIZE>)) {
        }

        void simulate() {
        }

        int get_population_size(){
        }

    private:
        double (*fitness_function)(std::bitset<BITSTRING_SIZE>);
        std::vector<Individual> population;

        void test_fitness() {
        }

        std::vector<std::pair<int, int>> select() {
        }

        double get_random_fitness(double max_fitness) {
        }

        int get_pos_from_fitness(double fitness) {
        }

        std::vector<Individual> crossover(std::vector<std::pair<int, int>> selected_pairs) {
        }

        std::vector<Individual> mutate(std::vector<Individual> children) {
        }
};

int main() {
    auto data = read_csv("data.csv");
    std::cout << "read file " << data.size() << " " << data[0].size() << std::endl;

    SimpleGA SGA(10, &fitness::sin);
    std::cout << "population size " << SGA.get_population_size() << std::endl;
    std::cout << fitness::sin(std::bitset<BITSTRING_SIZE>().set()) << std::endl;
    std::cout << fitness::sin(std::bitset<BITSTRING_SIZE>()) << std::endl;
    std::cout << fitness::sin(std::bitset<BITSTRING_SIZE>("100000")) << std::endl;

}
