#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <regex>
#include <utility>
#include <math.h>
#include <limits>

// TODO: representation of individuals.
// Ideas:
//  - initialize customer to closest depot (obs. be carefull of edgecase where it is not illegal to have all at for one depot or similar)
//  - elitism
//  - variation operators with different kinds of "jumps" - small perfections and larger jumps
//  - niching


struct Individual {

};

class GA {

};

class Problem {
    private:
        int num_customers, num_depots;
        std::vector<double> max_length, max_load;
        std::vector<std::pair<double, double>> positions;
        std::vector<double> cust_serv_dur, cust_demand;
        std::vector<std::vector<double>> distances;

        void calculate_distances() {
            distances.resize(positions.size(), std::vector<double>(positions.size()));
            for (int from=0; from<positions.size(); ++from) {
                for (int to=0; to<from; ++to) {
                    std::cout << "from " << from << " to " << to << std::endl;
                    double dist = euclidian_dist(positions[from], positions[to]);
                    distances[from][to]=dist;
                    distances[to][from]=dist;
                }
            }
        }

        double euclidian_dist(std::pair<double, double> p1, std::pair<double, double> p2) {
            return sqrt(pow(p1.first-p2.first, 2)+pow(p1.second-p2.second, 2));
        }

    public:

    Problem(
        std::vector<double> max_length,
        std::vector<double> max_load,
        std::vector<std::pair<double, double>> positions,
        std::vector<double> cust_serv_dur,
        std::vector<double> cust_demand
    ) {
        this->max_length = max_length;
        this->max_load = max_load;
        this->positions = positions;
        this->cust_serv_dur = cust_serv_dur;

        std::cout << "Calculating distance." << std::endl;
        calculate_distances();
        std::cout << "Calculated distance." << std::endl;
    }

    double get_distance(int from, int to) {
        return distances[from][to];
    }
};

namespace file {
    std::vector<std::vector<double>> read_flat(std::string file_name) {
        std::ifstream fs(file_name);
        std::vector<std::vector<double>> file;
        std::string line;
        std::cout << "Opening file " << file_name << std::endl;
        while (std::getline(fs, line)) {
            // Remove leading and trailing spaces.
            line = std::regex_replace(line, std::regex("^ +| +$|( ) +"), "$1");
            std::istringstream s(line);
            std::string val;
            std::vector<double> row;
            while (std::getline(s, val, ' '))
                if (val.compare(" ")!=0)
                    row.push_back(std::stod(val));
            file.push_back(row);
        }
        fs.close();
        return file;
    }

    Problem load_problem(std::string file_name) {
        std::vector<std::vector<double>> prob = read_flat(file_name);

        int max_vhcl_pr_depot = prob[0][0];
        int num_cust = prob[0][1];
        int num_depots = prob[0][2];

        std::vector<double> max_length;
        std::vector<double> max_load;
        for (int row = 1; row < num_depots+1; ++row) {
            max_length.push_back(prob[row][0]);
            max_load.push_back(prob[row][1]);
        }

        std::vector<std::pair<double, double>> positions;
        std::vector<double> cust_serv_dur, cust_demand;
        for (int row = num_depots+1; row < num_cust+num_depots+1; ++row) {
            positions.push_back(std::make_pair(prob[row][1], prob[row][2]));
            cust_serv_dur.push_back(prob[row][3]);
            cust_demand.push_back(prob[row][4]);
        }

        std::cout << "Finished cleaning data." << std::endl;

        return Problem(max_length, max_load, positions, cust_serv_dur, cust_demand);
    }
}



int main() {
    std::cout << "Hello World!" << std::endl;
    std::string file_name = "Data files project 2/Testing Data/Data Files/p01";
    std::vector<std::vector<double>> problem_file = file::read_flat(file_name);
    std::cout << file_name << std::endl;
    for (auto line:problem_file) {
        for (auto val:line)
            std::cout << val << " ";
        std::cout << std::endl;
    }

    Problem curr_prob = file::load_problem(file_name);
    std::cout.precision(17);
    std::cout << "distance " <<  curr_prob.get_distance(1,2) << std::endl;
    std::cout << "distance " <<  curr_prob.get_distance(3,2) << std::endl;
}
