#include <algorithm>
#include <numeric>
#include <random>
#include <functional>
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


class Problem {
    private:
        int num_customers, num_depots, vhcl_pr_depot;
        std::vector<double> max_length, max_load;
        std::vector<std::pair<double, double>> positions;
        std::vector<double> cust_serv_dur, cust_demand;
        std::vector<std::vector<double>> distances;
        std::vector<std::vector<double>> depot_distances;

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

        void find_depot_distances() {
            for (int customer=0; customer<num_customers; ++customer) {
                depot_distances.push_back(std::vector<double>(distances[customer].begin()+num_customers, distances[customer].end()));
        }

        double euclidian_dist(std::pair<double, double> p1, std::pair<double, double> p2) {
            return sqrt(pow(p1.first-p2.first, 2)+pow(p1.second-p2.second, 2));
        }

    public:

    Problem(
        int num_customers,
        int num_depots,
        int vhcl_pr_depot,
        std::vector<double> max_length,
        std::vector<double> max_load,
        std::vector<std::pair<double, double>> positions,
        std::vector<double> cust_serv_dur,
        std::vector<double> cust_demand
    ): num_customers{num_customers},
       num_depots{num_depots},
       vhcl_pr_depot{vhcl_pr_depot},
       max_length{max_length},
       max_load{max_load},
       positions{positions},
       cust_serv_dur{cust_serv_dur},
       cust_demand{cust_demand} {

        std::cout << "Calculating distances." << std::endl;
        calculate_distances();
        std::cout << "Calculating depot distances." << std::endl;
        find_depot_distances();
        std::cout << "Finished constructing problem." << std::endl;
    }

    double get_distance(int from, int to) {return distances[from][to];}

    std::vector<std::vector<double>> get_distances(){return distances;}

    std::vector<std::vector<double>> get_depot_distances(){return depot_distances;}

    int get_num_depots(){return num_depots;}

    int get_num_customers(){return num_customers;}

    int get_vhcl_pr_depot(){return vhcl_pr_depot;}

    std::vector<double> get_max_length(){return max_length;}
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

        int max_vhcl_pr_depot = (int)prob[0][0];
        int num_cust = (int)prob[0][1];
        int num_depots = (int)prob[0][2];

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
        for (int row = num_cust+num_depots+1; row < num_cust+2*num_depots+1; ++row) {
            positions.push_back(std::make_pair(prob[row][1], prob[row][2]));
        }


        std::cout << "Finished cleaning data." << std::endl;

        return Problem(num_cust, num_depots, max_vhcl_pr_depot, max_length, max_load, positions, cust_serv_dur, cust_demand);
    }
}


struct Individual {
    std::vector<std::vector<int>> chromosomes;
    std::vector<std::vector<int>> cust_on_depots;
    std::vector<std::vector<int>> chromosome_trip_ends;
    double fitness;

    Individual (Problem pr) {
        // TODO: Initialize individual.
        /**
        * Ways to initialize:
        *  - Deterministically/stochastically assign points to depots based on distance.
        *    - Create routes based on closest point deterministically/stochastically.
        **/

        initialize_chromosomes(pr);
    }

    void initialize_chromosomes(Problem &pr) {
        chromosomes.resize(pr.get_num_depots());
        // Assign customers to depots.
        cust_on_depots = assign_customers_to_depots(pr);
        // Assign routes to depots.
        for (int depot=0; depot<pr.get_num_depots(); ++depot) {
            chromosomes[depot] = order_depot(depot, cust_on_depots[depot], pr);
        }
    }

    std::vector<std::vector<int>> assign_customers_to_depots(Problem &pr) {
        // Stochastically assign customers to depots based on inverse distanse to depot.
        std::vector<std::vector<int>> customers_on_depot(pr.get_num_depots());
        for (int cust=0; cust<pr.get_num_customers(); ++cust) {
            std::vector<std::vector<double>> depot_dist = pr.get_depot_distances();
            std::vector<double> depot_prob;
            for (double dist:depot_dist[cust]) {
                depot_prob.push_back(1/std::pow(dist, 3));
            }
            double sum=0;
            std::for_each(depot_prob.begin(), depot_prob.end(), [&] (double n) {sum += n;});
            for (int depot=0; depot<pr.get_num_depots(); ++depot) {
                depot_prob[depot]=depot_prob[depot]/sum;
            }
            double prob = (double)(rand()) / (double)(RAND_MAX);
            double cum_sum = 0;
            for (int pos=0; pos<depot_prob.size(); ++pos) {
                cum_sum+=depot_prob[pos];
                if (cum_sum>=prob) {
                    customers_on_depot[pos].push_back(cust);
                    break;
                }
            }
        }
        return customers_on_depot;
    }

    template <class T>
    static std::vector<T> get_subset(std::vector<T> &vals, std::vector<int> idxs) {
        std::vector<T> subset(idxs.size());
        for (int i=0; i<idxs.size(); ++i)
            subset[i]=vals[i];
        return subset;
    }

    std::vector<int> order_depot(int depot_num, std::vector<int> customers, Problem &pr) {
        std::vector<int> order;
        double trip_dist = 0;
        int cur_cust = pr.get_num_customers() + depot_num;
        int num_vhcl = 1;
        std::vector<std::vector<double>> distances = pr.get_distances();
        while (customers.size() > 0) {
            std::vector<double> cust_distances = Individual::get_subset<double>(pr.get_distances()[cur_cust], customers);
            int closest_cust_pos = std::min_element(cust_distances.begin(), cust_distances.end()) - cust_distances.begin();
            int closest_cust = customers[closest_cust_pos];
            // Check that the current trip is not too long.
            if ((pr.get_max_length()[depot_num]==0) || trip_dist+distances[cur_cust][closest_cust]+distances[closest_cust][pr.get_num_customers()+depot_num] <= pr.get_max_length()[depot_num]){
                order.push_back(closest_cust);
                trip_dist += distances[cur_cust][closest_cust];
                customers.erase(customers.begin()+closest_cust_pos);
            } else if (trip_dist==0) {
                std::cout << "Trip could not be created as first point is too far away: " << closest_cust << " with distance "  << distances[closest_cust][pr.get_num_customers()+depot_num] << std::endl;
            } else {
                cur_cust = pr.get_num_customers() + depot_num;
                trip_dist = 0;
                ++num_vhcl;
            }
        }
        if (num_vhcl>pr.get_vhcl_pr_depot()) {
            std::cout << "Too many vehicles for depot: " << depot_num << " with "  << num_vhcl << " vehicles." << std::endl;
        }
        return order;
    }

    void find_direct_mapping() {
        chromosome_trip_ends.resize(chromosomes.resize());
        for (int num_chrm; num_chrm<chromosomes.size(); ++num_chrm) {

        }
    }
};

class GA {

};


int main() {
    std::cout << "Hello World!" << std::endl;
    std::string file_name = "Data files project 2/Testing Data/Data Files/p01";
    std::vector<std::vector<double>> problem_file = file::read_flat(file_name);
    std::cout << file_name << std::endl;

    Problem curr_prob = file::load_problem(file_name);
    std::cout.precision(17);
    std::cout << "distance " <<  curr_prob.get_distance(1,2) << std::endl;
    std::cout << "distance " <<  curr_prob.get_distance(3,2) << std::endl;
    Individual ind(curr_prob);
}
