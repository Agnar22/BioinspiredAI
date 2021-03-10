#include "problem.h"

Problem::Problem(
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

void Problem::calculate_distances() {
    std::cout << "positions size " << positions.size() << std::endl;
    distances.resize(positions.size(), std::vector<double>(positions.size()));
    for (int from=0; from<positions.size(); ++from) {
        for (int to=0; to<from; ++to) {
            double dist = euclidian_dist(positions[from], positions[to]);
            distances[from][to]=dist;
            distances[to][from]=dist;
        }
    }
}

void Problem::find_depot_distances() {
    for (int customer=0; customer<num_customers; ++customer) {
        depot_distances.push_back(std::vector<double>(distances[customer].begin()+num_customers, distances[customer].end()));
        std::cout << "Depot distance: " << depot_distances[customer].size() << " " << distances[customer].size() << std::endl; 
        std::cout << num_customers << std::endl;
    }
}

double Problem::euclidian_dist(std::pair<double, double> p1, std::pair<double, double> p2) {
    return sqrt(pow(p1.first-p2.first, 2)+pow(p1.second-p2.second, 2));
}

double Problem::get_distance(int from, int to) {return distances[from][to];}

std::vector<std::vector<double>> Problem::get_distances(){return distances;}

std::vector<std::vector<double>> Problem::get_depot_distances(){return depot_distances;}

int Problem::get_num_depots(){return num_depots;}

int Problem::get_num_customers(){return num_customers;}

int Problem::get_vhcl_pr_depot(){return vhcl_pr_depot;}

double Problem::get_customer_load(int cust){return cust_demand[cust];}

std::vector<double> Problem::get_max_lengths(){return max_length;}

double Problem::get_max_length(int depot){return max_length[depot];}

double Problem::get_max_load(int depot){return max_load[depot];}

