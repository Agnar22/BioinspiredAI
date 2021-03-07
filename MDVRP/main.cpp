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

    double get_customer_load(int cust){return cust_demand[cust];}

    std::vector<double> get_max_length(){return max_length;}
};

struct Individual {
    std::vector<std::vector<int>> chromosomes;
    std::vector<std::vector<int>> cust_on_depots;
    std::vector<std::vector<std::vector<int>>> chromosome_trips;
    std::vector<std::vector<double>> trip_dists;
    std::vector<std::vector<int>> trip_loads;
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
        chromosome_trips.resize(chromosomes.size());
        trip_dists.resize(chromosomes.size());
        trip_loads.resize(chromosomes.size());
        for (int depot=0; depot<pr.get_num_depots(); ++depot) {
            setup_trips(depot, cust_on_depots[depot], pr);
        }
    }

    std::vector<std::vector<int>> assign_customers_to_depots(Problem &pr) {
        // Stochastically assign customers to depots based on inverse distanse to depot.
        // TODO: Implement deterministically.
        std::vector<std::vector<int>> customers_on_depot(pr.get_num_depots());
        for (int cust=0; cust<pr.get_num_customers(); ++cust) {
            std::vector<std::vector<double>> depot_dist = pr.get_depot_distances();
            std::vector<double> depot_prob;
            for (double dist:depot_dist[cust]) {
                depot_prob.push_back(1/std::pow(dist, 10));
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
    static std::vector<T> get_subset(std::vector<T> &vals, std::vector<int> &idxs) {
        /**
         * Returns a subset of vals corresponding to the values of idxs.
         */
        std::vector<T> subset(idxs.size());
        for (int i=0; i<idxs.size(); ++i)
            subset[i]=vals[idxs[i]];
        return subset;
    }

    void setup_trips_forward(int depot_num, std::vector<int> customers, Problem &pr) {
        std::vector<int> order;
        std::vector<int> cur_trip;
        double trip_dist = 0;
        double trip_load = 0;
        double max_trip_dist = pr.get_max_length()[depot_num];
        int cur_cust = pr.get_num_customers() + depot_num;
        int num_vhcl = 1;
        std::vector<std::vector<double>> distances = pr.get_distances();
        while (customers.size() > 0) {
            std::vector<double> cust_distances = Individual::get_subset<double>(pr.get_distances()[cur_cust], customers);
            int closest_cust_pos = std::min_element(cust_distances.begin(), cust_distances.end()) - cust_distances.begin();
            int closest_cust = customers[closest_cust_pos];
            double closest_cust_dist = distances[cur_cust][closest_cust];
            double trip_and_ret_dist = trip_dist+closest_cust_dist+distances[closest_cust][pr.get_num_customers()+depot_num];
            double closest_cust_depot_dist = distances[closest_cust][pr.get_num_customers()+depot_num];
            // Check that the current trip is not too long.
            if (/*closest_cust_dist <= closest_cust_depot_dist && */(trip_and_ret_dist <= max_trip_dist || max_trip_dist==0)){
                order.push_back(closest_cust);
                cur_trip.push_back(closest_cust);
                trip_dist += distances[cur_cust][closest_cust];
                trip_load += pr.get_customer_load(closest_cust);
                customers.erase(customers.begin()+closest_cust_pos);
                cur_cust = closest_cust;
            } else if (trip_dist==0) {
                std::cout << "Trip could not be created as first point is too far away: " << closest_cust << " with distance "  << distances[closest_cust][pr.get_num_customers()+depot_num] << std::endl;
            } else {
                cur_cust = pr.get_num_customers() + depot_num;

                chromosome_trips[depot_num].push_back(cur_trip);
                cur_trip.clear();
                double ret_dist = distances[cur_cust][pr.get_num_customers()+depot_num];
                trip_dist+= ret_dist;
                trip_dists[depot_num].push_back(trip_dist);
                trip_dist = 0;

                trip_loads[depot_num].push_back(trip_load);
                trip_load = 0;

                ++num_vhcl;
            }
            if (trip_dist>0 && customers.size()==0) {
                chromosome_trips[depot_num].push_back(cur_trip);
                cur_trip.clear();
                double ret_dist = distances[cur_cust][pr.get_num_customers()+depot_num];
                trip_dist+= ret_dist;
                trip_dists[depot_num].push_back(trip_dist);
                trip_dist = 0;

                trip_loads[depot_num].push_back(trip_load);
                trip_load = 0;
            }
        }
        if (num_vhcl>pr.get_vhcl_pr_depot()) {
            std::cout << "Too many vehicles for depot: " << depot_num << " with "  << num_vhcl << " vehicles." << std::endl;
        }
    }

    void setup_trips_backward(int depot_num, std::vector<int> customers, Problem &pr) {
        std::vector<std::vector<double>> distances = pr.get_distances();
        // Is it better that the back is a new trip?
        if (chromosome_trips[depot_num].size()<pr.get_vhcl_pr_depot()) {
            int l_cust_lt = chromosome_trips[depot_num].back().back();
            int sl_cust_lt = chromosome_trips[depot_num].back()[chromosome_trips[depot_num].back().size()-2];
            double dist_l_dep_pt = distances[l_cust_lt][pr.get_num_customers()+depot_num];
            double dist_sl_l_pt = distances[sl_cust_lt][l_cust_lt];
            double dist_sl_dep_pt = distances[sl_cust_lt][pr.get_num_customers()+depot_num];

            if (2*dist_l_dep_pt+dist_sl_dep_pt<dist_sl_l_pt+dist_l_dep_pt) {
                chromosome_trips[depot_num].back().erase(chromosome_trips[depot_num].back().end()-1);
                chromosome_trips[depot_num].push_back(std::vector<int>{l_cust_lt});
                trip_dists[depot_num].back()+=dist_sl_dep_pt-dist_sl_l_pt-dist_l_dep_pt;
                trip_dists[depot_num].push_back(2*dist_l_dep_pt);
            }
        }
        for (int trip_num=chromosome_trips[depot_num].size()-1; trip_num>0; --trip_num) {
            /**
             * pt=previous trip, ct=current trip, l=last, sl=second last, f=first
             */
            // TODO: Edge case if pt consist of only one node.
            int l_cust_pt = chromosome_trips[depot_num][trip_num-1].back();
            int sl_cust_pt = chromosome_trips[depot_num][trip_num-1][chromosome_trips[depot_num][trip_num-1].size()-2];
            int f_cust_ct = chromosome_trips[depot_num][trip_num][0];
            double dist_l_dep_pt = distances[l_cust_pt][pr.get_num_customers()+depot_num];
            double dist_sl_l_pt = distances[sl_cust_pt][l_cust_pt];
            double dist_sl_dep_pt = distances[sl_cust_pt][pr.get_num_customers()+depot_num];
            double dist_l_pt_f_ct = distances[l_cust_pt][f_cust_ct];
            double dist_dep_f_ct = distances[pr.get_num_customers()+depot_num][f_cust_ct];

            if (dist_sl_dep_pt+dist_l_pt_f_ct<dist_sl_l_pt+dist_dep_f_ct) {
                chromosome_trips[depot_num][trip_num-1].erase(chromosome_trips[depot_num][trip_num-1].end()-1);
                chromosome_trips[depot_num][trip_num].insert(chromosome_trips[depot_num][trip_num].begin(), l_cust_pt);
                trip_dists[depot_num][trip_num-1]+=dist_sl_dep_pt-dist_sl_l_pt-dist_l_dep_pt;
                trip_dists[depot_num][trip_num]+=dist_l_dep_pt+dist_l_pt_f_ct-dist_dep_f_ct;
                ++trip_num;
            }
        }

    }

    void setup_trips(int depot_num, std::vector<int> customers, Problem &pr) {
        /**
         * Initializes the trips for a depot by greedily adding the customer that is closest.
         * A trip is ended if it is shorter to start a new or the maximum trip distance is exceeded.
         * The trip distance and trip load is calculated for each trip.
         */
        setup_trips_forward(depot_num, customers, pr);
        setup_trips_backward(depot_num, customers, pr);
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
