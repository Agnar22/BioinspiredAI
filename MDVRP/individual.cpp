#include "individual.h"

Individual::Individual(Problem pr) {
    /**
    * Ways to initialize:
    *  - Deterministically/stochastically assign points to depots based on distance.
    *    - Create routes based on closest point deterministically/stochastically.
    **/

    initialize_chromosomes(pr);
}

void Individual::initialize_chromosomes(Problem &pr) {
    chromosomes.resize(pr.get_num_depots());
    // Assign customers to depots.
    std::cout << "assigning cust to depots" << std::endl;
    cust_on_depots = assign_customers_to_depots(pr, false);
    std::cout << "cust on depots size: " << cust_on_depots.size() << std::endl;
    for (int depot=0; depot<cust_on_depots.size(); ++depot) {
        for (int cust=0; cust<cust_on_depots[depot].size(); ++cust) {
            std::cout << cust_on_depots[depot][cust] << " ";
        }
        std::cout << std::endl;
    }
    // Assign routes to depots.
    chromosome_trips.resize(chromosomes.size());
    trip_dists.resize(chromosomes.size());
    trip_loads.resize(chromosomes.size());
    for (int depot=0; depot<pr.get_num_depots(); ++depot) {
        setup_trips(depot, cust_on_depots[depot], pr);
    }
}

std::vector<std::vector<int>> Individual::assign_customers_to_depots(Problem &pr, bool stoch) {
    // Stochastically assign customers to depots based on inverse distanse to depot.
    std::vector<std::vector<int>> customers_on_depot(pr.get_num_depots());
    for (int cust=0; cust<pr.get_num_customers(); ++cust) {
        std::vector<std::vector<double>> depot_dist = pr.get_depot_distances();
        if (stoch) {
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
        } else {
            int closest_depot = std::min_element(depot_dist[cust].begin(), depot_dist[cust].end()) - depot_dist[cust].begin();
            customers_on_depot[closest_depot].push_back(cust);
        }
    }
    return customers_on_depot;
}

void Individual::setup_trips(int depot_num, std::vector<int> customers, Problem &pr) {
    /**
     * Initializes the trips for a depot by greedily adding the customer that is closest.
     * A trip is ended if it is shorter to start a new or the maximum trip distance is exceeded.
     * The trip distance and trip load is calculated for each trip.
     */
    setup_trips_forward(depot_num, customers, pr);
    setup_trips_backward(depot_num, customers, pr);
}


void Individual::setup_trips_forward(int depot_num, std::vector<int> customers, Problem &pr) {
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
        double ret_dist = distances[cur_cust][pr.get_num_customers()+depot_num];
        // Check that the current trip is not too long.
        if (closest_cust_dist <= ret_dist+closest_cust_depot_dist && (trip_and_ret_dist <= max_trip_dist || max_trip_dist==0)){
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

void Individual::setup_trips_backward(int depot_num, std::vector<int> customers, Problem &pr) {
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
        if (chromosome_trips[depot_num][trip_num-1].size()!=1) {
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
        } else {
            int l_cust_pt = chromosome_trips[depot_num][trip_num-1].back();
            int f_cust_ct = chromosome_trips[depot_num][trip_num][0];
            double dist_l_dep_pt = distances[l_cust_pt][pr.get_num_customers()+depot_num];
            double dist_l_pt_f_ct = distances[l_cust_pt][f_cust_ct];
            double dist_dep_f_ct = distances[pr.get_num_customers()+depot_num][f_cust_ct];

            if (dist_l_pt_f_ct < dist_l_dep_pt+dist_dep_f_ct) {
                chromosome_trips[depot_num][trip_num].insert(chromosome_trips[depot_num][trip_num].begin(), l_cust_pt);
                trip_dists[depot_num][trip_num]+=dist_l_dep_pt+dist_l_pt_f_ct-dist_dep_f_ct;
                chromosome_trips[depot_num].erase(chromosome_trips[depot_num].begin()+(trip_num)-1);
                trip_dists.erase(trip_dists.begin()+(trip_num-1));
                ++trip_num;
            }
        }
    }
}
template <class T>
std::vector<T> Individual::get_subset(std::vector<T> &vals, std::vector<int> &idxs) {
    std::vector<T> subset(idxs.size());
    for (int i=0; i<idxs.size(); ++i)
        subset[i]=vals[idxs[i]];
    return subset;
}