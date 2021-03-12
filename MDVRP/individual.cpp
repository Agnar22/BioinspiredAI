#include "individual.h"

// TODO: Implement remove_customer(int depot, int trip, int cust_pos, Problem &pr)
// and use it where remove_customers is not suited.

Individual::Individual(Problem &pr) {
    /**
    * Ways to initialize:
    *  - Deterministically/stochastically assign points to depots based on distance.
    *    - Create routes based on closest point deterministically/stochastically.
    **/
    tot_dist=0;
    initialize_chromosomes(pr);
}

void Individual::initialize_chromosomes(Problem &pr) {
    // Assign customers to depots.
    cust_on_depots = assign_customers_to_depots(pr, false);
    for (int depot=0; depot<cust_on_depots.size(); ++depot) {
        for (int cust=0; cust<cust_on_depots[depot].size(); ++cust) {
            std::cout << cust_on_depots[depot][cust] << " ";
        }
        std::cout << std::endl;
    }
    // Assign routes to depots.
    chromosome_trips.resize(pr.get_num_depots());
    trip_dists.resize(pr.get_num_depots());
    trip_loads.resize(pr.get_num_depots());
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
    double max_trip_dist = pr.get_max_length(depot_num);
    int cur_cust = pr.get_num_customers() + depot_num;
    int num_vhcl = 1;
    std::vector<std::vector<double>> distances = pr.get_distances();
    while (customers.size() > 0) {
        std::vector<double> cust_distances = Individual::get_subset(pr.get_distances()[cur_cust], customers);
        int closest_cust_pos = std::min_element(cust_distances.begin(), cust_distances.end()) - cust_distances.begin();
        int closest_cust = customers[closest_cust_pos];
        double closest_cust_dist = distances[cur_cust][closest_cust];
        double trip_and_ret_dist = trip_dist+closest_cust_dist+distances[closest_cust][pr.get_num_customers()+depot_num];
        double closest_cust_depot_dist = distances[closest_cust][pr.get_num_customers()+depot_num];
        double ret_dist = distances[cur_cust][pr.get_num_customers()+depot_num];
        // Check that the current trip is not too long.
        if (closest_cust_dist <= ret_dist+closest_cust_depot_dist &&
            (trip_and_ret_dist <= max_trip_dist || max_trip_dist==0) &&
            trip_load+pr.get_customer_load(closest_cust) <= pr.get_max_load(depot_num)){
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
            tot_dist += trip_dists[depot_num].back();

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
            tot_dist += trip_dists[depot_num].back();
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

double Individual::calculate_trip_distance(std::vector<int> &customers, int depot, Problem &pr) {
    if (customers.size()==0)
        throw std::runtime_error("No customers when calculating trip distance.");

    int first_cust = customers[0];
    int last_cust = customers[customers.size()-1];
    int num_cust = pr.get_num_customers();
    double sum_dist = pr.get_distance(num_cust+depot, first_cust)+pr.get_distance(last_cust, num_cust+depot);

    for (int pos=0; pos<customers.size()-1; ++pos)
        sum_dist+=pr.get_distance(customers[pos], customers[pos+1]);
    return sum_dist;
}

double Individual::get_fitness() const {
    return 1.0/(double)(tot_dist);
}

void Individual::remove_customers(std::vector<int> &custs, Problem &pr) {
    // X chromosomes - not used, remove?
    // X cust_on_depot
    // X chromosome_trips
    // X trip_dists
    // X trip_loads
    // X fitness
    for (int cust:custs) {
        int rmd_pos = remove_from_2d_vector(cust_on_depots, cust, false);
        int num_trips_on_depot = chromosome_trips[rmd_pos].size();
        int trip_pos = remove_from_2d_vector(chromosome_trips[rmd_pos], cust, true);
        if (num_trips_on_depot==chromosome_trips[rmd_pos].size()) {
            double trip_dist = calculate_trip_distance(chromosome_trips[rmd_pos][trip_pos], rmd_pos, pr);
            double diff_trip_dist = trip_dist-trip_dists[rmd_pos][trip_pos];
            trip_dists[rmd_pos][trip_pos] = trip_dist;
            trip_loads[rmd_pos][trip_pos] -= pr.get_customer_load(cust);
            tot_dist += diff_trip_dist;
        } else {
            tot_dist -= trip_dists[rmd_pos][trip_pos];
            trip_dists[rmd_pos].erase(trip_dists[rmd_pos].begin()+trip_pos);
            trip_loads[rmd_pos].erase(trip_loads[rmd_pos].begin()+trip_pos);
        }
    }
}

void Individual::reversal_mutation(int depot, Problem &pr) {
    // TODO: The original paper proposes mutation over trips, this is not done here.
    // TODO: Can be optimized by only calculating the distances at the end, as these are the ones that are changing.
    // Reverses a random part of a trip. Ensures that the new trip does not violate trip length constraints.

    int from, to, trip = rand()%chromosome_trips[depot].size();
    do {
        from = rand()%chromosome_trips[depot].size();
        to = rand()%chromosome_trips[depot].size();
        if (from>to)
            std::swap(from, to);
    } while (from==to);

    double prev_trip_dist = trip_dists[depot][trip];
    std::reverse(chromosome_trips[depot][trip].begin()+from, chromosome_trips[depot][trip].begin()+to);
    double new_trip_dist = calculate_trip_distance(chromosome_trips[depot][trip], depot, pr);
    if (new_trip_dist>pr.get_max_length(depot))
        std::reverse(chromosome_trips[depot][trip].begin()+from, chromosome_trips[depot][trip].begin()+to);
    else {
        trip_dists[depot][trip]=new_trip_dist;
        tot_dist+=new_trip_dist-prev_trip_dist;
    }
}

void Individual::re_routing_mutation(int depot, Problem &pr) {
    int trip = rand()%chromosome_trips[depot].size();
    int cust_pos = rand()%chromosome_trips[depot][trip].size();
    int cust = chromosome_trips[depot][trip][cust];
    auto insert_costs_and_position = find_insert_costs(cust, depot, pr);
    std::vector<double> insert_costs = insert_costs_and_position.first;
    std::vector<std::pair<int, int>> positions = insert_costs_and_position.second;

    int best_idx = std::min_element(insert_costs.begin(), insert_costs.end())-insert_costs.begin();
    std::pair<int, int> best_pos = positions[best_idx];
    std::vector<int> cust_vec = {cust};
    remove_customers(cust_vec, pr);
    insert_customer(depot, trip, cust_pos, cust, pr);
}

void Individual::swapping_mutation(int depot, Problem &pr) {
    int trip1 = rand()&chromosome_trips[depot].size();
    int trip2 = rand()&chromosome_trips[depot].size();
    int cust1_pos, cust2_pos;
    do {
        cust1_pos = rand()%chromosome_trips[depot][trip1].size();
        cust2_pos = rand()%chromosome_trips[depot][trip2].size();
    } while (trip1==trip2 && cust1_pos==cust2_pos);
    std::vector<int> cust1 = {chromosome_trips[depot][trip1][cust1_pos]};
    std::vector<int> cust2 = {chromosome_trips[depot][trip2][cust2_pos]};

    remove_customers(cust1, pr);
    remove_customers(cust2, pr);

    insert_customer(depot, trip2, cust2_pos, cust1[0], pr);
    insert_customer(depot, trip1, cust1_pos, cust2[0], pr);

    if (chromosome_trips[depot][trip1][cust1_pos]!=cust2[0] || chromosome_trips[depot][trip2][cust2_pos]!=cust1[0]) {
        std::cout << "Swapping mutation was not possible, undoing the mutation." << std::endl;
        remove_customers(cust1, pr);
        remove_customers(cust2, pr);

        insert_customer(depot, trip1, cust1_pos, cust1[0], pr);
        insert_customer(depot, trip2, cust2_pos, cust2[0], pr);
    }
}

std::pair<std::vector<double>, std::vector<std::pair<int, int>>> Individual::find_insert_costs(int cust, int depot, Problem &pr){
    std::vector<std::pair<int, int>> insert_positions;
    std::vector<double> insert_costs;
    double max_length = pr.get_max_length(depot);
    double max_load = pr.get_max_load(depot);
    int num_cust = pr.get_num_customers();

    for (int trip=0; trip<chromosome_trips[depot].size(); ++trip) {
        for (int pos_in_trip=0; pos_in_trip<chromosome_trips[depot].size()+1; ++pos_in_trip) {
            std::vector<int> cur_trip = chromosome_trips[depot][trip];
            int cust_before = pos_in_trip==0 ? depot+num_cust : cur_trip[pos_in_trip-1];
            int cust_after = pos_in_trip==chromosome_trips[depot].size() ? depot+num_cust : cur_trip[pos_in_trip];

            double insert_cost = marginal_cost(cust_before, cust_after, cust, pr);
            double trip_length = insert_cost+trip_dists[depot][trip];
            double trip_load = pr.get_customer_load(cust)+trip_loads[depot][trip];

            if (trip_length <= max_length && trip_load <= max_load) {
                insert_positions.push_back(std::make_pair(trip, pos_in_trip));
                insert_costs.push_back(insert_cost);
            }
        }
    }
    return std::make_pair(insert_costs, insert_positions);
}

void Individual::insert_stochastically(int cust, double prob_greedy, int depot, Problem &pr) {
    std::vector<std::pair<int, int>> insert_positions;
    std::vector<double> insert_costs;

    auto insert_costs_and_position = find_insert_costs(cust, depot, pr);
    insert_costs = insert_costs_and_position.first;
    insert_positions = insert_costs_and_position.second;

    if (insert_positions.size()==0) {
        if (chromosome_trips[depot].size()>=pr.get_vhcl_pr_depot())
            throw std::runtime_error("Could not create new tour from customer.");
        // TODO: Add new route...
    } else {
        int chosen_insert_pos;
        if (true) { // TODO: Generate random number.
            // TODO: Can optimize by keeping the minimum stored as a variable.
            chosen_insert_pos = std::min_element(insert_costs.begin(), insert_costs.end())-insert_costs.begin();
        } else {
            chosen_insert_pos = rand()%insert_positions.size();
        }
        std::pair<int, int> insert_pos = insert_positions[chosen_insert_pos];
        insert_customer(depot, insert_pos.first, insert_pos.second, cust, pr);
    }
}

double Individual::marginal_cost(int cust_from, int cust_to, int cust_added_between, Problem &pr) {
    return pr.get_distance(cust_from, cust_added_between) + pr.get_distance(cust_added_between, cust_to) - pr.get_distance(cust_from, cust_to);
}

void Individual::insert_customer(int depot, int trip, int pos_in_trip, int cust, Problem &pr) {
    // TODO: Check constraints on trip length and load.
    // TODO: Check that arguments are valid.
    int num_cust = pr.get_num_customers();
    std::vector<int> cur_trip = chromosome_trips[depot][trip];
    int cust_before = pos_in_trip==0 ? depot+num_cust : cur_trip[pos_in_trip-1];
    int cust_after = pos_in_trip==chromosome_trips[depot].size() ? depot+num_cust : cur_trip[pos_in_trip];

    cust_on_depots[depot].push_back(cust);
    chromosome_trips[depot][trip].insert(chromosome_trips[depot][trip].begin()+pos_in_trip, cust);
    trip_dists[depot][trip] += marginal_cost(cust_before, cust_after, cust, pr);
    trip_loads[depot][trip] += pr.get_customer_load(cust);
    tot_dist += marginal_cost(cust_before, cust_after, cust, pr);
}

int Individual::remove_from_2d_vector(std::vector<std::vector<int>> &nested_vec, int cust, bool remove_if_empty) {
    // TODO: Can be optimised to search the most promising depot first.
    for (int idx=0; idx<nested_vec.size(); ++idx) {
        auto pos = std::find(nested_vec[idx].begin(), nested_vec[idx].end(), cust);
        if (pos!=nested_vec[idx].end()) {
            nested_vec[idx].erase(pos);
            if (nested_vec[idx].size() == 0 && remove_if_empty) {
                nested_vec.erase(nested_vec.begin()+idx);
            }
            return idx;
        }
    }
    throw std::runtime_error("Customer not found in 2d_vector when removing.");
}

std::vector<double> Individual::get_subset(std::vector<double> &vals, std::vector<int> &idxs) {
    std::vector<double> subset(idxs.size());
    for (int i=0; i<idxs.size(); ++i)
        subset[i]=vals[idxs[i]];
    return subset;
}