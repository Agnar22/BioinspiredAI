#include "individual.h"
#include "file.h"
#include "config.h"


Individual::Individual(Problem &pr) {
    tot_dist=0;
    initialize_chromosomes(pr);
}

void Individual::initialize_chromosomes(Problem &pr) {
    bool initialized = false;
    while (!initialized) {
        try {
            cust_on_depots = assign_customers_to_depots(pr, STOCHASTIC_ASSIGNMENT);

            // Assign routes to depots.
            chromosome_trips.resize(pr.get_num_depots());
            trip_dists.resize(pr.get_num_depots());
            trip_loads.resize(pr.get_num_depots());
            for (int depot=0; depot<pr.get_num_depots(); ++depot) {
                setup_trips(depot, cust_on_depots[depot], pr);
            }
            initialized = true;
        } catch (std::invalid_argument e) {
            std::cout << e.what() << std::endl;
            file::write_solution(*this, "p01.res");
            reset();
        }
    }
}

void Individual::reset() {
    cust_on_depots.clear();
    chromosome_trips.clear();
    trip_dists.clear();
    trip_loads.clear();
    tot_dist=0;
}

std::vector<std::vector<int>> Individual::assign_customers_to_depots(Problem &pr, bool stoch) {
    std::vector<std::vector<int>> customers_on_depot(pr.get_num_depots());
    std::vector<int> depot_loads(pr.get_num_depots(), 0);
    for (int cust=0; cust<pr.get_num_customers(); ++cust) {
        std::vector<std::vector<double>> depot_dist = pr.get_depot_distances();
        // Stochastically assign customers to depots based on inverse distanse to depot.
        if (stoch) {
            std::vector<double> depot_prob;
            for (double dist:depot_dist[cust]) {
                depot_prob.push_back(1/std::pow(dist, POW_INV_DEPOT_DIST));
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
                    depot_loads[pos] += pr.get_customer_load(cust);
                    break;
                }
            }
        } else {
            // Deterministically assign customer to closest depot.
            int closest_depot = std::min_element(depot_dist[cust].begin(), depot_dist[cust].end()) - depot_dist[cust].begin();
            customers_on_depot[closest_depot].push_back(cust);
            depot_loads[closest_depot] += pr.get_customer_load(cust);
        }
    }

    // Make sure that no depot gets assigned more total load than its vehicles can carry.
    bool sufficient_capacity = false;
    while (!sufficient_capacity) {
        sufficient_capacity = true;
        for (int depot=0; depot<pr.get_num_depots(); ++depot) {
            if (depot_loads[depot]<=pr.get_max_load(depot)*pr.get_vhcl_pr_depot())
                continue;
            sufficient_capacity = false;

            std::vector<int> depot_for_cust(pr.get_num_customers());
            for (int cur_depot=0; cur_depot<pr.get_num_depots(); ++cur_depot) {
                for (int cur_cust:customers_on_depot[cur_depot]) {
                    depot_for_cust[cur_cust] = cur_depot;
                }
            }
            std::vector<int> all_customers(pr.get_num_customers());
            std::iota(all_customers.begin(), all_customers.end(), 0);
            std::vector<int> other_custs = remove_subset<int>(all_customers, customers_on_depot[depot]);
            std::vector<double> inv_dist;
            std::vector<int> chosen_depot;
            for (int cust:customers_on_depot[depot]) {
                std::vector<double> distances_to_other_custs = get_subset<double>(pr.get_distances()[cust], other_custs);
                int closest_cust_pos = std::min_element(distances_to_other_custs.begin(), distances_to_other_custs.end())-distances_to_other_custs.begin();
                int closest_cust = other_custs[closest_cust_pos];
                inv_dist.push_back(1.0/std::pow(pr.get_distance(cust, closest_cust), POW_RE_ASSIGN_INV_DEPOT_DIST));
                chosen_depot.push_back(depot_for_cust[closest_cust]);
            }

            while (depot_loads[depot]>pr.get_max_load(depot)*pr.get_vhcl_pr_depot()) {
                int chosen_pos = random_choice(inv_dist);
                int chosen_cust = customers_on_depot[depot][chosen_pos];
                int next_depot = chosen_depot[chosen_pos];

                customers_on_depot[depot].erase(customers_on_depot[depot].begin()+chosen_pos);
                inv_dist.erase(inv_dist.begin()+chosen_pos);
                chosen_depot.erase(chosen_depot.begin()+chosen_pos);
                customers_on_depot[next_depot].push_back(chosen_cust);
                depot_loads[depot] -= pr.get_customer_load(chosen_cust);
                depot_loads[next_depot] += pr.get_customer_load(chosen_cust);
            }
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
    if (GREEDY_PROB>(double)(rand()) / (double)(RAND_MAX))
        setup_trips_forward(0, depot_num, customers, pr);
    else
        c_and_w_algorithm(0, depot_num, customers, pr);
}


void Individual::setup_trips_forward(int attempt_num, int depot_num, std::vector<int> customers_to_order, Problem &pr) {
    std::vector<int> order;
    std::vector<int> cur_trip;
    std::vector<int> customers = customers_to_order;
    double trip_dist = 0;
    double trip_load = 0;
    double max_trip_dist = pr.get_max_length(depot_num);
    int cur_cust = pr.get_num_customers() + depot_num;
    int num_vhcl = 1;
    std::vector<std::vector<double>> distances = pr.get_distances();
    int same_counter = 0;
    int not_same_counter = 0;
    while (customers.size() > 0) {
        std::vector<double> cust_distances = get_subset(pr.get_distances()[cur_cust], customers);
        int closest_cust_pos_deterministically = std::min_element(cust_distances.begin(), cust_distances.end()) - cust_distances.begin();
        std::vector<double> unscaled_cust_probs;
        for (double dist:cust_distances)
            unscaled_cust_probs.push_back(1.0/std::pow(dist, POW_INV_CUST_DIST));
        int closest_cust_pos = random_choice(unscaled_cust_probs);
        same_counter+= closest_cust_pos == closest_cust_pos_deterministically;
        not_same_counter+= closest_cust_pos != closest_cust_pos_deterministically;
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
            throw std::invalid_argument("Trip could not be created as first point is too far away");
        } else {
            if (trip_load+pr.get_customer_load(closest_cust) > pr.get_max_load(depot_num))
                std::cout << "Too much load: " << trip_load+pr.get_customer_load(closest_cust) << " " << pr.get_max_load(depot_num) << std::endl;
            else
                std::cout << "Too long: " << trip_and_ret_dist << " " << max_trip_dist << std::endl;

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
    //std::cout << same_counter << " " << not_same_counter << std::endl;
    if (num_vhcl>pr.get_vhcl_pr_depot()) {
        std::cout << "Too many vehicles for depot: " << depot_num << " with "  << num_vhcl << " vehicles." << std::endl;
        if (attempt_num==RETRY_ATTEMPTS)
            throw std::invalid_argument("Failed more than three times to set up. Giving up...");
        chromosome_trips[depot_num].clear();
        tot_dist -= std::accumulate(trip_dists[depot_num].begin(), trip_dists[depot_num].end(), 0.0);
        trip_dists[depot_num].clear();
        trip_loads[depot_num].clear();
        setup_trips_forward(attempt_num+1, depot_num, customers_to_order, pr);
    } else {
        std::cout << "Completed setting up for depot: " << depot_num << " with "  << num_vhcl << " vehicles." << std::endl;
    }
}

std::pair<int, int> find_in_2d_vector(std::vector<std::vector<int>> &vec, int val) {
    for (int x=0; x<vec.size(); ++x)
        for (int y=0; y<vec[x].size(); ++y)
            if (vec[x][y]==val)
                return std::make_pair(x, y);
    std::invalid_argument("Val was not found in vector.");
}

void Individual::c_and_w_algorithm(int attempt_num, int depot_num, std::vector<int> customers_to_order, Problem &pr) {
    std::vector<int> num_ends(pr.get_num_customers(), 2);
    std::vector<double> dists;
    std::vector<std::vector<int>> trips;
    std::vector<int> load_for_trips;
    std::priority_queue<std::pair<double, std::pair<int, int>>> pq;
    for (int cust:customers_to_order) {
        dists.push_back(2*pr.get_distance(cust, depot_num+pr.get_num_customers()));
        tot_dist += 2*pr.get_distance(cust, depot_num+pr.get_num_customers());
        trips.push_back(std::vector<int>{cust});
        load_for_trips.push_back(pr.get_customer_load(cust));
        for (int cust2:customers_to_order) {
            if (cust>=cust2)
                continue;
            int dep = depot_num + pr.get_num_customers();
            double savings = pr.get_distance(cust, dep) + pr.get_distance(cust2, dep) - pr.get_distance(cust, cust2);
            pq.push(std::make_pair(savings, std::make_pair(cust, cust2)));
        }
    }
    while (!pq.empty() && (trips.size()>pr.get_vhcl_pr_depot() || pq.top().first>=0)) {
        auto top = pq.top();
        int first_cust = top.second.first;
        int second_cust = top.second.second;
        pq.pop();
        // Check that both are ends of trips.
        if (!(num_ends[first_cust]>0 && num_ends[second_cust]>0))
            continue;
        auto first_pos = find_in_2d_vector(trips, first_cust);
        auto second_pos = find_in_2d_vector(trips, second_cust);
        if (first_pos.first==second_pos.first)
            continue;
        // If it is possible to merge top.first and top.second wrt. load and length: merge
        if (load_for_trips[first_pos.first]+load_for_trips[second_pos.first]>pr.get_max_load(depot_num) ||
            (dists[first_pos.first]+dists[second_pos.first]-top.first>pr.get_max_length(depot_num) && pr.get_max_length(depot_num) > 0))
            continue;
        if (DROP_PROB>(double)(rand()) / (double)(RAND_MAX))
            continue;
        if (first_pos.second==0)
            std::reverse(trips[first_pos.first].begin(), trips[first_pos.first].end());
        if (second_pos.second!=0)
            std::reverse(trips[second_pos.first].begin(), trips[second_pos.first].end());
        num_ends[first_cust]-=1;
        num_ends[second_cust]-=1;
        trips[first_pos.first].insert(trips[first_pos.first].end(), trips[second_pos.first].begin(), trips[second_pos.first].end());
        trips.erase(trips.begin()+second_pos.first);
        dists[first_pos.first] += dists[second_pos.first]-top.first;
        dists.erase(dists.begin()+second_pos.first);
        load_for_trips[first_pos.first] += load_for_trips[second_pos.first];
        load_for_trips.erase(load_for_trips.begin()+second_pos.first);
        tot_dist -= top.first;
    }
    if (trips.size()>pr.get_vhcl_pr_depot()) {
        if (attempt_num==RETRY_ATTEMPTS)
            throw std::invalid_argument("Failed more than three times to set up. Giving up...");
        tot_dist -= std::accumulate(dists.begin(), dists.end(), 0.0);
        c_and_w_algorithm(attempt_num+1, depot_num, customers_to_order, pr);
    } else {
        chromosome_trips[depot_num] = trips;
        trip_dists[depot_num] = dists;
        trip_loads[depot_num] = load_for_trips;
    }
}



int Individual::random_choice(std::vector<double> &unscaled_probs) {
    double sum=0;
    std::for_each(unscaled_probs.begin(), unscaled_probs.end(), [&] (double n) {sum += n;});
    double prob = (double)(rand()) / (double)(RAND_MAX);
    double cum_sum = 0;
    for (int pos=0; pos<unscaled_probs.size(); ++pos) {
        cum_sum+=unscaled_probs[pos]/sum;
        if (cum_sum>=prob) {
            return pos;
        }
    }
    throw std::invalid_argument("Not possible too choose stochastically.");
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
    // TODO: Implement remove_customer(int depot, int trip, int cust_pos, Problem &pr)
    // and use it where remove_customers is not suited.
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
    if (chromosome_trips[depot][trip].size()==1)
        return;
    do {
        from = rand()%chromosome_trips[depot][trip].size();
        to = rand()%chromosome_trips[depot][trip].size();
        if (from>to)
            std::swap(from, to);
    } while (from==to);

    double prev_trip_dist = trip_dists[depot][trip];
    std::reverse(chromosome_trips[depot][trip].begin()+from, chromosome_trips[depot][trip].begin()+to+1);
    double new_trip_dist = calculate_trip_distance(chromosome_trips[depot][trip], depot, pr);
    if (new_trip_dist>pr.get_max_length(depot) && pr.get_max_length(depot)!=0)
        std::reverse(chromosome_trips[depot][trip].begin()+from, chromosome_trips[depot][trip].begin()+to+1);
    else {
        trip_dists[depot][trip]=new_trip_dist;
        tot_dist+=new_trip_dist-prev_trip_dist;
    }
}

void Individual::re_routing_mutation(int depot, Problem &pr) {
    int trip = rand()%chromosome_trips[depot].size();
    int cust_pos = rand()%chromosome_trips[depot][trip].size();
    int cust = chromosome_trips[depot][trip][cust_pos];

    std::vector<int> cust_vec = {cust};
    remove_customers(cust_vec, pr);

    auto insert_costs_and_position = find_insert_costs(cust, depot, pr);
    std::vector<double> insert_costs = insert_costs_and_position.first;
    std::vector<std::pair<int, int>> positions = insert_costs_and_position.second;

    if (insert_costs.size()!=0) {
        int best_idx = std::min_element(insert_costs.begin(), insert_costs.end())-insert_costs.begin();
        std::pair<int, int> best_pos = positions[best_idx];
        insert_customer(depot, best_pos.first, best_pos.second, cust, pr);
    } else {
        cust_on_depots[depot].push_back(cust);
        chromosome_trips[depot].insert(chromosome_trips[depot].begin()+trip, std::vector<int>{cust});
        double trip_dist = 2*pr.get_distance(cust, pr.get_num_customers()+depot);
        trip_dists[depot].insert(trip_dists[depot].begin()+trip, trip_dist);
        trip_loads[depot].insert(trip_loads[depot].begin()+trip, pr.get_customer_load(cust));
        tot_dist+=trip_dist;
    }
}

void Individual::swapping_mutation(int depot, Problem &pr) {
    int trip1 = rand()%chromosome_trips[depot].size();
    int trip2 = rand()%chromosome_trips[depot].size();
    int cust1_pos, cust2_pos;
    if (trip1 == trip2 && chromosome_trips[depot][trip1].size()<=2)
        return;
    if (chromosome_trips[depot][trip1].size()==1 || chromosome_trips[depot][trip2].size()==1)
        return;

    do {
        cust1_pos = rand()%chromosome_trips[depot][trip1].size();
        cust2_pos = rand()%chromosome_trips[depot][trip2].size();
    } while (trip1==trip2 && cust1_pos==cust2_pos);
    if (trip1==trip2 && cust1_pos > cust2_pos)
        std::swap(cust1_pos, cust2_pos);

    std::vector<int> cust1 = {chromosome_trips[depot][trip1][cust1_pos]};
    std::vector<int> cust2 = {chromosome_trips[depot][trip2][cust2_pos]};

    remove_customers(cust1, pr);
    remove_customers(cust2, pr);

    bool inserted2=false, inserted1=false;
    inserted2 = insert_customer(depot, trip1, cust1_pos, cust2[0], pr);
    if (inserted2)
        inserted1 = insert_customer(depot, trip2, cust2_pos, cust1[0], pr);
    if (!inserted1) {
        if (inserted2)
            remove_customers(cust2, pr);
        bool undone1 = insert_customer(depot, trip1, cust1_pos, cust1[0], pr);
        bool undone2 = insert_customer(depot, trip2, cust2_pos, cust2[0], pr);
        if (!undone1 || !undone2)
            throw std::runtime_error("Error when inserting customers back in.");
    }
}

void Individual::inter_depot_mutation(int depot, Problem &pr) {
    double lowest_dist=INFINITY;
    int lowest_cust, to_depot, to_trip, insert_pos;
    for (int from_trip=0; from_trip<chromosome_trips[depot].size(); ++from_trip) {
        std::vector<int> cur_trip = chromosome_trips[depot][from_trip];
        for (int from_pos=0; from_pos<chromosome_trips[depot][from_trip].size(); ++from_trip) {
            int cust = cur_trip[from_pos];
            int cust_before = from_pos==0 ? depot+pr.get_num_customers() : cur_trip[from_pos-1];
            int cust_after = from_pos==chromosome_trips[depot][from_trip].size() ? depot+pr.get_num_customers() : cur_trip[from_pos];
            double marg_cost = marginal_cost(cust_before, cust_after, cust, pr);
            for (int other_depot=0; other_depot<pr.get_num_depots(); ++other_depot) {
                if (other_depot!=depot) {
                    auto insert_costs_and_pos = find_insert_costs(cust, other_depot, pr);
                    auto insert_costs = insert_costs_and_pos.first;
                    if (insert_costs.size()==0)
                        continue;
                    int best_pos = std::min_element(insert_costs.begin(), insert_costs.end())-insert_costs.begin();
                    double cur_lowest_dist = insert_costs[best_pos] - marg_cost;

                    if (cur_lowest_dist<lowest_dist) {
                        lowest_dist = cur_lowest_dist;
                        lowest_cust = cust;
                        to_depot = other_depot;
                        to_trip = insert_costs_and_pos.second[best_pos].first;
                        insert_pos = insert_costs_and_pos.second[best_pos].second;
                        if (INTER_DEPOT_RANDOM_BEST_PROB > (double)(rand()) / (double)(RAND_MAX))
                            goto finished_searching;
                    }
                }
            }
        }
    }
    finished_searching:
    if (lowest_dist!=INFINITY) {
        std::vector<int> chosen_cust{lowest_cust};
        remove_customers(chosen_cust, pr);
        insert_customer(to_depot, to_trip, insert_pos, lowest_cust, pr);
    }
}

std::pair<std::vector<double>, std::vector<std::pair<int, int>>> Individual::find_insert_costs(int cust, int depot, Problem &pr){
    std::vector<std::pair<int, int>> insert_positions;
    std::vector<double> insert_costs;
    double max_length = pr.get_max_length(depot);
    double max_load = pr.get_max_load(depot);
    int num_cust = pr.get_num_customers();

    for (int trip=0; trip<chromosome_trips[depot].size(); ++trip) {
        for (int pos_in_trip=0; pos_in_trip <= chromosome_trips[depot][trip].size(); ++pos_in_trip) {
            std::vector<int> cur_trip = chromosome_trips[depot][trip];
            int cust_before = pos_in_trip==0 ? depot+num_cust : cur_trip[pos_in_trip-1];
            int cust_after = pos_in_trip==chromosome_trips[depot][trip].size() ? depot+num_cust : cur_trip[pos_in_trip];

            double insert_cost = Individual::marginal_cost(cust_before, cust_after, cust, pr);
            double trip_length = insert_cost+trip_dists[depot][trip];
            double trip_load = pr.get_customer_load(cust)+trip_loads[depot][trip];

            if ((trip_length <= max_length || max_length==0) && trip_load <= max_load) {
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
        cust_on_depots[depot].push_back(cust);
        chromosome_trips[depot].push_back(std::vector<int>{cust});
        trip_dists[depot].push_back(2*pr.get_depot_distances()[cust][depot]);
        trip_loads[depot].push_back(pr.get_customer_load(cust));
        tot_dist += trip_dists[depot].back();
    } else {
        int chosen_insert_pos;
        if (prob_greedy>(double)(rand()) / (double)(RAND_MAX)) {
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

bool Individual::insert_customer(int depot, int trip, int pos_in_trip, int cust, Problem &pr) {
    int num_cust = pr.get_num_customers();
    std::vector<int> cur_trip = chromosome_trips[depot][trip];
    int cust_before = pos_in_trip==0 ? depot+num_cust : cur_trip[pos_in_trip-1];
    int cust_after = pos_in_trip==cur_trip.size() ? depot+num_cust : cur_trip[pos_in_trip];

    double marginal_cost = Individual::marginal_cost(cust_before, cust_after, cust, pr);

    if (marginal_cost + trip_dists[depot][trip] > pr.get_max_length(depot) && pr.get_max_length(depot)!=0)
        return false;
    if (pr.get_customer_load(cust) + trip_loads[depot][trip] > pr.get_max_load(depot))
        return false;

    cust_on_depots[depot].push_back(cust);
    if (pos_in_trip < chromosome_trips[depot][trip].size())
        chromosome_trips[depot][trip].insert(chromosome_trips[depot][trip].begin()+pos_in_trip, cust);
    else
        chromosome_trips[depot][trip].push_back(cust);
    trip_dists[depot][trip] += marginal_cost;
    trip_loads[depot][trip] += pr.get_customer_load(cust);
    tot_dist += marginal_cost;
    return true;
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
