#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

#include <algorithm>
#include <vector>
#include <utility>
#include <numeric>
#include <queue>
#include "problem.h"

class Individual {
    public:
        std::vector<std::vector<int>> cust_on_depots;
        std::vector<std::vector<std::vector<int>>> chromosome_trips;
        std::vector<std::vector<double>> trip_dists;
        std::vector<std::vector<int>> trip_loads;
        double tot_dist;

        Individual(){};
        Individual(Problem&);
        void initialize_chromosomes(Problem&);
        void reset();
        std::vector<std::vector<int>> assign_customers_to_depots(Problem&, bool);
        void setup_trips(int, std::vector<int>, Problem&);
        void setup_trips_forward(int, int, std::vector<int>, Problem&);
        void c_and_w_algorithm(int, int, std::vector<int>, Problem&);
        int random_choice(std::vector<double>&);
        static double calculate_trip_distance(std::vector<int>&, int, Problem&);
        double get_fitness() const;
        void remove_customers(std::vector<int>&, Problem&);
        int remove_from_2d_vector(std::vector<std::vector<int>>&, int, bool);
        void reversal_mutation(int, Problem&);
        void re_routing_mutation(int, Problem&);
        void swapping_mutation(int, Problem&);
        void inter_depot_mutation(int, Problem&);
        std::pair<std::vector<double>, std::vector<std::pair<int, int>>> find_insert_costs(int, int, Problem&);
        void insert_stochastically(int, double, int, Problem&);
        static double marginal_cost(int, int, int, Problem&);
        bool insert_customer(int, int, int, int, Problem&);
};

inline bool operator<(const Individual &a, const Individual &b) {
    if (a.tot_dist<1 || a.tot_dist>100000)
        throw std::runtime_error("Individual fitness is not within plausible range.");
    if (b.tot_dist<1 || b.tot_dist>100000)
        throw std::runtime_error("Individual fitness is not within plausible range.");
    return a.get_fitness()<b.get_fitness();
};

inline bool operator==(const Individual &a, const Individual &b) {
    if (a.chromosome_trips.size() != b.chromosome_trips.size())
        return false;

    for (int depot=0; depot<a.chromosome_trips.size(); ++depot) {
        if (a.chromosome_trips[depot].size() != b.chromosome_trips[depot].size())
            return false;
        for (int trip=0; trip<a.chromosome_trips[depot].size(); ++trip) {
            if (a.chromosome_trips[depot][trip].size() != b.chromosome_trips[depot][trip].size())
                return false;
            for (int pos=0; pos<a.chromosome_trips[depot][trip].size(); ++pos)
                if (a.chromosome_trips[depot][trip][pos] != b.chromosome_trips[depot][trip][pos])
                    return false;
        }
    }
    return true;
}

template<typename T>
std::vector<T> get_subset(std::vector<T> &vals, std::vector<int> &idxs) {
    std::vector<T> subset(idxs.size());
    for (int i=0; i<idxs.size(); ++i)
        subset[i]=vals[idxs[i]];
    return subset;
}


template<typename T>
std::vector<T> remove_subset(std::vector<T> vals, std::vector<int> &idxs) {
    std::sort(idxs.begin(), idxs.end(), std::greater<T>());
    for (int idx:idxs)
        vals.erase(vals.begin()+idx);
    return vals;
}

#endif