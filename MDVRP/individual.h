#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

#include <algorithm>
#include <vector>
#include <utility>
#include "problem.h"

class Individual {
    public:
        std::vector<std::vector<int>> chromosomes;
        std::vector<std::vector<int>> cust_on_depots;
        std::vector<std::vector<std::vector<int>>> chromosome_trips;
        std::vector<std::vector<double>> trip_dists;
        std::vector<std::vector<int>> trip_loads;
        double tot_dist;

        Individual(Problem&);
        void initialize_chromosomes(Problem&);
        std::vector<std::vector<int>> assign_customers_to_depots(Problem&, bool);
        void setup_trips(int, std::vector<int>, Problem&);
        void setup_trips_forward(int, std::vector<int>, Problem&);
        void setup_trips_backward(int, std::vector<int>, Problem&);
        static double calculate_trip_distance(std::vector<int>&, int, Problem&);
        double get_fitness();
        void remove_customers(std::vector<int>&, Problem&);
        int remove_from_2d_vector(std::vector<std::vector<int>>&, int);
        void reversal_mutation(int, Problem&);
        void re_routing_mutation(int, Problem&);
        void swapping_mutation(int, Problem&);
        std::pair<std::vector<double>, std::vector<std::pair<int, int>>> find_insert_costs(int, int, Problem&);
        void insert_stochastically(int, double, int, Problem&);
        double marginal_cost(int, int, int, Problem&);
        void insert_customer(int, int, int, int, Problem&);
        template <class T>
        static std::vector<T> get_subset(std::vector<T>&, std::vector<int>&);
};

inline bool operator<(const Individual &a, const Individual&b) {
    if (a.tot_dist<1 || a.tot_dist>100000)
        throw std::runtime_error("Individual fitness is not within plausible range.");
    if (b.tot_dist<1 || b.tot_dist>100000)
        throw std::runtime_error("Individual fitness is not within plausible range.");
    return a.get_fitness()<b.get_fitness();
};

#endif