#ifndef PROBLEM_H
#define PROBLEM_H

#include <vector>
#include <iostream>
#include <math.h>


class Problem {
    private:
        int num_customers, num_depots, vhcl_pr_depot;
        std::vector<double> max_length, max_load;
        std::vector<std::pair<double, double>> positions;
        std::vector<double> cust_serv_dur, cust_demand;
        std::vector<std::vector<double>> distances;
        std::vector<std::vector<double>> depot_distances;

        void calculate_distances();
        void find_depot_distances();
        double euclidian_dist(std::pair<double, double> p1, std::pair<double, double> p2);

    public:
        Problem(int,int,int,std::vector<double>,std::vector<double>,std::vector<std::pair<double, double>>,std::vector<double>,std::vector<double>);
        double get_distance(int, int);
        std::vector<std::vector<double>> get_distances();
        std::vector<std::vector<double>> get_depot_distances();
        int get_num_depots();
        int get_num_customers();
        int get_vhcl_pr_depot();
        double get_customer_load(int);
        std::vector<double> get_max_length();
};

#endif