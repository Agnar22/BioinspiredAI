#include <cmath>
#include <utility>
#include "gtest/gtest.h"
#include "../file.h"
#include "../individual.h"
#include "../problem.h"
#include "../file.h"

struct TestProblem: public testing::Test {
    public:
        Problem pr;

        void SetUp() {
            std::string problem = "p01";
            std::string file_name = "../../Data files project 2/Testing Data/Data Files/"+problem;
            pr = file::load_problem(file_name);
        };
};

struct TestIndividual: public testing::Test {
    public:
        Problem pr;
        Individual ind;

        void SetUp() {
            std::string problem = "p01";
            std::string file_name = "../../Data files project 2/Testing Data/Data Files/"+problem;
            pr = file::load_problem(file_name);
            ind = Individual(pr);
        };
};

std::vector<int> dim_flat (std::vector<int> const & v) { return v; }

template <typename T>
std::vector<int> dim_flat(std::vector<std::vector<T>> const & v) {
    std::vector<int> ret;
    for ( auto const & e : v ) {
        auto s = dim_flat(e);
        ret.reserve( ret.size() + s.size() );
        ret.insert( ret.end(), s.cbegin(), s.cend() );
    }
    return ret;
}

TEST_F(TestProblem, load_problem) {
    EXPECT_EQ(pr.get_num_depots(), 4);
    EXPECT_EQ(pr.get_num_customers(), 50);
    EXPECT_EQ(pr.get_vhcl_pr_depot(), 4);
    EXPECT_EQ(pr.get_customer_load(10), 19);
    EXPECT_EQ(pr.get_max_length(2), 0);
    EXPECT_EQ(pr.get_max_load(2), 80);
    EXPECT_EQ(pr.get_depot_distances().size(), 50);
    EXPECT_EQ(pr.get_depot_distances()[0].size(), 4);
    EXPECT_LT(std::abs(pr.get_depot_distances()[2][2]-34.05877273), 0.001);
    EXPECT_LT(std::abs(pr.get_distances()[2][52]-34.05877273), 0.001);
    EXPECT_LT(std::abs(pr.get_depot_distances()[4][1]-14.14213562), 0.001);
    EXPECT_LT(std::abs(pr.get_distances()[4][51]-14.14213562), 0.001);
    EXPECT_LT(std::abs(pr.get_depot_distances()[49][0]-39.81205847), 0.001);
    EXPECT_LT(std::abs(pr.get_distances()[49][50]-39.81205847), 0.001);
}

TEST_F(TestIndividual, calculate_trip_distance) {
    std::vector<int> customers = {0, 10, 30};
    int depot = 1;
    double trip_dist = Individual::calculate_trip_distance(customers, depot, pr);
    EXPECT_LT(std::abs(trip_dist - 84.25123), 0.001);
}

TEST_F(TestIndividual, get_fitness) {
    std::cout << "get_fitness total_distance: " << ind.tot_dist << std::endl;
    ind.tot_dist=100.0;
    EXPECT_EQ(1.0/100.0, ind.get_fitness());
}

TEST_F(TestIndividual, remove_customers) {
    try {
        std::vector<int> cust_to_remove{pr.get_num_customers()};
        ind.remove_customers(cust_to_remove, pr);
        FAIL();
    } catch (std::runtime_error r) {
        SUCCEED();
    }
    for (int cust=0; cust<pr.get_num_customers(); ++cust) {
        std::vector<int> cust_on_depots_f = dim_flat(ind.cust_on_depots);
        EXPECT_NE(std::find(cust_on_depots_f.begin(), cust_on_depots_f.end(), cust), cust_on_depots_f.end());
        std::vector<int> chromosome_trips_f = dim_flat(ind.chromosome_trips);
        EXPECT_NE(std::find(chromosome_trips_f.begin(), chromosome_trips_f.end(), cust), chromosome_trips_f.end());

        std::vector<int> cust_to_remove{cust};
        ind.remove_customers(cust_to_remove, pr);
        cust_on_depots_f = dim_flat(ind.cust_on_depots);
        EXPECT_EQ(std::find(cust_on_depots_f.begin(), cust_on_depots_f.end(), cust), cust_on_depots_f.end());
        chromosome_trips_f = dim_flat(ind.chromosome_trips);
        EXPECT_EQ(std::find(chromosome_trips_f.begin(), chromosome_trips_f.end(), cust), chromosome_trips_f.end());
        // TODO: Check that the change in: trip_dists, trip_loads, tot_dist is correct.
    }
    for (int depot=0; depot<pr.get_num_depots(); ++depot) {
        EXPECT_EQ(ind.cust_on_depots[depot].size(), 0);
        EXPECT_EQ(ind.chromosome_trips[depot].size(), 0);
        EXPECT_EQ(ind.trip_dists[depot].size(), 0);
        EXPECT_EQ(ind.trip_loads[depot].size(), 0);
    }
    EXPECT_LT(std::abs(ind.tot_dist), 1e-8);
}

TEST(Individual, reversal_mutation) {}

TEST(Individual, re_routing_mutation) {}

TEST(Individual, swapping_mutation) {}

TEST(Individual, find_insert_costs) {}

TEST(Individual, insert_stochastically) {}

TEST(Individual, marginal_cost) {}

TEST(Individual, insert_customer) {}

TEST_F(TestIndividual, get_subset) {
    std::vector<double> values{10, 20, 30, 25, 15, 7};
    std::vector<double> subset_of_values{20, 30, 15, 7};
    std::vector<int> idxs{1, 2, 4, 5};

    EXPECT_EQ(Individual::get_subset(values, idxs), subset_of_values);
}

TEST_F(TestIndividual, comparison) {
    Individual ind2;
    ind.tot_dist=100;
    ind2.tot_dist=200;

    EXPECT_LT(ind2, ind);
    std::swap(ind.tot_dist, ind2.tot_dist);
    EXPECT_LT(ind, ind2);

    std::vector<std::pair<double, double>> illegal_distances = {std::make_pair(0, 100), std::make_pair(20, 0), std::make_pair(100001, 100), std::make_pair(100, 100001)};
    for (std::pair<double, double> pair_dist: illegal_distances) {
        try {
            ind.tot_dist=pair_dist.first;
            ind2.tot_dist=pair_dist.second;
            ind<ind2;
            FAIL();
        } catch (std::runtime_error r) {
            SUCCEED();
        }
    }
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}