#include <cmath>
#include <utility>
#include "gtest/gtest.h"
#include "../file.h"

struct TestProblem: public testing::Test {
    public:
        Problem pr;

        virtual void SetUp() {
            std::string problem = "p01";
            std::string file_name = "../../Data files project 2/Testing Data/Data Files/"+problem;
            pr = file::load_problem(file_name);
        };
};

struct TestIndividual: public testing::Test {
    public:
        Problem pr;
        Individual ind;

        virtual void SetUp() {
            std::string problem = "p01";
            std::string file_name = "../../Data files project 2/Testing Data/Data Files/"+problem;
            pr = file::load_problem(file_name);
            ind = Individual(pr);
        };
};

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
    ind.tot_dist=100.0;
    EXPECT_EQ(1.0/100.0, ind.get_fitness());
}

TEST(Individual, remove_customers) {}

TEST(Individual, reversal_mutation) {}

TEST(Individual, re_routing_mutation) {}

TEST(Individual, swapping_mutation) {}

TEST(Individual, find_insert_costs) {}

TEST(Individual, insert_stochastically) {}

TEST(Individual, marginal_cost) {}

TEST(Individual, insert_customer) {}

TEST(Individual, get_subset) {}

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