#include <cmath>
#include "gtest/gtest.h"
#include "../file.h"

TEST(File, load_problem) {
    std::string problem = "p01";
    std::string file_name = "../../Data files project 2/Testing Data/Data Files/"+problem;
    Problem pr = file::load_problem(file_name);

    ASSERT_EQ(pr.get_num_depots(), 4);
    ASSERT_EQ(pr.get_num_customers(), 50);
    ASSERT_EQ(pr.get_vhcl_pr_depot(), 4);
    ASSERT_EQ(pr.get_customer_load(10), 19);
    ASSERT_EQ(pr.get_max_length(2), 0);
    ASSERT_EQ(pr.get_max_load(2), 80);
    ASSERT_EQ(pr.get_depot_distances().size(), 50);
    ASSERT_EQ(pr.get_depot_distances()[0].size(), 4);
    ASSERT_LT(std::abs(pr.get_depot_distances()[2][2]-34.05877273), 0.001);
    ASSERT_LT(std::abs(pr.get_distances()[2][52]-34.05877273), 0.001);
    ASSERT_LT(std::abs(pr.get_depot_distances()[4][1]-14.14213562), 0.001);
    ASSERT_LT(std::abs(pr.get_distances()[4][51]-14.14213562), 0.001);
    ASSERT_LT(std::abs(pr.get_depot_distances()[49][0]-39.81205847), 0.001);
    ASSERT_LT(std::abs(pr.get_distances()[49][50]-39.81205847), 0.001);
}

TEST(Individual, calculate_trip_distance) {}

TEST(Individual, get_fitness) {}

TEST(Individual, remove_customers) {}

TEST(Individual, reversal_mutation) {}

TEST(Individual, re_routing_mutation) {}

TEST(Individual, swapping_mutation) {}

TEST(Individual, find_insert_costs) {}

TEST(Individual, insert_stochastically) {}

TEST(Individual, marginal_cost) {}

TEST(Individual, insert_customer) {}

TEST(Individual, get_subset) {}

TEST(Individual, comparison) {}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}