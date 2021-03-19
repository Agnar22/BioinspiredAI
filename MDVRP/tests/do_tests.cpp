#include <cmath>
#include <utility>
#include <numeric>
#include "gtest/gtest.h"
#include "../file.h"
#include "../individual.h"
#include "../ga.h"
#include "../problem.h"
#include "../file.h"

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

        static void test_individual(Individual ind, Problem pr, bool test_cust_on_depot) {
            // Check that cust_on_depot is correct.
            std::vector<int> cust_on_depot_f = dim_flat(ind.cust_on_depots);
            if (test_cust_on_depot) {
                EXPECT_EQ(cust_on_depot_f.size(), pr.get_num_customers());
                for (int cust=0; cust<pr.get_num_customers(); ++cust)
                    EXPECT_NE(std::find(cust_on_depot_f.begin(), cust_on_depot_f.end(), cust), cust_on_depot_f.end());
            }

            // Check that chromosome_trips is correct.
            if (test_cust_on_depot)
                EXPECT_EQ(dim_flat(ind.chromosome_trips).size(), pr.get_num_customers());
            for (int depot=0; depot<pr.get_num_depots(); ++depot) {
                EXPECT_LE(ind.chromosome_trips[depot].size(), pr.get_vhcl_pr_depot());
                EXPECT_EQ(dim_flat(ind.chromosome_trips[depot]).size(), ind.cust_on_depots[depot].size());
                for (int cust:ind.cust_on_depots[depot]) {
                    std::vector<int> chr_trip_dep_f = dim_flat(ind.chromosome_trips[depot]);
                    EXPECT_NE(std::find(chr_trip_dep_f.begin(), chr_trip_dep_f.end(), cust), chr_trip_dep_f.end());
                }
            }

            // Check that trip_dists is correct.
            for (int depot=0; depot<pr.get_num_depots(); ++depot) {
                for (int trip=0; trip<ind.chromosome_trips[depot].size(); ++trip) {
                    double correct_trip_dist = Individual::calculate_trip_distance(ind.chromosome_trips[depot][trip], depot, pr);
                    std::cout << depot << " " << trip << " " << ind.trip_dists[depot][trip] << " " << correct_trip_dist << std::endl;
                    EXPECT_LT(std::abs(ind.trip_dists[depot][trip]-correct_trip_dist), 1e-8);
                    if (pr.get_max_length(depot)!=0)
                        EXPECT_LE(ind.trip_dists[depot][trip], pr.get_max_length(depot));
                }
            }

            // Check that trip_loads is correct.
            for (int depot=0; depot<pr.get_num_depots(); ++depot) {
                for (int trip=0; trip<ind.chromosome_trips[depot].size(); ++trip) {
                    double tot_trip_load = 0;
                    for (int pos=0; pos<ind.chromosome_trips[depot][trip].size(); ++pos) {
                        tot_trip_load += pr.get_customer_load(ind.chromosome_trips[depot][trip][pos]);
                    }
                    EXPECT_EQ(ind.trip_loads[depot][trip], tot_trip_load);
                    EXPECT_LE(ind.trip_loads[depot][trip], pr.get_max_load(depot));
                }
            }

            // Check that tot_dist is correct.
            double sum_trip_dists = 0;
            for (std::vector<double> depot_trip_dists:ind.trip_dists)
                sum_trip_dists += std::accumulate(depot_trip_dists.begin(), depot_trip_dists.end(), 0.0);
            EXPECT_LT(std::abs(ind.tot_dist-sum_trip_dists), 1e-8);
        }
};

struct TestGA: public testing::Test {

    public:
        Problem pr;
        GA ga;
        int num_individuals;

        void SetUp() {
            srand(42);
            std::string problem = "p06";
            std::string file_name = "../../Data files project 2/Testing Data/Data Files/"+problem;
            pr = file::load_problem(file_name);
            num_individuals = 10;
            ga = GA(pr, num_individuals);
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

TEST_F(TestIndividual, initializer) {
    TestIndividual::test_individual(ind, pr, true);
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

TEST_F(TestIndividual, remove_customers) {
    try {
        std::vector<int> cust_to_remove{pr.get_num_customers()};
        ind.remove_customers(cust_to_remove, pr);
        FAIL();
    } catch (std::runtime_error r) {
        SUCCEED();
    }
    for (int cust=0; cust<pr.get_num_customers(); ++cust) {
        // TODO: Some of this code can be removed as we are using TestIndividual::test_individual(ind, pr)
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
        TestIndividual::test_individual(ind, pr, false);
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

TEST_F(TestIndividual, marginal_cost) {
    EXPECT_LT(std::abs(Individual::marginal_cost(10, 5, 15, pr)-19.73497714), 1e-8);
    EXPECT_LT(std::abs(Individual::marginal_cost(51, 15, 10, pr)-0.01887903), 1e-8);
}

TEST(Individual, insert_customer) {}

TEST_F(TestIndividual, get_subset) {
    std::vector<double> values{10, 20, 30, 25, 15, 7};
    std::vector<double> subset_of_values{20, 30, 15, 7};
    std::vector<int> idxs{1, 2, 4, 5};

    EXPECT_EQ(get_subset(values, idxs), subset_of_values);
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

TEST_F(TestGA, best_cost_route_crossover) {
    Individual l = ga.get_individual(1);
    Individual r = ga.get_individual(8);
    auto parents = std::make_pair(l, r);

    auto children = GA::best_cost_route_crossover(parents, pr);

    TestIndividual::test_individual(l, pr, true);
    TestIndividual::test_individual(r, pr, true);
    TestIndividual::test_individual(children.first, pr, true);
    TestIndividual::test_individual(children.second, pr, true);
}


TEST_F(TestGA, mutate) {
    for (int pos=0; pos<3; ++pos) {
        Individual ind = ga.get_individual(8);
        Individual ind_cp = ind;
        int attempts = 0;
        while (ind==ind_cp && attempts<30) {
            ga.mutate(ind, 1.0*(pos==0), 1.0*(pos==1), 1.0*(pos==2), 0, pr);
            ++attempts;
        }
        TestIndividual::test_individual(ind, pr, true);
        EXPECT_FALSE(ind==ind_cp);
    }
}

TEST_F(TestGA, initialized_population) {
    for (int ind=0; ind<num_individuals; ++ind) {
        TestIndividual::test_individual(ga.get_individual(ind), pr, true);
    }
}

TEST_F(TestGA, simulate) {
    num_individuals = 100;
    ga = GA(pr, num_individuals);
    // TODO: Add test to check that the population/the best improves.
    ga.simulate(5, 0.5, 1.0, 1.0, 1.0, 10000);

    for (int ind=0; ind<num_individuals; ++ind) {
        TestIndividual::test_individual(ga.get_individual(ind), pr, true);
    }
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
