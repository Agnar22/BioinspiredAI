#ifndef CONFIG_H
#define CONFIG_H

#include <string>

/**
 *  GA
 */
static int POPULATION_SIZE = 300;
static int NUM_GENERATIONS = 300;
static int TOURNAMENT_SIZE = 6;
static double ELITE_PERCENTAGE = 3.0; // Try to tune
static double STOCH_TOURNAMENT_PROB = 0.2; // Might try to tune
static double REV_MUT_PROB = 0.3; // Try to rune
static double RE_ROUTING_PROB = 0.9; // Try to tune
static double SWAPPING_PROB = 0.5; // Try to tune
static int INTER_DEPOT_SWAPPING_INTERVAL = 30;
static double INTER_DEPOT_RANDOM_BEST_PROB = 0.2;

static int BEST_COST_ROUTE_TRIES = 3;
static double GREEDY_INSERT_PROB = 0.95; // Try to tune

/**
 *  Individual
 */
// Assigning customers to depot.
static bool STOCHASTIC_ASSIGNMENT = true;
static int POW_INV_DEPOT_DIST = 10;
static int POW_RE_ASSIGN_INV_DEPOT_DIST = 7;

// Ordering of customers within depots.
static double GREEDY_PROB = 0.1; // Might try to tune
static int RETRY_ATTEMPTS = 3;
// - Setup trips forward.
static int POW_INV_CUST_DIST = 8;
// - c_and_w_algorithm.
static double DROP_PROB = 0.03; // Try to tune

/**
 * Other
 */
static std::string PROBLEM_DIR = "../Data files project 2/Testing Data/Data Files/";
static std::string PROBLEM_NAME = "p01";

#endif