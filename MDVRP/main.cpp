#include "config.h"
#include "ga.h"
#include "individual.h"
#include "problem.h"
#include "file.h"

// TODO: Implement merge route.
// TODO: Implement split route.
// TODO: Implement intra depot mutation.
// TODO: Create test to check all problemts.
// TODO: Multithreading of different problems and same problem.
// TODO: Setup github/workflow.

int main() {
    Problem curr_prob = file::load_problem(PROBLEM_DIR+PROBLEM_NAME);
    GA ga(curr_prob, POPULATION_SIZE, PROBLEM_NAME+".res");
    ga.simulate(TOURNAMENT_SIZE, STOCH_TOURNAMENT_PROB, REV_MUT_PROB, RE_ROUTING_PROB, SWAPPING_PROB, INTER_DEPOT_SWAPPING_INTERVAL);
}
