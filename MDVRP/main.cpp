#include "config.h"
#include "ga.h"
#include "individual.h"
#include "problem.h"
#include "file.h"

// TODO: Implement merge route.
// TODO: Implement split route.
// TODO: Setup github/workflow.
// TODO: Hyperparameter search.
// TODO: Multithreading of different problems and same problem.
// Ideas:
// x- initialize customer to closest depot (obs. be carefull of edgecase where it is not illegal to have all at for one depot or similar)
// x- elitism
// x- variation operators with different kinds of "jumps" - small perfections and larger jumps

int main() {
    Problem curr_prob = file::load_problem(PROBLEM_DIR+PROBLEM_NAME);
    GA ga(curr_prob, POPULATION_SIZE, PROBLEM_NAME+".res");
    ga.simulate(TOURNAMENT_SIZE, STOCH_TOURNAMENT_PROB, REV_MUT_PROB, RE_ROUTING_PROB, SWAPPING_PROB, INTER_DEPOT_SWAPPING_INTERVAL);
}
