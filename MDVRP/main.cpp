#include "ga.h"
#include "individual.h"
#include "problem.h"
#include "file.h"

// TODO: Finish implementing neccessary methods for simulation.
// TODO: Check that all edge cases is handled. This should concern the Individual for the most part.
// TODO: Create "integration tests" on the different problems and unit tests.
// TODO: Create namespace for different acceptance criterias.
// TODO: Setup github/workflow.
// TODO: Hyperparameter search.
// TODO: Multithreading of different problems and same problem.
// TODO: Implement different kinds of initialization of individuals.
// Ideas:
//  - initialize customer to closest depot (obs. be carefull of edgecase where it is not illegal to have all at for one depot or similar)
//  - elitism
//  - variation operators with different kinds of "jumps" - small perfections and larger jumps
//  - niching


int main() {
    /**
     *   std::cout << "Hello World!" << std::endl;
     *   std::string problem = "p11";
     *   std::string file_name = "../Data files project 2/Testing Data/Data Files/"+problem;
     *   std::vector<std::vector<double>> problem_file = file::read_flat(file_name);
     *   std::cout << file_name << std::endl;
     */


    std::string problem = "p03";
    std::string file_name = "../Data files project 2/Testing Data/Data Files/"+problem;
    Problem curr_prob = file::load_problem(file_name);
    Individual ind(curr_prob);
    file::write_solution(ind, problem+".res");
    system("pause");
    GA ga(curr_prob, 400);
    ga.simulate(5, 0.1, 0.1, 0.1, 0.1, 10000);
}
