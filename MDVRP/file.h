#ifndef FILE_H
#define FILE_H

#include <vector>
#include <string>
#include <fstream>
#include <regex>
#include "individual.h"
#include "problem.h"

namespace file {
    std::vector<std::vector<double>> read_flat(std::string);
    Problem load_problem(std::string);
    void write_solution(Individual&, std::string);
}

#endif