#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <regex>

// TODO: representation of individuals.
// Ideas:
//  - initialize customer to closest depot (obs. be carefull of edgecase where it is not illegal to have all at for one depot or similar)
//  - elitism
//  - variation operators with different kinds of "jumps" - small perfections and larger jumps
//  - niching


struct Individual {

};

class GA {

};

class Problem {

};

namespace file {
    std::vector<std::vector<double>> read_flat(std::string file_name) {
        std::ifstream fs(file_name);
        std::vector<std::vector<double>> file;
        std::string line;
        std::cout << "Opening file" << std::endl;
        while (std::getline(fs, line)) {
            // Remove leading and trailing spaces.
            line = std::regex_replace(line, std::regex("^ +| +$|( ) +"), "$1");
            std::istringstream s(line);
            std::string val;
            std::vector<double> row;
            while (std::getline(s, val, ' '))
                if (val.compare(" ")!=0)
                    row.push_back(std::stod(val));
            file.push_back(row);
        }
        fs.close();
        return file;
    }
}



int main() {
    std::cout << "Hello World!" << std::endl;
    std::string file_name = "Data files project 2/Testing Data/Data Files/p01";
    std::vector<std::vector<double>> problem_file = file::read_flat(file_name);
    std::cout << file_name << std::endl;
    for (auto line:problem_file) {
        for (auto val:line)
            std::cout << val << " ";
        std::cout << std::endl;
    }
}
