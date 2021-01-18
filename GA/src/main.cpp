#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include <fstream>

std::vector<std::vector<double>> read_csv(std::string file_name) {
    std::ifstream fs(file_name);
    std::vector<std::vector<double>> file;
    std::string line;
    while (std::getline(fs, line)) {
        std::istringstream s(line);
        std::string val;
        std::vector<double> row;
        while (std::getline(s, val, ','))
            row.push_back(std::stod(val));
        file.push_back(row);
    }
    fs.close();
    return file;
}

int main() {
   auto data = read_csv("data.csv"); 
   std::cout << "read file " << data.size() << " " << data[0].size() << std::endl; 
}
