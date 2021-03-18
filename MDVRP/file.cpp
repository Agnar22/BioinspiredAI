#include "file.h"

std::vector<std::vector<double>> file::read_flat(std::string file_name) {
    std::ifstream fs(file_name);
    std::vector<std::vector<double>> file;
    std::string line;
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

Problem file::load_problem(std::string file_name) {
    std::vector<std::vector<double>> prob = read_flat(file_name);

    int max_vhcl_pr_depot = (int)prob[0][0];
    int num_cust = (int)prob[0][1];
    int num_depots = (int)prob[0][2];

    std::vector<double> max_length;
    std::vector<double> max_load;
    for (int row = 1; row < num_depots+1; ++row) {
        max_length.push_back(prob[row][0]);
        max_load.push_back(prob[row][1]);
    }

    std::vector<std::pair<double, double>> positions;
    std::vector<double> cust_serv_dur, cust_demand;
    for (int row = num_depots+1; row < num_cust+num_depots+1; ++row) {
        positions.push_back(std::make_pair(prob[row][1], prob[row][2]));
        cust_serv_dur.push_back(prob[row][3]);
        cust_demand.push_back(prob[row][4]);
    }
    for (int row = num_cust+num_depots+1; row < num_cust+2*num_depots+1; ++row) {
        positions.push_back(std::make_pair(prob[row][1], prob[row][2]));
    }
    return Problem(num_cust, num_depots, max_vhcl_pr_depot, max_length, max_load, positions, cust_serv_dur, cust_demand);
}

void file::write_solution(Individual &ind, std::string file_name) {
    std::ofstream sol_file;
    sol_file.open(file_name);
    double tot_dist=0;
    for (auto depot:ind.trip_dists)
        for (auto trip:depot)
            tot_dist+=trip;
    sol_file << tot_dist << std::endl;
    for (int depot=0; depot<ind.chromosome_trips.size(); ++depot) {
        for (int trip=0; trip<ind.chromosome_trips[depot].size(); ++trip) {
            sol_file << (depot+1) << " " << (trip+1) << " " << ind.trip_dists[depot][trip] << " " << ind.trip_loads[depot][trip] << " " << (depot+1);
            for (int idx=0; idx<ind.chromosome_trips[depot][trip].size(); ++idx) {
                sol_file << " " << (ind.chromosome_trips[depot][trip][idx]+1);
            }
            sol_file << std::endl;
        }
    }
    sol_file.close();
}