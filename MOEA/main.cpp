#include "file.h"
#include "individual.h"
#include "objective.h"
#include "nsga.h"
#include "ga.h"
#include <iostream>
#include <chrono>
#include <queue>
#include <stdlib.h>
#include <string>


void segment_and_display_image(cv::Mat orig_image, std::vector<Dir> genes, int width, int height, int treshold) {
    cv::Mat img(height, width, CV_8UC3, cv::Scalar(0, 0, 0));
    std::vector<bool> visited(genes.size(), false);
    std::queue<std::pair<int, cv::Vec3b>> q;

    for (int gene_pos = 0; gene_pos<genes.size(); ++gene_pos) {
        if (!visited[gene_pos])
            q.push(std::make_pair(gene_pos, cv::Vec3b{(unsigned char)(rand()%256), (unsigned char)(rand()%256), (unsigned char)(rand()%256)}));
        while (!q.empty()) {
            auto top = q.front();
            visited[top.first]=true;
            q.pop();

            int x = top.first%width;
            int y = top.first/width;
            img.at<cv::Vec3b>(y, x)=cv::Vec3b{top.second};

            int parent_pos = find_pos(top.first, genes[top.first], width, height);
            if (!visited[parent_pos]) {
                cv::Vec3b c = top.second;
                if (euc_dist(orig_image.at<cv::Vec3b>(y, x), orig_image.at<cv::Vec3b>(parent_pos/width, parent_pos%width))>treshold)
                    c = cv::Vec3b{(unsigned char)(rand()%256), (unsigned char)(rand()%256), (unsigned char)(rand()%256)};
                q.push(std::make_pair(parent_pos, c));
                visited[parent_pos] = true;
            }
            std::vector<Dir> directions{Dir::u, Dir::r, Dir::d, Dir::l};
            for (Dir dir:directions) {
                try {
                    int child_pos = find_pos(top.first, dir, width, height);
                    if (find_pos(child_pos, genes[child_pos], width, height)==top.first && !visited[child_pos]) {
                        cv::Vec3b c = top.second;
                        if (euc_dist(orig_image.at<cv::Vec3b>(y, x), orig_image.at<cv::Vec3b>(child_pos/width, child_pos%width))>treshold)
                            c = cv::Vec3b{(unsigned char)(rand()%256), (unsigned char)(rand()%256), (unsigned char)(rand()%256)};
                        q.push(std::make_pair(child_pos, c));
                        visited[child_pos] = true;
                    }
                } catch (std::exception e) {
                    std::cout << e.what() << std::endl;
                }
            }
        }
    }
    for (int x=0; x<5; ++x) {
        for (int y=0; y<5; ++y)
            std::cout<<img.at<cv::Vec3b>(y, x) << " ";
        std::cout<<std::endl;
    }
    std::string name="Image";
    cv::namedWindow(name);
    cv::imshow(name, img);
    cv::waitKey(0);
}

bool is_edge(std::vector<int> roots, int x, int y, int width, int height) {
    int pos=y*width+x;
    if (
        (y>0 && x>0 && roots[pos] != roots[pos-width-1]) || // Up-left
        (y>0 && roots[pos] != roots[pos-width]) || // Up
        (y>0 && x<width-1 && roots[pos] != roots[pos-width+1]) || // Up-right
        (x<width-1 && roots[pos] != roots[pos+1]) || // Right
        (y<height-1 && x<width-1 && roots[pos] != roots[pos+width+1]) || // Down-right
        (y<height-1 && roots[pos] != roots[pos+width]) || // Down
        (y<height-1 && x>0 && roots[pos] != roots[pos+width-1]) || // Down-left
        (x>0 && roots[pos] != roots[pos+width-1])// Left
    ) {
        return true;
    }
    return false;
}

std::vector<std::vector<int>> create_type_2_seg(std::vector<int> roots, int width, int height) {
    std::vector<std::vector<int>> type_2_seg;
    for (int y=0; y<height; ++y) {
        std::vector<int> row;
        for (int x=0; x<width; ++x) {
            if (is_edge(roots, x, y, width, height))
                row.push_back(0);
            else
                row.push_back(255);
        }
        type_2_seg.push_back(row);
    }
    return type_2_seg;
}

void display_2d_vector(std::vector<std::vector<int>> borders, cv::Mat img, bool type_1, std::string name) {
    if (type_1) { // Type 1 segmentation.
        for (int y=0; y<borders.size(); ++y)
            for (int x=0; x<borders[y].size(); ++x)
                if (borders[y][x]==0)
                    img.at<cv::Vec3b>(y, x)=cv::Vec3b{0, 255, 0};
        cv::imshow("test", img);
        cv::waitKey(0);

    } else { // Type 2 segmentation.
        cv::Mat img(borders.size(), borders[0].size(), CV_8UC3, cv::Scalar(0, 0, 0));
        for (int y=0; y<borders.size(); ++y)
            for (int x=0; x<borders[y].size(); ++x)
                img.at<cv::Vec3b>(y, x)=cv::Vec3b{(unsigned char)borders[y][x], (unsigned char)borders[y][x], (unsigned char)borders[y][x]};
        cv::imwrite(name, img);
        std::cout << "Saved image: " << name << std::endl;
    }
}

int main() {
    std::cout << "Loading image" << std::endl;
    auto img = file::read_image_to_vector("/home/agnar/Git/BioinspiredAI/MOEA/training_images/353013/test_image.jpg");

    std::cout << "Loaded image" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    Individual ind(img);
    auto stop = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(stop-start).count() << std::endl;
    std::cout << "Finished" << std::endl;
    segment_and_display_image(img, ind.genes, img.cols, img.rows, 2000);
    std::cout << "Edge value: " << obj::edge_value(img, ind.genes, ind.root, img.cols, img.rows) << std::endl;
    std::cout << "Connectivity: " << obj::connectivity(ind.root, img.cols, img.rows) << std::endl;
    std::cout << "Overall deviation: " << obj::overall_deviation(img, ind.root, img.cols, img.rows) << std::endl;
    std::vector<Individual> population;
    for (int ind=0; ind<30; ++ind) {
        population.push_back(Individual(img));
    }
    auto sorted = nsga::fast_nondominated_sort(population, img);
    for (int front=0; front<sorted.size(); ++front) {
        std::cout << "Front " << front << std::endl;
        for (int ind=0; ind<sorted[front].size(); ++ind) {
            std::cout << sorted[front][ind].edge_value << " " << sorted[front][ind].connectivity << " "<< sorted[front][ind].overall_deviation << std::endl;
        }
    }

    std::cout << std::endl;
    auto crowding_sorted = nsga::crowding_sort(sorted[0]);
    for (int ind=0; ind<crowding_sorted.size(); ++ind) {
        std::cout << crowding_sorted[ind].edge_value << " " << crowding_sorted[ind].connectivity << " "<< crowding_sorted[ind].overall_deviation << std::endl;
    }

    GA nsga_ii(30, true, img);
    nsga_ii.simulate();

    for (int num=0; num<nsga_ii.population.size(); ++num) {
        auto type_2_seg = create_type_2_seg(nsga_ii.population[num].root, img.cols, img.rows);
        display_2d_vector(type_2_seg, img, false, "../img/Student_Segmentation_Files/"+std::to_string(num)+".jpg");
    }
    
}