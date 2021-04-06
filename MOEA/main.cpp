#include "file.h"
#include "individual.h"
#include "objective.h"
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
            q.push(std::make_pair(0, cv::Vec3b{(unsigned char)(rand()%256), (unsigned char)(rand()%256), (unsigned char)(rand()%256)}));
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

int main() {
    std::cout << "Loading image" << std::endl;
    auto img = file::read_image_to_vector("/home/agnar/Git/BioinspiredAI/MOEA/training_images/353013/test_image.jpg");

    std::cout << "Loaded image" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    Individual ind(img);
    auto stop = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(stop-start).count() << std::endl;
    std::cout << "Finished" << std::endl;
    segment_and_display_image(img, ind.genes, img.cols, img.rows, 20);
    std::cout << edge_value(img, ind.genes, ind.root, img.rows, img.cols) << std::endl;
    std::cout << connectivity(ind.root, img.rows, img.cols) << std::endl;
    std::cout << overall_deviation(img, ind.root, img.rows, img.cols) << std::endl;
}