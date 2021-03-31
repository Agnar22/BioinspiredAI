#include "file.h"
#include "individual.h"
#include <iostream>
#include <chrono>

int main() {
    std::cout << "Loading image" << std::endl;
    auto img = file::read_image_to_vector("C:\\Users\\Agnar\\OneDrive - NTNU\\8. semester, 2021\\Bioinspirert AI\\BioinspiredAI\\MOEA\\training_images\\118035\\Test image.jpg");
    std::cout << "Loaded image" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    Individual ind(img);
    auto stop = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(stop-start).count() << std::endl;
    std::cout << "Finished" << std::endl;
}