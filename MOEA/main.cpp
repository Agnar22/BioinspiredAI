#include "file.h"
#include "individual.h"
#include <iostream>

int main() {
    std::cout << "Loading image" << std::endl;
    auto img = file::read_image_to_vector("C:\\Users\\Agnar\\OneDrive - NTNU\\8. semester, 2021\\Bioinspirert AI\\BioinspiredAI\\MOEA\\training_images\\118035\\Test image.jpg");
    std::cout << "Loaded image" << std::endl;
    Individual ind(img);
    std::cout << "Finished" << std::endl;
}