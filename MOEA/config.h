#ifndef CONFIG_H
#define CONFIG_H

#include <string>

static std::string IMAGE_PATH = "/home/agnar/Git/BioinspiredAI/MOEA/training_images/353013/Test image.jpg";
static std::string TYPE_1_SAVE_PATH = "../img/Student_Segmentation_Files_Type_1/";
static std::string TYPE_2_SAVE_PATH = "../img/Student_Segmentation_Files/";

static bool NSGA_ii = true;
static int POPULATION_SIZE = 100;
static int GENERATIONS = 15;
static int CONNECTIVITY_MIN = 50;

static int MIN_SEGMENTS = 4;
static int MAX_SEGMENTS = 20;

static int SPLIT_MUTATION = 0.3;

#endif