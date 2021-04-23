#include "gtest/gtest.h"
#include "../individual.h"
#include "../nsga.h"
#include "../dir.h"
#include <opencv2/opencv.hpp>
#include <iostream>

struct TestImage: public testing::Test {
    public:
        cv::Mat test_image = {3, 3, CV_8UC3, cv::Scalar(0,0,0)};

        void SetUp() {
            std::vector<cv::Vec3b> pxls = {
                cv::Vec3b{0,0,0},
                cv::Vec3b{1,2,3},
                cv::Vec3b{10,20,30},
                cv::Vec3b{100,100,100},
                cv::Vec3b{110,120,130},
                cv::Vec3b{110,120,130},
                cv::Vec3b{200,200,200},
                cv::Vec3b{200,200,200},
                cv::Vec3b{200,200,200},
            };
            for (int row=0; row<3; ++row)
                for (int col=0; col<3; ++col)
                    test_image.at<cv::Vec3b>(row, col)=pxls[row*3+col];
        };
};

TEST_F(TestImage, simple) {
    test_image.at<cv::Vec3b>(0, 0)=cv::Vec3b{1, 1, 1};
    std::cout << test_image << std::endl;
    EXPECT_EQ(true, true);
}

TEST_F(TestImage, init_individual) {
    Individual ind(test_image);
    EXPECT_EQ(ind.genes[0], Dir::r);
    EXPECT_EQ(ind.genes[1], Dir::s);
    EXPECT_EQ(ind.genes[2], Dir::s);
    EXPECT_EQ(ind.genes[3], Dir::s);
    EXPECT_EQ(ind.genes[4], Dir::s);
    EXPECT_EQ(ind.genes[5], Dir::l);
    EXPECT_EQ(ind.genes[6], Dir::r);
    EXPECT_EQ(ind.genes[7], Dir::s);
    EXPECT_EQ(ind.genes[8], Dir::l);
}

TEST_F(TestImage, objectives) {
    Individual ind(test_image);
    ind.calculate_objectives(test_image);

    EXPECT_DOUBLE_EQ(ind.edge_value, 4842.3314038349436); // Value after calculating was 2*2421.168415
    EXPECT_DOUBLE_EQ(ind.connectivity, 32.0/8.0);
    EXPECT_DOUBLE_EQ(ind.overall_deviation, 3.7416573867739413);
}

TEST_F(TestImage, roots) {
    Individual ind(test_image);
    EXPECT_EQ(ind.root[0], 0);
    EXPECT_EQ(ind.root[1], 0);
    EXPECT_EQ(ind.root[2], 2);
    EXPECT_EQ(ind.root[3], 3);
    EXPECT_EQ(ind.root[4], 4);
    EXPECT_EQ(ind.root[5], 4);
    EXPECT_EQ(ind.root[6], 7);
    EXPECT_EQ(ind.root[7], 7);
    EXPECT_EQ(ind.root[8], 7);
}

TEST_F(TestImage, fast_nondominated_sort) {
    std::vector<Individual> inds = {
        Individual(test_image), Individual(test_image), Individual(test_image), Individual(test_image),
        Individual(test_image), Individual(test_image), Individual(test_image), Individual(test_image)
    };

    std::vector<std::vector<double>> objectives = {
        {9, 1.2, 2}, {9.5, 1.3, 1.2}, {8.5, 2.1, 0.92}, // Second front.
        {10, 1, 1}, {9, 2, 0.9}, // First front.
        {0, 10, 10}, // Last front.
        {7, 3, 3}, {6, 2, 5} // Third front.
    };

    for (int pos=0; pos<8; ++pos) {
        inds[pos].edge_value = objectives[pos][0];
        inds[pos].connectivity = objectives[pos][1];
        inds[pos].overall_deviation = objectives[pos][2];
    }

    std::vector<int> answers = {1, 1, 1, 0, 0, 3, 2, 2};

    auto sorted = nsga::fast_nondominated_sort(inds, test_image, false);

    EXPECT_EQ(sorted.size(), 4);

    for (int front=0; front<4; ++front) {
        for (int ind=0; ind<sorted[front].size(); ++ind){
            int pos = std::find(inds.begin(), inds.end(), sorted[front][ind])-inds.begin();
            EXPECT_EQ(answers[pos], front);
        }
    }
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}