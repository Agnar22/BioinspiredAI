#ifndef OBJECTIVE_H
#define OBJECTIVE_h

#include <vector>
#include <opencv2/opencv.hpp>
#include <map>
#include "dir.h"
#include "individual.h"

static float edge_value(cv::Mat &orig_img, std::vector<Dir> &genes, std::vector<int> &root, int width, int height) {
    float sum_edge_value=0;
    for (int row=0; row<height; ++row) {
        for (int col=0; col<width; ++col) {
            int curr_gene = row*width+col;
            int neigh_down = curr_gene+width;
            int neigh_right = curr_gene+1;
            if (row<height-1 && root[curr_gene]!=root[neigh_down])
                sum_edge_value += euc_dist(orig_img.at<cv::Vec3b>(row, col), orig_img.at<cv::Vec3b>(row+1, col));
            if (col<width-1 && root[curr_gene]!=root[neigh_right])
                sum_edge_value += euc_dist(orig_img.at<cv::Vec3b>(row, col), orig_img.at<cv::Vec3b>(row, col+1));
        }
    }
    return 2*sum_edge_value;
}

static float connectivity(std::vector<int> &root, int width, int height) {
    float sum_connectivity=0;
    for (int row=0; row<height; ++row) {
        for (int col=0; col<width; ++col) {
            int curr_gene = row*width+col;
            int neigh_down = curr_gene+width;
            int neigh_right = curr_gene+1;
            if (row<height-1 && root[curr_gene]!=root[neigh_down])
                sum_connectivity += 1/8;
            if (col<width-1 && root[curr_gene]!=root[neigh_right])
                sum_connectivity += 1/8;
        }
    }
    return 2*sum_connectivity;
}

static std::map<int, cv::Vec3d> find_segment_means(cv::Mat &img, std::vector<int> &root) {
    std::map<int, cv::Vec3d> segment_sum;
    std::map<int, int> pxls_in_segment;
    for (int pos=0; pos<root.size(); ++pos) {
        int row = pos/img.rows;
        int col = pos%img.rows;
        auto cur_pxl = img.at<cv::Vec3b>(row, col);
        if (segment_sum.find(pos)==segment_sum.end()) {
            segment_sum[pos] = cur_pxl;
            pxls_in_segment[pos] = 1;
        } else {
            segment_sum[pos][0] += cur_pxl[0];
            segment_sum[pos][1] += cur_pxl[1];
            segment_sum[pos][2] += cur_pxl[2];
            pxls_in_segment[pos] += 1;
        }
    }

    std::map<int, cv::Vec3d> segment_means;
    for (std::pair<int, cv::Vec3b> elem:segment_sum) {
        int key = elem.first;
        cv::Vec3b value = elem.second;
        cv::Vec3b segment_sum = value;
        segment_means[key] = cv::Vec3b{
            segment_sum[0]/pxls_in_segment[key],
            segment_sum[1]/pxls_in_segment[key],
            segment_sum[2]/pxls_in_segment[key]
        };
    }
    return segment_means;
}

static float overall_deviation(cv::Mat &orig_img, std::vector<int> &root, int width, int height) {
    float sum_overall_deviation=0;
    std::map<int, cv::Vec3d> segment_means=find_segment_means(orig_img, root);
    for (int row=0; row<height; ++row) {
        for (int col=0; col<width; ++col) {
            int curr_gene = row*width+col;
            sum_overall_deviation += euc_dist(orig_img.at<cv::Vec3b>(row, col), segment_means[root[curr_gene]]);
        }
    }
    return sum_overall_deviation;
}

#endif