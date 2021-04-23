#include "objective.h"

double obj::edge_value(cv::Mat &orig_img, std::vector<Dir> &genes, std::vector<int> &root, int width, int height) {
    double sum_edge_value=0;
    for (int row=0; row<height; ++row) {
        for (int col=0; col<width; ++col) {
            int curr_gene = row*width+col;
            int neigh_down = curr_gene+width;
            int neigh_right = curr_gene+1;
            int neigh_down_right = curr_gene+width+1;
            int neigh_down_left = curr_gene+width-1;
            if (row<height-1 && root[curr_gene]!=root[neigh_down]) {
                sum_edge_value += euc_dist(orig_img.at<cv::Vec3b>(row, col), orig_img.at<cv::Vec3b>(row+1, col));
            }
            if (col<width-1 && root[curr_gene]!=root[neigh_right]) {
                sum_edge_value += euc_dist(orig_img.at<cv::Vec3b>(row, col), orig_img.at<cv::Vec3b>(row, col+1));
            }
            if (row<height-1 && col<width-1 && root[curr_gene]!=root[neigh_down_right]) {
                sum_edge_value += euc_dist(orig_img.at<cv::Vec3b>(row, col), orig_img.at<cv::Vec3b>(row+1, col+1));
            }
            if (row<height-1 && col>0 && root[curr_gene]!=root[neigh_down_left]) {
                sum_edge_value += euc_dist(orig_img.at<cv::Vec3b>(row, col), orig_img.at<cv::Vec3b>(row+1, col-1));
            }
        }
    }
    return 2*sum_edge_value;
}

double obj::connectivity(std::vector<int> &root, int width, int height) {
    double sum_connectivity=0.0;
    for (int row=0; row<height; ++row) {
        for (int col=0; col<width; ++col) {
            int curr_gene = row*width+col;
            int neigh_down = curr_gene+width;
            int neigh_right = curr_gene+1;
            int neigh_down_right = curr_gene+width+1;
            int neigh_down_left = curr_gene+width-1;
            if (row<height-1 && root[curr_gene]!=root[neigh_down])
                sum_connectivity += 1.0/8.0;
            if (col<width-1 && root[curr_gene]!=root[neigh_right])
                sum_connectivity += 1.0/8.0;
            if (row<height-1 && col<width-1 && root[curr_gene]!=root[neigh_down_right])
                sum_connectivity += 1.0/8.0;
            if (row<height-1 && col>0 && root[curr_gene]!=root[neigh_down_left])
                sum_connectivity += 1.0/8.0;
        }
    }
    return 2.0*sum_connectivity;
}

std::map<int, std::vector<double>> obj::find_segment_means(cv::Mat &img, std::vector<int> &root) {
    std::map<int, std::vector<double>> segment_sum;
    std::map<int, int> pxls_in_segment;
    for (int pos=0; pos<root.size(); ++pos) {
        int row = pos/img.cols;
        int col = pos%img.cols;
        int curr_root = root[pos];
        auto cur_pxl = img.at<cv::Vec3b>(row, col);
        if (segment_sum.find(curr_root)==segment_sum.end()) {
            segment_sum[curr_root] = {(double)cur_pxl[0], (double)cur_pxl[1], (double)cur_pxl[2]};
            pxls_in_segment[curr_root] = 1;
        } else {
            segment_sum[curr_root][0] += cur_pxl[0];
            segment_sum[curr_root][1] += cur_pxl[1];
            segment_sum[curr_root][2] += cur_pxl[2];
            pxls_in_segment[curr_root]++;
        }
    }

    std::map<int, std::vector<double>> segment_means;
    for (std::pair<int, std::vector<double>> elem:segment_sum) {
        int key = elem.first;
        std::vector<double> segment_sum = elem.second;
        segment_means[key] = std::vector<double>{
            (double)(segment_sum[0])/(double)(pxls_in_segment[key]),
            (double)(segment_sum[1])/(double)(pxls_in_segment[key]),
            (double)(segment_sum[2])/(double)(pxls_in_segment[key])
        };
    }
    return segment_means;
}

double obj::overall_deviation(cv::Mat &orig_img, std::vector<int> &root, int width, int height) {
    double sum_overall_deviation=0;
    std::map<int, std::vector<double>> segment_means=find_segment_means(orig_img, root);
    for (int i=0; i<9; ++i) {
        auto pxl = orig_img.at<cv::Vec3b>(i/width, i%width);
        auto mean = segment_means[root[i]];
    }
    for (int row=0; row<height; ++row) {
        for (int col=0; col<width; ++col) {
            int curr_gene = row*width+col;
            sum_overall_deviation += euc_dist(orig_img.at<cv::Vec3b>(row, col), segment_means[root[curr_gene]]);
        }
    }
    return sum_overall_deviation;
}