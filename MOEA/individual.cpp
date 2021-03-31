#include "individual.h"

Individual::Individual(cv::Mat img) {
    initialize_genes(img);
}

void Individual::initialize_genes(cv::Mat img) {
    genes.clear();
    genes.assign(img.rows*img.cols, Dir::s);

    int src=rand()%(img.rows*img.cols);
    Graph prim(img.rows*img.cols, img.cols, img.rows);
    
    for (int row=0; row<img.rows; ++row) {
        for (int col=0; col<img.cols; ++col) {
            int cur_pos = row*img.cols+col;
            if (row!=0)
                prim.add_edge(cur_pos, cur_pos-img.cols, euc_dist(img.at<cv::Vec3b>(cur_pos),img.at<cv::Vec3b>(cur_pos-img.cols)));
            if (col!=0)
                prim.add_edge(cur_pos, cur_pos-1, euc_dist(img.at<cv::Vec3b>(cur_pos),img.at<cv::Vec3b>(cur_pos-1)));
        }
    }
    genes = prim.prim_mst(src);
}