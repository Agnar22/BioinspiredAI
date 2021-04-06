#include "individual.h"

Individual::Individual(cv::Mat img) {
    initialize_genes(img);
    find_roots();
}

void Individual::calculate_objectives(cv::Mat &img) {
    // TODO: Can be optimized to only calculate where something has changed.
    edge_value = obj::edge_value(img, genes, root, img.cols, img.rows);
    connectivity = obj::connectivity(root, img.cols, img.rows);
    overall_deviation = obj::overall_deviation(img, root, img.cols, img.rows);
}

void Individual::initialize_genes(cv::Mat img) {
    genes.clear();
    genes.assign(img.rows*img.cols, Dir::s);
    width = img.cols;
    height = img.rows;

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

void Individual::find_roots() {
    // FIXME: Must be carefull of "loops".
    root.clear();
    root.resize(genes.size(), -1);
    for (int gene=0; gene<genes.size(); ++gene)
        find_root(gene);
}

int Individual::find_root(int gene) {
    if (root[gene]!=-1)
        return root[gene];
    if (get_actual_dir(genes[gene], gene, width, height)==Dir::s)
        return root[gene]=gene;
    int parent_gene = find_pos(gene, genes[gene], width, height);
    return root[gene]=find_root(parent_gene);
}