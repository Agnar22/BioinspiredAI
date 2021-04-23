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
    visited.clear();
    visited.resize(genes.size(), false);
    return root_search(gene, -1);
}

int Individual::root_search(int gene, int max_visited) {
    if (root[gene]!=-1) // FIXME: Is this correct? Shouldn't it be gene?
        return root[gene];
    if (visited[gene])
        return root[gene] = max_visited;
    visited[gene] = true;
    if (get_actual_dir(genes[gene], gene, width, height)==Dir::s) {
        root[gene]=gene;
        return gene;
    }
    int parent_gene = find_pos(gene, genes[gene], width, height);
    // Two neighbours are pointing at each other. Choosing the one with the highest
    // gene as root.
    if (reverse_dir(genes[parent_gene])==genes[gene] && gene>parent_gene) {
        return root[gene] = gene;
    }
    return root[gene]=root_search(parent_gene, std::max(gene, max_visited));

}

void Individual::set_descendants(int node, int val) {
    std::queue<int> bfs;
    bfs.push(node);

    while (!bfs.empty()) {
        int curr_node = bfs.front();
        bfs.pop();
        root[curr_node] = val;

        for (int child:find_children(curr_node)) {
            bfs.push(child);
        }
    }
}

std::vector<int> Individual::find_children(int pos) {
    std::vector<int> children;
    // Up from pos.
    if (pos/width!=0 && genes[pos-width] == Dir::d)
        children.push_back(pos-width);
    // Down from pos.
    if (pos/width!=height-1 && genes[pos+width] == Dir::u)
        children.push_back(pos+width);
    // Left from pos.
    if (pos%width!=0 && genes[pos-1] == Dir::r)
        children.push_back(pos-1);
    // Rifht from pos.
    if (pos%width!=width-1 && genes[pos+1] == Dir::l)
        children.push_back(pos+1);
    return children;
}
}