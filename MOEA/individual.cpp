#include "individual.h"

Individual::Individual(cv::Mat img, int treshold) {
    initialize_genes(img, treshold);
    find_roots();
}

Individual::Individual(Individual &l, Individual &r, std::vector<int> &segments, int h, int w) {
    width = w;
    height = h;
    genes.clear();
    for (int pos=0; pos<width*height; ++pos) {
        if (segments[pos]==0)
            genes.push_back(l.genes[pos]);
        else
            genes.push_back(r.genes[pos]);
    }
    //genes.insert(genes.end(), l.genes.begin(), l.genes.begin()+crossover_pos);
    //genes.insert(genes.end(), r.genes.begin()+crossover_pos, r.genes.end());
    //std::cout << "finding roots" << std::endl;
    find_roots();
    //std::cout << "found roots" << std::endl;
}

Individual::Individual(Individual &l, Individual &r, int crossover_pos, cv::Mat &img) {
    width = img.cols;
    height = img.rows;
    genes.clear();
    genes.insert(genes.end(), l.genes.begin(), l.genes.begin()+crossover_pos);
    genes.insert(genes.end(), r.genes.begin()+crossover_pos, r.genes.end());
    //std::cout << "finding roots" << std::endl;
    find_roots();
    //std::cout << "found roots" << std::endl;
}

void Individual::calculate_objectives(cv::Mat &img) {
    // TODO: Can be optimized to only calculate where something has changed.
    edge_value = obj::edge_value(img, genes, root, img.cols, img.rows);
    connectivity = obj::connectivity(root, img.cols, img.rows);
    overall_deviation = obj::overall_deviation(img, root, img.cols, img.rows);
}

void Individual::initialize_genes(cv::Mat img, int segments) {
    genes.clear();
    genes.assign(img.rows*img.cols, Dir::s);
    width = img.cols;
    height = img.rows;


    UnionFind uf(img.rows*img.cols);
    auto rg = uf.kruskal_prim_mst(segments, img);
    root = rg.first;
    genes = rg.second;
    /*
    int src=rand()%(img.rows*img.cols);
    //std::cout << "src: " << src << std::endl;
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
    genes = prim.prim_mst(src, treshold);
    */
}

void Individual::find_roots() {
    // FIXME: Must be carefull of "loops".
    root.clear();
    root.resize(genes.size(), -1);
    for (int gene=0; gene<genes.size(); ++gene) {
        //std::cout << "gene: " << gene << " genes.size() " << genes.size() << " root.size() " << root.size() << std::endl;
        // Crashing at gene=188.
        find_root(gene);
        //std::cout << "found root for " << gene << std::endl;
    }
}

int Individual::find_root(int gene) {
    visited.clear();
    visited.resize(genes.size(), false);
    return root_search(gene, -1);
}

int Individual::root_search(int gene, int max_visited) {
    //if (print)
    //std::cout << "gene: " << gene << std::endl;
    if (root[gene]!=-1) // FIXME: Is this correct? Shouldn't it be gene?
        return root[gene];
    //std::cout << "visited: " << gene << std::endl;
    if (visited[gene])
        return root[gene] = max_visited;
    //std::cout << "get_actual " << gene << std::endl;
    visited[gene] = true;
    if (get_actual_dir(genes[gene], gene, width, height)==Dir::s) {
        //std::cout << "pointing to itself " << gene << " " << root.size() << std::endl;
        root[gene]=gene;
        return gene;
    }
    //std::cout << "finding pos " << gene << std::endl;
    int parent_gene = find_pos(gene, genes[gene], width, height);
    // Two neighbours are pointing at each other. Choosing the one with the highest
    // gene as root.
    //if (print)
    //std::cout << "parent_gene: " << parent_gene << std::endl;
    if (reverse_dir(genes[parent_gene])==genes[gene] && gene>parent_gene) {
        //if (print)
        //std::cout << "parent pointing back and I am older " << std::endl;
        return root[gene] = gene;
    }
    //if (print)
    //std::cout << "returning, max_visited " << max_visited << " gene " << gene  << std::endl;
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

void Individual::mutate(int pos) {
    Dir new_dir = directions[rand()%directions.size()];
    while (new_dir == genes[pos])
        new_dir = directions[rand()%directions.size()];

    int new_root = find_root(pos);
    //std::cout << "mutated " << std::endl;
}