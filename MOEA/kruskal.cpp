#include "kruskal.h"

UnionFind::UnionFind ( int N) {
    setSize.assign (N, 1) ; // er kun 1 i hver gruppe
    numSets = N; //N grupper
    rank.assign (N, 0) ;
    p.assign (N, 0) ;
    for (int i = 0; i < N; i++)
        p[i] = i;
}
int UnionFind::findSet ( int i) {
    return (p[i] == i) ? i : (p[i] = findSet (p[i]) ) ;
    //p[i]= findSet (p[i]) er kun for aa optimalisere ,
    // slik at den husker hvem som er lederen til neste gang
}
bool UnionFind::isSameSet ( int i, int j) { return findSet (i) == findSet (j); }
void UnionFind::unionSet ( int i, int j) {
    if (! isSameSet (i, j)) {
        numSets--;
        int x = findSet (i) , y = findSet (j);
        // rank holder avstanden fra lederen til sine medlemer kort
        if ( rank [x] > rank [y]) {
            p[y] = x;
            setSize [x] += setSize [y];
        } else {
            p[x] = y;
            setSize [y] += setSize [x];
            if (rank [x] == rank [y])
                rank [y ]++;
        }
    }
}
int UnionFind::numDisjointSets () { return numSets ; }
int UnionFind::sizeOfSet ( int i) { return setSize [ findSet (i) ]; }

void UnionFind::add_edge(int u, int v, int w) {
    sEdge e1;
    e1.u = u;
    e1.v = v;
    e1.weight = w;
    edges.push_back(e1);
    sEdge e2;
    e2.u = v;
    e2.v = u;
    e2.weight = w;
    edges.push_back(e2);
}

std::vector<int> UnionFind::kruskal_mst(int min_tree_size, int max_tree_size, int num_trees) {
    std::sort(edges.begin(), edges.end());
    for (int i = 0; i < edges.size() && numDisjointSets()>num_trees; i ++) {
        int size_u = sizeOfSet(edges[i].u);
        int size_v = sizeOfSet(edges[i].v);
        if (!isSameSet(edges[i].u, edges[i].v) && size_u + size_v <= max_tree_size) { // vil ikke danne en sykel
            unionSet(edges[i].u, edges[i].v);
            //minSpanTree.push_back(edges[i]);
        }
    }
    std::cout << "number of sets: " << numDisjointSets() << std::endl;
    for (int i = 0; i < edges.size(); i ++) {
        int size_u = sizeOfSet(edges[i].u);
        int size_v = sizeOfSet(edges[i].v);
        if (!isSameSet(edges[i].u, edges[i].v) && size_u < min_tree_size || size_v < min_tree_size) { // vil ikke danne en sykel
            std::cout << "unified" << std::endl;
            unionSet(edges[i].u, edges[i].v);
        }
    }
    for (int i=0; i<p.size(); ++i)
        findSet(i);
    return p;
}

std::pair<std::vector<int>, std::vector<Dir>> UnionFind::kruskal_prim_mst(int subtrees, cv::Mat &img) {
    int h = img.rows;
    int w = img.cols;
    std::vector<Dir> genes(h*w, Dir::s);
    std::priority_queue<sEdge> pq;
    for (int i=0; i<subtrees; ++i) {
        sEdge e;
        int start = rand()%(h*w);
        e.u = start;
        e.v = start;
        e.weight = 0;
        pq.push(e);
    }

    while (numDisjointSets()>subtrees) {
        sEdge top = pq.top();
        pq.pop();
        if (sizeOfSet(top.v)==1) {
            unionSet(top.u, top.v);
            genes[top.v] = find_direction(top.v, top.u, w, h);
            //std::cout << numDisjointSets() << std::endl;
            int i = top.v/w;
            int j = top.v%w;
            for (Dir dir:directions) {
                if (j==0 && dir==Dir::l || j==w-1 && dir==Dir::r ||i==0 && dir==Dir::u || i==h-1 && dir==Dir::d || dir==Dir::s)
                    continue;
                int neigh = find_pos(top.v, dir, img.cols, img.rows);
                if (sizeOfSet(neigh)!=1)
                    continue;
                int weight = euc_dist(img.at<cv::Vec3b>(i, j), img.at<cv::Vec3b>(neigh/img.cols, neigh%img.cols));
                sEdge e;
                e.u = top.v;
                e.v = neigh;
                e.weight = -weight;
                pq.push(e);
            }
        }
    }
    return std::make_pair(p, genes);
}