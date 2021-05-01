#ifndef KRUSKAL_H
#define KRUSKAL_H

#include <vector>
#include <iostream>
#include <algorithm>
#include <queue>
#include <opencv2/opencv.hpp>
#include "individual.h"
#include "dir.h"

struct sEdge {
    int u, v; // kant mellom u og v
    int weight ; // vekttallet
    sEdge () {}
    bool operator<( const sEdge &other ) const { // sortere med hensyn paa vekttallet
        return weight < other.weight;
    }
};

class UnionFind {
    private :
        std::vector <int > p, rank, setSize;
        std::vector<sEdge> edges;
        int numSets;
    public:
    UnionFind(int);
    int findSet(int);
    bool isSameSet(int, int);
    void unionSet(int, int);
    int numDisjointSets();
    int sizeOfSet(int);
    void add_edge(int, int, int);
    std::vector<int> kruskal_mst(int, int, int);
    std::pair<std::vector<int>, std::vector<Dir>> kruskal_prim_mst(int, cv::Mat&);
};

#endif