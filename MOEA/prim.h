#ifndef PRIM_H
#define PRIM_H

#include<bits/stdc++.h>
#include "dir.h"
#define INF 0x3f3f3f3f

// iPair ==>  Integer Pair
typedef std::pair<int, int> iPair;

class Graph {
    int V, width, height;
    std::list<std::pair<int, int>> *adj;

public:
    Graph(int, int, int);
    void add_edge(int, int, int);
    std::vector<Dir> prim_mst(int, int);
};

#endif