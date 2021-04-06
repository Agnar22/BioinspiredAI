#include "prim.h"

Graph::Graph(int V, int w, int h) {
    this->V = V;
    adj = new std::list<iPair> [V];
    width = w;
    height = h;
}

void Graph::add_edge(int u, int v, int w) {
    adj[u].push_back(std::make_pair(v, w));
    adj[v].push_back(std::make_pair(u, w));
}

std::vector<Dir> Graph::prim_mst(int src) {
    std::priority_queue<iPair, std::vector<iPair>, std::greater<iPair>> pq;
    std::vector<int> key(V, INF);
    std::vector<Dir> parent(V, Dir::s);
    std::vector<bool> inMST(V, false);
    pq.push(std::make_pair(0, src));
    key[src] = 0;

    while (!pq.empty()) {
        int u = pq.top().second;
        pq.pop();
        inMST[u] = true;

        std::list<std::pair<int, int> >::iterator i;
        for (i = adj[u].begin(); i != adj[u].end(); ++i) {
            int v = (*i).first;
            int weight = (*i).second;

            if (inMST[v] == false && key[v] > weight) {
                key[v] = weight;
                pq.push(std::make_pair(key[v], v));
                parent[v] = find_direction(v, u, width, height);
            }
        }
    }
    return parent;
}