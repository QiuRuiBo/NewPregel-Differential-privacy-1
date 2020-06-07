//
// Created by Jerry Zhang on 2019-09-16.
//
#include "../Distributed_graph.h"
#ifndef GEODISTRIBUTEDGRAPH_COUNT_SUBGRAPHS_H
#define GEODISTRIBUTEDGRAPH_COUNT_SUBGRAPHS_H


class Count_subgraphs {
public:
	int operator()(Graph&, map<int, int>&);
	int operator_byboost(Graph&);
private:
};


#endif //GEODISTRIBUTEDGRAPH_COUNT_SUBGRAPHS_H
