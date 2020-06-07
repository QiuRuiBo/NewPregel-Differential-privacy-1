//
// Created by Jerry Zhang on 2019-06-25.
//

#ifndef GEODISTRIBUTEDGRAPH_UTILS_H
#define GEODISTRIBUTEDGRAPH_UTILS_H

#include <vector>
#include <ctime>
#include <random>
#include <map>
#include "Distributed_graph.h"
namespace utils {
    /*
     * This function generates random permutation from 1 to size.
     */
    std::vector<int> generateRandomPermutation(int size);
	
    /*
     * This function generates random number from 0 to max( inclusive )
     */
    int generateRandomNumber(int max);
	int generateRandomNumberFrom0TOmax(int);
	int findLocation(int a[], int size, int which, Graph& g);  //created by ribo
	vector<int> adjacent_vertices(int,Graph&);  //created by ribo
	float feibonaqi_sum(int n);  //created by ribo
	float feibonaqi_n(float n, float n_front, int k_num);  //created by ribo
	vector<Edge> selfloop_edges(Graph&);
	int* getInPIdistribution(Graph* g);
	int* getOutPIdistribution(Graph* g);
	//std::map<string, float> getPIdistribution(int* in_degree_array, int* out_degree_array, int num_edges, int num_vertices);
    /**
     *
     */
    int* getIn1kSeries(Graph *g);
    int* getOut1kSeries(Graph *g);
}
#endif //GEODISTRIBUTEDGRAPH_UTILS_H
