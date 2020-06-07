//
// Created by Jerry Zhang on 2019-09-16.
//

#include "Count_subgraphs.h"
#include <iostream>
#include <stack>
#include <boost/graph/connected_components.hpp>
#include "../utils.h"
//#include <unistd.h>
using namespace std;
using namespace boost;

int Count_subgraphs::operator_byboost(Graph& g)  //Only for undirected graph
{
	std::vector<int> component(num_vertices(g));
	std::cout << "Connected components = " << connected_components(g, &component[0]) << std::endl;
	return connected_components(g, &component[0]);
}
int Count_subgraphs::operator()(Graph &g,map<int,int> &record) {  //Use for directed or undirected graph
    int vertex_count = num_vertices(g);
    bool *visit = new bool[vertex_count];
    for (int i = 0;i < vertex_count;i++) {
        visit[i] = false;
    }
    stack<int> to_visit;
    to_visit.push(0);
    int subgraph_count = 0;
    bool done_flag = false;
	
	int count = 0,count1=0;
    while (!done_flag) {
        subgraph_count++;
        while (!to_visit.empty()) {
            int v = to_visit.top();
            to_visit.pop();
			if (visit[v])
			{
				continue;
			}
            visit[v] = true;
			(record)[v] = subgraph_count;
			//printf("Vertex %d divided into subgraph %d.\n", v, record[v]);
			size_t vertex;
			vertex = boost::vertex(v, g);
			vector<int> adjacent_vertices = utils::adjacent_vertices(v, g);
			for (int i = 0; i<adjacent_vertices.size(); i++) {
				size_t vd = adjacent_vertices[i];
                to_visit.push(vd);
				/*
				if (edge(vd, v, g).second == true)
					printf("%d is in edge.\n",vd);
				if (edge(v, vd, g).second == true)
					printf("%d is out edge.\n",vd);
				*/
            }
        }
        done_flag = true;
        for (int i = 0;i < vertex_count;i++) {
            if (!visit[i]) {
                to_visit.push(i);
                done_flag = false;
                break;
            }
        }
    }
    return subgraph_count;
}
