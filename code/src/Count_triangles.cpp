#include "Count_triangles.h"


//Here count the number of triangles of the input graph G.
long int Triangles::countTriangles(Graph1::Graph g)
{
	//TODO 
	Graph1::Subgraph3Counter SC;
	return (SC.count_subgraphs(g).getTriangles());
	
	//return (*(SC->count_subgraphs(g)).getTriangles);
	//return 100;
}


