#include "FindVertexRandomNeighbor.h"
#include <string>
#include "src/Utils.h"
/*According to sign
we can determine whether we want to find the outgoing neighbor vertex 
or the incoming neighbor vertex of the transmitted vertex
*/
int findRandomVertexNeighbor::FindRandomVertexNeighbor(string sign,int vertexNum,Graph g)
{
	string in = "in",out="out";
	int NeighborVertex = 0;
	int finInNeighbor(int,Graph);
	int findOutNeighbor(int,Graph);

	if (strcmp(sign.c_str(), in.c_str()) == 0) {
		/*Find another random vertex if the vertex's out/in edge list is empty*/
		if ((NeighborVertex = finInNeighbor(vertexNum, g)) == -1)  
		{
			std::cout << "Vertex " << vertexNum << "'s In-edge is null." << std::endl;
			//return -1;
		}
	}

	else if (strcmp(sign.c_str(), out.c_str()) == 0) {
		/*Find another random vertex if the vertex's out/in edge list is empty*/
		if ((NeighborVertex = findOutNeighbor(vertexNum, g)) == -1)
		{
			std::cout << "Vertex " << vertexNum << "'s Out-edge is null." << std::endl;
		}
	}

	return NeighborVertex;
}

int finInNeighbor(int vertexNum, Graph g) {
	Graph::in_edge_iterator inedgeIt, inedgeEnd;
	boost::tie(inedgeIt, inedgeEnd) = boost::in_edges(vertexNum,g);  //find the target vertex's in_edge list
	auto findrandomEdge = inedgeIt;
	int vertexIndegree = 0, i = 0,selectEdge=0;
	vertexIndegree = boost::in_degree(vertexNum,g);
	if (vertexIndegree == 0)
		return -1;
	selectEdge = utils::generateRandomNumberFrom0TOmax(vertexIndegree);
	for (findrandomEdge = inedgeIt; i < selectEdge; i++, findrandomEdge++);
	return source(*findrandomEdge, g);

}

int findOutNeighbor(int vertexNum, Graph g) {
	Graph::out_edge_iterator outedgeIt, outedgeEnd;
	boost::tie(outedgeIt, outedgeEnd) = boost::out_edges(vertexNum, g);  //find the target vertex's out_edge list
	auto findrandomEdge = outedgeIt;
	int vertexOutdegree = 0, i = 0, selectEdge=0;
	vertexOutdegree = boost::out_degree(vertexNum,g);
	if (vertexOutdegree == 0)
		return -1;
	selectEdge = utils::generateRandomNumberFrom0TOmax(vertexOutdegree);
	for (findrandomEdge = outedgeIt; i < selectEdge; i++, findrandomEdge++);
	return target(*findrandomEdge,g);
}



