#include "Distributed_graph.h"
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <map>
#include <vector>
#include "Utils.h"
#include <string>
#include "Count_triangles.h"
#include "../FindVertexRandomNeighbor.h"
#include "Check_connectivity/Count_subgraphs.h"

#define N 99  
#define GIVEUP_TIMES 10
#define CUTCoefficient_1 1
#define CUTCoefficient_0_9 0.9
#define CUTCoefficient_0_6 0.6
#define CUTCoefficient_0_7 0.7
#define CUTCoefficient_0_5 0.5
#define CUTCoefficient_0_3 0.3
#define CUTCoefficient_0_1 0.1
#define CUTCoefficient_0 0

// [vertex_id1] [vertex_id2]
// NOTE: vertex id should start from 1.
extern int numofdcs;
bool distributed_graph::line_parser(const std::string& filename,const std::string& textline, int ecount) {
	std::stringstream strm(textline);
	size_t source = 0;
	size_t target = 0;
	float weight = 0.0;
	strm >> source;             //The source value is the first row of data(which vertex it connected to) in pawer-law graph
	strm.ignore(1);             //ignore the empty row
	strm >> target;            //The target value is the second row of data (ID) in pawer-law graph
	strm.ignore(1);
	strm >> weight;
	//std::cout << "Jerry: " << source << " " << target <<std::endl;
	
	//if the edge is not self-loop,add this edge into the graph
	if (source != target){
		std::pair<Edge, bool> cur_edge = add_edge(source, target, *(this->g));
		//int weight = boost::get(boost::edge_weight_t(), *(this->g), cur_edge.first);
		boost::put(boost::edge_weight_t(), *(this->g), cur_edge.first, weight);
		//printf("weight=%f\n", weight);
		random_access_edges.emplace(std::make_pair(ecount,cur_edge.first));
	}
	return true;
}





void distributed_graph::load_from_file(std::string filename, int vertices, int edges){
	this->g = new Graph();
	//load vertices	
	std::ifstream graph_file(filename);     //Opens the file *filename in a read-only manner
	std::string line;
	this->num_vertices = vertices;
	this->num_edges = edges;
	std::vector<float> rndlocs;
	std::ifstream rndfile("rndlog");        //Opens the file rndlog in a read-only manner
	std::string rndline;

	//Read the contents of the rndlog file line by line and place them in the rndlocs vector
	for(int i=0; i<1000000; i++){           
		std::getline(rndfile,rndline);      //getline():Read a line from the document rndfile and store it in the rndline  
		std::stringstream strm(rndline);    //stringstream:Convert a string to a number
		float loc;
		strm >> loc;
		rndlocs.push_back(loc);
	}

	//add vertices
	for(int i = 0; i < num_vertices; i++){
		MyVertex* v = new MyVertex();
		v->vertex_id = i;
		//////////////////
		// v->data->input_size = 10;
		//////////////////
		int loc = std::floor(rndlocs[i%1000000]*(float)numofdcs);    //floor(x):Returns the largest integer not greater than x
		//std::cout<< " Jerry: loc: " << loc << std::endl;
		v->location_id = (Amazon_EC2_regions)loc;
		//std::cout << "Qiu:" << "Vertex's location id is:" << v->location_id << std::endl;
		// v->data->routed_nodes.push_back(v->vertex_id);
		add_vertex(*v, *this->g);
	}
	
	//add edges
	std::getline(graph_file, line);
	int ecount = 0;
	while (!graph_file.eof()){  //eof():Determines whether the end of the file has been reached
		line_parser(filename, line, ecount);  //add edge
		std::getline(graph_file, line);
		ecount++;
		if(ecount % 1000000 == 0)
			std::cout << ecount << " edges inserted\n";
	}
	//for(*this->g->operator[i]!=*this->g->out_edge_list)
	/*
	for (int i = 0; i < ecount; i++)
	{

		std::cout << "Qiu:The edges are:" << random_access_edges[i] << std::endl;
	}
	*/
	graph_file.close();
	std::cout << "should have " << edges << " edges," << "added "<< ecount <<" edges." << std::endl;
	//printf("should have %d edges, added %d edges\n",edges,ecount);
	this->num_edges = ecount;
	std::cout << "load graph from file done.\n";
}
//end of load_from_file









void distributed_graph::load_from_file_1k_series(std::string filename, int vertex_num, int edge_num) {
	this->g = new Graph();
	this->original_g = new Graph();
	this->temp_g = new Graph();
	//load vertex_num
	std::ifstream graph_file(filename);
	std::string line;
	this->num_vertices = vertex_num;
	this->num_edges = edge_num;
	std::vector<float> rndlocs;
	std::ifstream rndfile("rndlog");
	std::string rndline;
	for (int i = 0; i < 1000000; i++) {
		std::getline(rndfile, rndline);
		std::stringstream strm(rndline);
		float loc;
		strm >> loc;
		rndlocs.push_back(loc);
	}
	//add vertex_num
	for (int i = 0; i < num_vertices; i++) {
		MyVertex* v = new MyVertex();
		v->vertex_id = i;
		//////////////////
		// v->data->input_size = 10;
		//////////////////
		int loc = std::floor(rndlocs[i % 1000000] * (float)numofdcs);
		//std::cout<< " Jerry: loc: " << loc << std::endl;
		v->location_id = (Amazon_EC2_regions)loc;
		// v->data->routed_nodes.push_back(v->vertex_id);
		add_vertex(*v, *this->g);
	}
	//add edge_num
	std::getline(graph_file, line);
	int ecount = 0;
	while (!graph_file.eof()) {
		line_parser(filename, line, ecount);
		std::getline(graph_file, line);
		ecount++;
		if (ecount % 1000000 == 0)
			std::cout << ecount << " edge_num inserted\n";
	}
	/*
	// Print out original graph's 1k-series
	std::cout << "Original graph in degree distribution: " << std::endl;
	int *original_in = utils::getIn1kSeries(this->g);
	for (int i = 0;i < num_vertices + 1;i++) {
		if (original_in[i] != 0) {
			printf("%d,%d\n", i, original_in[i]);
		}
	}
	// Print out original graph's 1k-series
	std::cout << std::endl << "Original graph out degree distribution: " << std::endl;
	int *original_out = utils::getIn1kSeries(this->g);
	for (int i = 0;i < num_vertices + 1;i++) {
		if (original_in[i] != 0) {
			printf("%d,%d\n", i, original_out[i]);
		}
	}
	*/


	cout << "Original graph has " << boost::num_edges(*g) << " edges" << endl;

	int dc0_numVertices = 0;
	//count the number of vertices in dc0
	for (int index = 0; index < num_vertices; index++) {
		if ((*this->g)[index].location_id == 0)
			dc0_numVertices++;
	}
	*original_g = *g; //save the original graph.
	Count_subgraphs countSubGraph;     //Check whether the generated graph is a connected graph
	// First, iterate through edge_num and remove all the edges that is not in the 0th dc, to obtain the subgraph stored in dc0
	int removedEdgeNum = 0;
	// Stores degree information for other DCs in the original graph
	int* origin_indegree = new int[num_vertices];
	int* origin_outdegree = new int[num_vertices];
	// First int is vertex id, second is the degree.Store the degree information after the disturbance
	std::map<int, int> outDegreeMap;
	std::map<int, int> inDegreeMap;



	/*initial*/
	for (int i = 0; i < num_vertices; i++) {
		origin_indegree[i] = 0;
		origin_indegree[i] = 0;
		outDegreeMap[i] = 0;
		inDegreeMap[i] = 0;
	}


	// Stores degree information in the original graph
	for (int i = 0; i < num_vertices; i++) {
		Vertex v = vertex(i, *this->g);
		int in_index = in_degree(v, *this->g);
		int out_index = out_degree(v, *this->g);
		origin_indegree[i] = in_index;
		origin_outdegree[i] = out_index;
		inDegreeMap[i] = in_index;
		outDegreeMap[i] = out_index;
	}

	/*
	** remove edges that are not in DC0
	*/
	edge_iter ei, ei_end, next, source_next;
	std::tie(ei, ei_end) = edges(*this->g);
	for (next = ei; next != ei_end; ei = next) {
		++next;
		if ((*this->g)[source(*ei, *this->g)].location_id == (*this->g)[target(*ei, *this->g)].location_id
			&& (*this->g)[source(*ei, *this->g)].location_id == 0);  //do nothing
		else
		{
			remove_edge(*ei, *this->g);
			removedEdgeNum++;
		}
	}

	//Update the remaining degree information of the original graph
	for (int i = 0; i < num_vertices; i++) {
		Vertex v = vertex(i, *this->g);
		int in_index = in_degree(v, *this->g);
		int out_index = out_degree(v, *this->g);
		origin_indegree[i] -= in_index;
		origin_outdegree[i] -= out_index;
		inDegreeMap[i] -= in_index;
		outDegreeMap[i] -= out_index;
	}

	//Scrambles the vertices' id of other DCs
	srand((unsigned)time(NULL));
	std::cout << "Jerry : " << removedEdgeNum << " edges removed, " << ecount - removedEdgeNum << " edges left" << std::endl;
	std::vector<int> permutation = utils::generateRandomPermutation(vertex_num - dc0_numVertices); //permutation is a random vector which size from 1 to n.
	for (int i = 0; i < vertex_num - dc0_numVertices; i++)
	{
		int new_index = utils::findLocation(origin_indegree, vertex_num, i, *this->g);
		int old_index = utils::findLocation(origin_indegree, vertex_num, permutation[i], *this->g);
		if (new_index == -1 || old_index == -1)
		{
			std::cout << "Failed to find disturbing position" << std::endl;
		}
		inDegreeMap[new_index] = origin_indegree[old_index];
		outDegreeMap[new_index] = origin_outdegree[old_index];
	}

	// Assignation done, Wire the edges.
	int count = 0;
	int sum = 0;
	for (auto out_pair : outDegreeMap) {
		sum += out_pair.second;
	}
	std::cout << "Jerry: total out degree = " << sum << std::endl;
	count = 0;
	int failed_count = 0;
	// New Wiring method
	srand((unsigned)time(NULL));//Set random number seed, make each generation of random sequences different
	auto out_it = outDegreeMap.begin();
	auto out_end = outDegreeMap.end();
	auto in_it = inDegreeMap.begin();
	auto in_end = inDegreeMap.end();

	/*
	*Save the backup data to restore the original state of the diagram
	*if the composition fails
	*/
	map<int, int> outDegreeMap_backup, inDegreeMap_backup;
	outDegreeMap_backup.clear();
	for (int i = 0; i < outDegreeMap.size(); i++)
		outDegreeMap_backup[i] = outDegreeMap[i];
	inDegreeMap_backup.clear();
	for (int i = 0; i < inDegreeMap.size(); i++)
		inDegreeMap_backup[i] = inDegreeMap[i];
	temp_g->clear();
	*this->temp_g = *this->g;
	int num_edges = ecount;
	int generate_1k_graph_count = 0;
	int numOfSubGraph = 0;
	//the first int is vertex's id,the second is the location of subgraph
	std::map<int, int>* record_subGraph;
	std::map<int, int> record;  //Application address space to record_subGraph,don't use it!!!

	//initial
	record_subGraph = &record;
	for (int i = 0; i < num_vertices; i++)
		record_subGraph->insert(make_pair<int, int>(0, 0));


	std::cout << "The number of subgraph of the synthetic graph in last time is:" << numOfSubGraph << std::endl;
	/*
	*Restore the picture state before 1k composition
	*/
	g->clear();
	*g = *temp_g;
	outDegreeMap.clear();
	for (int i = 0; i < outDegreeMap_backup.size(); i++)
		outDegreeMap[i] = outDegreeMap_backup[i];
	inDegreeMap.clear();
	for (int i = 0; i < inDegreeMap_backup.size(); i++)
		inDegreeMap[i] = inDegreeMap_backup[i];
	out_it = outDegreeMap.begin();
	out_end = outDegreeMap.end();
	in_it = inDegreeMap.begin();
	in_end = inDegreeMap.end();
	num_edges = ecount;
	while (out_it != out_end) {
		if (out_it->second == 0) {
			auto temp = out_it;
			out_it++;
			outDegreeMap.erase(temp);
			if (out_it == out_end)
				break;
			continue;
		}
		in_it = inDegreeMap.begin();
		while (1) {
			if (in_it == in_end)
			{
				in_it = inDegreeMap.begin();
				if (in_it == inDegreeMap.end())
				{
					out_it++;
					break;
				}
			}
			if (in_it->second == 0) {
				auto temp = in_it;
				in_it++;
				inDegreeMap.erase(temp);
				continue;
			}
			if (out_it->second == 0) {
				//out_it++;
				break;
			}
			//ADD PI distribution to wire edges.
			double p = (double)(((double)out_it->second * in_it->second / (double)num_edges));
			double random = rand() % (N + 1) / (double)(N + 1);


			//std::cout <<"P=" <<p <<"Random=" << random << std::endl;
			//std::cout << out_it->first << std::endl;
			//std::cout << edge(out_it->first, in_it->first, *this->g).first << edge(out_it->first, in_it->first, *this->g).second<< std::endl;


			//if (p > random)    //Run the following code with probability p
			{
				//The adding edge is already exit if return false.
				bool success = add_edge(out_it->first, in_it->first, *this->g).second;
				if (success) {
					count++;
					num_edges--;
					//std::cout << "add " << count << " edges." << std::endl;
					//std::cout << "ADD success." << std::endl;
					out_it->second--;
					in_it->second--;
					if (out_it->second < 0 || in_it->second < 0) {
						std::cout << "ERROR:" << std::endl;
					}
				}
				else  //rewire edge
				{
					size_t randomVertex, secondVertex;
					std::tie(ei, ei_end) = edges(*this->g);
					findRandomVertexNeighbor  FRVN;
					while (1) {
						//srand((unsigned)time(NULL));
						randomVertex = utils::generateRandomNumberFrom0TOmax(boost::num_edges(*this->g));
						//rewire the edge.
						int i = 0;
						//find the target vertex,and find the vertex's random out edge
						for (next = ei, i = 0; i < randomVertex; i++)next++;
						randomVertex = source(*next, *this->g);
						secondVertex = FRVN.FindRandomVertexNeighbor("out", randomVertex, *this->g);
						if (secondVertex == -1)  //vertex has not any out edge
							continue;
						if ((*this->g)[randomVertex].location_id == 0 && (*this->g)[secondVertex].location_id == 0)
						{
							//printf("Selected edge in DC0,retry.........\n");
							continue;
						}
						else if (edge(randomVertex, secondVertex, *this->g).second == false)
							continue;
						else if (edge(out_it->first, secondVertex, *this->g).second == true)
							continue;
						else if (edge(randomVertex, in_it->first, *this->g).second == true)
							continue;
						else
							break;
					}

					remove_edge(randomVertex, secondVertex, *this->g);
					add_edge(out_it->first, secondVertex, *this->g);
					add_edge(randomVertex, in_it->first, *this->g);
					//std::cout << "Rewire edge sucess." << std::endl;
					count++;
					num_edges--;
					out_it->second--;
					in_it->second--;
					if (out_it->second < 0 || in_it->second < 0) {
						std::cout << "ERROR:" << std::endl;
						if (out_it == out_end)
							break;
					}
				}
				in_it++;
			}

		}
		generate_1k_graph_count++;
		//numOfSubGraph = countSubGraph.operator()(*this->g,*record_subGraph);
	}







	numOfSubGraph = countSubGraph.operator()(*this->g, *record_subGraph);
	std::cout << "The number of subgraph of PI distribution graph is:" << numOfSubGraph << std::endl;
	if(numOfSubGraph != 1)
		std::cout << "Now start reconnecting the edges so that the synthetic graph is a connected graph" << std::endl;
	findRandomVertexNeighbor  FRVN;
	//std::cout << "-----------------The number of subgraph of boost is:" << countSubGraph.operator_byboost(*this->g) << std::endl;
	//preprocessed
	map<int, int> subgraph_count;  //Record the number of vertices per subgraph
	for (int i = 0; i < num_vertices; i++)
	{
		subgraph_count[(*record_subGraph)[i]]++;
	}
	/*
	for (auto next = subgraph_count.begin(); next != subgraph_count.end(); next++)
	{
		printf("Subgraph %d has %d vertices.\n", next->first, next->second);
		if (next->second == 1)
		{
			//out_degree()
			printf("Degree distribution of orphan vertices:out=%d,in=%d.\n", out_degree(vertex(next->first, *g), *g), in_degree(vertex(next->first, *g), *g));
			printf("orphan vertex %d out edge=%d,in edge=%d\n", next->first, FRVN.FindRandomVertexNeighbor("out", next->first, *g), FRVN.FindRandomVertexNeighbor("in", next->first, *g));
		}
	}
	*/

	//rewire
	srand((unsigned)time(NULL));
	while (numOfSubGraph != 1)
	{
		int indexofsubgraph = utils::generateRandomNumberFrom0TOmax(num_vertices - 1);
		int numOfSubGraph_temp = 0;
		bool reset = false;
		while ((*record_subGraph)[indexofsubgraph] == 1)
		{
			indexofsubgraph++;
			if (indexofsubgraph >= num_vertices)  //Prevents arrays from crossing bounds
			{
				reset = true;
				break;
			}
		}
		if (reset)
			continue;
		size_t firVerInSub = boost::vertex(indexofsubgraph, *g);
		size_t secondVerInSub = FRVN.FindRandomVertexNeighbor("out", firVerInSub, *this->g);
		if (secondVerInSub == -1)  //vertex has not any out edge
		{
			continue;
		}
		if (edge(firVerInSub, secondVerInSub, *this->g).second == false)
		{
			continue;  //if the edge is not exist,choose again.
		}
		size_t randVerInOri = utils::generateRandomNumberFrom0TOmax(num_vertices - 1);
		if ((*record_subGraph)[randVerInOri] != 1)
			reset = false;
		while ((*record_subGraph)[randVerInOri] != 1)
		{
			randVerInOri++;
			if (randVerInOri >= num_vertices)  //Prevents arrays from crossing bounds
			{
				reset = true;
				break;
			}
		}
		if (reset)
			continue;
		randVerInOri = boost::vertex(randVerInOri, *g);  //change int to vertex
		size_t secondVerInOri = FRVN.FindRandomVertexNeighbor("out", randVerInOri, *this->g);
		if (secondVerInOri == -1)  //vertex has not any out edge
		{
			continue;
		}
		if (edge(firVerInSub, secondVerInSub, *this->g).second == false)
		{
			continue;  //if the edge is not exist,choose again.
		}

		if ((*this->g)[firVerInSub].location_id == 0 && (*this->g)[secondVerInSub].location_id == 0)
		{
			printf("selected dc0's edge when control subgraph.\n");
			continue;
		}

		if (edge(firVerInSub, secondVerInOri, *this->g).second == true)
		{
			//printf("edge(firVerInSub, secondVerInOri, *this->g) should not be exist");
			continue;  //if the edge is exist,choose again.
		}
		if (edge(randVerInOri, secondVerInSub, *this->g).second == true)
		{
			//printf("edge(randVerInOri, secondVerInSub, *this->g) should not be exist");
			continue;  //if the edge is exist,choose again.
		}
		///////////////////////Selection edge complete////////////////////////////////////

		remove_edge(firVerInSub, secondVerInSub, *this->g);
		remove_edge(randVerInOri, secondVerInOri, *this->g);

		add_edge(firVerInSub, secondVerInOri, *this->g).second;
		add_edge(randVerInOri, secondVerInSub, *this->g).second;

		numOfSubGraph_temp = countSubGraph.operator()(*this->g, *record_subGraph);


		if (numOfSubGraph_temp < numOfSubGraph)
		{
			numOfSubGraph = numOfSubGraph_temp;
			std::cout << "The number of subgraph is:" << numOfSubGraph << std::endl;
		}
		else  //if the rewire is useless,drop it
		{
			remove_edge(firVerInSub, secondVerInOri, *this->g);
			remove_edge(randVerInOri, secondVerInSub, *this->g);

			add_edge(firVerInSub, secondVerInSub, *this->g).second;
			add_edge(randVerInOri, secondVerInOri, *this->g).second;

			countSubGraph.operator()(*this->g, *record_subGraph);

			//printf("rewire is useless,drop it.\n");
		}
	}

	//std::cout << "The edges are " << boost::num_edges(*g) << " after control subgraph." << std::endl;


	//Rewiring edges make the number of triangles the same.
	Triangles* CT = new Triangles();
	int Original_TRIANGLE = 0, Synthetic_TRIANGLE = 0;
	Graph1::Graph gg;
	gg._g = *this->original_g;
	Original_TRIANGLE = CT->countTriangles(gg);
	gg._g = *this->g;
	Synthetic_TRIANGLE = CT->countTriangles(gg);
	std::cout << "The number of triangles of original graph is:" << Original_TRIANGLE << std::endl;
	std::cout << "The number of triangles of synthetic graph is:" << Synthetic_TRIANGLE << std::endl;
	std::cout << "Now begin to rewire the edges" << std::endl;

	std::cout << "Before rewire:num_edges=" << boost::num_edges(*this->g) << std::endl;
	int rewireEdge = 0;
	//create at 9_24 by Ribo to test the time consumer line.
	int executeTimes = 0, numSubgraph = 0;
	size_t randomVertex;
	size_t target_neighbor = 0, source_neighbor = 0, fourthVertex = 0, secondVertex = 0;
	bool little = false, remove = false; //bigger = false;
	int giveup = 0;
	float random, p;
	//while (Synthetic_TRIANGLE < Synthetic_TRIANGLE*1.2 && Synthetic_TRIANGLE< Original_TRIANGLE*0.86) {
	//while (Synthetic_TRIANGLE < Original_TRIANGLE || Synthetic_TRIANGLE > Original_TRIANGLE) {
	while (Synthetic_TRIANGLE < 0) {     //choose this line to represent the number of triangles not reconstructed
		std::tie(ei, ei_end) = edges(*this->g);  //the graph g will change every circle
		if (Synthetic_TRIANGLE < Original_TRIANGLE)
			little = true;
		else
			little = false;
		//select a random edge from g.
		randomVertex = utils::generateRandomNumberFrom0TOmax(boost::num_edges(*this->g));
		//If the neighbor of neighbor of a random vertex has an outdegree greater than 2,
		//rewire the edge.
		int i = 0;
		//find the target vertex,and find the vertex's random out edge
		for (next = ei, i = 0; i < randomVertex; i++)next++;

		randomVertex = source(*next, *this->g);

		secondVertex = FRVN.FindRandomVertexNeighbor("out", randomVertex, *this->g);
		if (secondVertex == -1)  //vertex has not any out edge
			continue;
		if (edge(randomVertex, secondVertex, *this->g).second == false)
			continue;

		target_neighbor = FRVN.FindRandomVertexNeighbor("out", secondVertex, *this->g);
		if (target_neighbor == -1)  //vertex has not any out edge
			continue;
		if (edge(secondVertex, target_neighbor, *this->g).second == false)
			continue;

		source_neighbor = FRVN.FindRandomVertexNeighbor("in", randomVertex, *this->g);
		if (source_neighbor == -1)  //vertex has not any out edge
			continue;
		if (edge(source_neighbor, randomVertex, *this->g).second == false)
			continue;

		if (little == true && out_degree(target_neighbor, *this->g) >= 2)
		{
			fourthVertex = FRVN.FindRandomVertexNeighbor("out", target_neighbor, *this->g);
			if (fourthVertex == -1)  //vertex has not any out edge
				continue;
			if (edge(target_neighbor, fourthVertex, *this->g).second == false)
				continue;

			/*
			*Exit this loop if the deleted edges needed to reconstruct the
			*number of triangles already exist
			*/
			if (edge(target_neighbor, randomVertex, *g).second == true
				|| edge(source_neighbor, fourthVertex, *g).second == true)
				continue;
			/*
			if ((*this->g)[target_neighbor].location_id == 0 && (*this->g)[fourthVertex].location_id == 0)
			{
				printf("Selected edge in DC0,retry.........\n");
				continue;
			}
			if ((*this->g)[source_neighbor].location_id == 0 && (*this->g)[randomVertex].location_id == 0)
			{
				printf("Selected edge in DC0,retry.........\n");
				continue;
			}
			*/
			//remove two edge
			remove_edge(target_neighbor, fourthVertex, *this->g);
			remove_edge(source_neighbor, randomVertex, *this->g);

			//add two edge
			bool addSucess_source = false, addSucess_target = false;
			addSucess_source = add_edge(target_neighbor, randomVertex,
				*this->g).second;
			addSucess_target = add_edge(source_neighbor, fourthVertex,
				*this->g).second;

			//if the number of triangles are more than the original graph,then save the rewiring
			//else undo the rewiring.
			int temp;
			gg._g = *this->g;
			temp = CT->countTriangles(gg);
			numSubgraph = countSubGraph.operator()(*this->g, *record_subGraph);
			if (temp > Synthetic_TRIANGLE && numSubgraph == 1)
			{
				Synthetic_TRIANGLE = temp;
				if (addSucess_source == true && addSucess_target == true)
					std::cout << "Rewire " << ++rewireEdge << " edge." << std::endl;
				std::cout << "The total number of times selected for this round of refactoring is:" << executeTimes << std::endl;
				std::cout << "The total number of triangles is:" << temp << std::endl;

				executeTimes = 0;
			}
			else
			{
				/*
				*Restore the graph
				*/

				/*
				if ((*this->g)[target_neighbor].location_id == 0 && (*this->g)[randomVertex].location_id == 0)
				{
					printf("Selected edge in DC0,retry.........\n");
					continue;
				}
				if ((*this->g)[source_neighbor].location_id == 0 && (*this->g)[fourthVertex].location_id == 0)
				{
					printf("Selected edge in DC0,retry.........\n");
					continue;
				}
				*/
				//remove two edge
				remove_edge(target_neighbor, randomVertex, *this->g);
				remove_edge(source_neighbor, fourthVertex, *this->g);

				//add two edge
				addSucess_source = add_edge(target_neighbor, fourthVertex,
					*this->g).second;
				addSucess_target = (add_edge(source_neighbor, randomVertex,
					*this->g)).second;
				executeTimes++;
			}
		}


		/*
		*If the selected vertex is a vertex in the triangle
		*then broke the triangle structural
		*/
		if (little == false)
		{
			//find the random next vertex and delete the edge
			fourthVertex = FRVN.FindRandomVertexNeighbor("in", target_neighbor, *this->g);
			if (fourthVertex == -1)  //vertex has not any out edge
				continue;
			if (edge(fourthVertex, target_neighbor, *this->g).second == false)
				continue;

			/*
			*Exit this loop if the deleted edges needed to reconstruct the
			*number of triangles already exist
			*/
			if (edge(fourthVertex, secondVertex, *g).second == true
				|| edge(randomVertex, target_neighbor, *g).second == true)
				continue;

			/*
			if ((*this->g)[randomVertex].location_id == 0 && (*this->g)[secondVertex].location_id == 0)
			{
				printf("Selected edge in DC0,retry.........\n");
				continue;
			}
			if ((*this->g)[fourthVertex].location_id == 0 && (*this->g)[target_neighbor].location_id == 0)
			{
				printf("Selected edge in DC0,retry.........\n");
				continue;
			}
			*/
			//remove two random edges
			remove_edge(randomVertex, secondVertex, *this->g);
			remove_edge(fourthVertex, target_neighbor, *this->g);

			//add two edges
			bool addSucess_source, addSucess_target;
			addSucess_source = add_edge(fourthVertex, secondVertex, *this->g).second;
			addSucess_target = add_edge(randomVertex, target_neighbor, *this->g).second;

			//if the number of triangles are more than the original graph,then save the rewiring
			//else undo the rewiring.
			int temp;
			gg._g = *this->g;
			temp = CT->countTriangles(gg);
			numSubgraph = countSubGraph.operator()(*this->g, *record_subGraph);
			if (temp < Synthetic_TRIANGLE && numSubgraph == 1)
			{
				Synthetic_TRIANGLE = temp;
				if (addSucess_source == true && addSucess_target == true)
					std::cout << "Rewire " << ++rewireEdge << " edge." << std::endl;
				std::cout << "The total number of times selected for this round of refactoring is:" << executeTimes << std::endl;
				std::cout << "The total number of triangles is:" << temp << std::endl;
				executeTimes = 0;
			}
			else
			{
				/*
				*Restore the graph
				*/

				//remove two edge
				remove_edge(randomVertex, target_neighbor, *this->g);
				remove_edge(fourthVertex, secondVertex, *this->g);

				//add two edge
				add_edge(randomVertex, secondVertex, *this->g);
				add_edge(fourthVertex, target_neighbor, *this->g);
				executeTimes++;
			}
		}
	}
	std::cout << "The number of triangles of synthetic graph is:" << CT->countTriangles(gg) << std::endl;

	cout << "synthetic graph has " << boost::num_edges(*g) << " edges" << endl;




	// Rewiring self loop
	
	std::vector<Edge> self_loop_edges;
	int index1 = 0;
	int index2 = 1;
	bool reset = false;
	self_loop_edges = utils::selfloop_edges(*g);  //find the self-loop edges in the graph
	printf("Initial have %d self-loop edges.\n", (int)self_loop_edges.size());
	while (!self_loop_edges.empty()) {
		if (self_loop_edges.size() == 1) {
			std::cout << " Only One self loop edge in the graph.Delete it directly." << std::endl;
			remove_edge(self_loop_edges[0], *g);
			std::vector<Edge>::iterator it = self_loop_edges.begin();
			self_loop_edges.erase(it + 0);
		}
		if (self_loop_edges.size() == 2) {
			std::cout << " Two self loop edge in the graph." << std::endl;
			bool failed = false;
			if (edge(source(self_loop_edges[0], *this->g), source(self_loop_edges[1], *this->g), *this->g).second)
				failed = true;
			if (edge(source(self_loop_edges[1], *this->g), source(self_loop_edges[0], *this->g), *this->g).second)
				failed = true;
			if (failed == false) {
				add_edge(source(self_loop_edges[0], *this->g), source(self_loop_edges[1], *this->g), *this->g);
				add_edge(source(self_loop_edges[1], *this->g), source(self_loop_edges[0], *this->g), *this->g);
			}
			//Regardless of whether the edge energy already exists, remove the edge from the loop
			remove_edge(self_loop_edges[0], *g);
			remove_edge(self_loop_edges[1], *g);
			std::vector<Edge>::iterator it = self_loop_edges.begin();
			self_loop_edges.erase(it + 0);
			self_loop_edges.erase(it + 1);
		}
		if (self_loop_edges.size() > 2) {
			index2 = 0;
			while (index2 < self_loop_edges.size()) {
				while (edge(source(self_loop_edges[index1], *this->g), source(self_loop_edges[index2], *this->g), *this->g).second)
				{
					index2++;  //If the target edge already exists, do index2 and go down
					if (index2 >= self_loop_edges.size())
					{
						reset = true;
						break;
					}
				}
				if (reset)
				{
					index1++;
					continue;
				}
				bool sucess = add_edge(source(self_loop_edges[index1], *this->g), source(self_loop_edges[index2], *this->g), *this->g).second;
				if (sucess)
					//printf("add edge sucess!\n");
				if(edge(source(self_loop_edges[index1], *this->g), source(self_loop_edges[index1], *this->g), *this->g).second)
					remove_edge(self_loop_edges[index1], *this->g);
				if (edge(source(self_loop_edges[index2], *this->g), source(self_loop_edges[index2], *this->g), *this->g).second)
					remove_edge(self_loop_edges[index2], *this->g);

				count++;
				index1++;
				index2++;
			}
			if (index1 >= self_loop_edges.size())
			{
				for (int i = 0; i < self_loop_edges.size(); i++)
				{
					remove_edge(self_loop_edges[i], *g);
				}
				printf("Self - loop resolution failed,Delete all self - looping edges directly.\n");
			}
			add_edge(source(self_loop_edges[index2 - 1], *this->g), source(self_loop_edges[0], *this->g), *this->g);
			count++;
		}
		self_loop_edges = utils::selfloop_edges(*g);  //renew the self-loop edges in the graph
	}

	printf("Now have %d self-loop edges.\n", (int)self_loop_edges.size());
	std::cout << "After deal with self-loop.Here are " << boost::num_edges(*g) << " edges." << std::endl;
		std::ofstream save("synthetic_g.txt");
		std::tie(ei, ei_end) = edges(*this->g);
		for (next = ei; next != ei_end; next++) {

			save << source(*next, *this->g) << " " << target(*next, *this->g) << std::endl;
		}
		save.close();
		graph_file.close();
		this->num_edges = ecount;
		std::cout << "load graph from file and split into 1k-series done.\n";

}//end of load_from_file_1k_series







void distributed_graph::load_synthetic_file(const char* filename){
	//load live journal graph from files
	if(strcmp(filename,"powerlaw-1000-2.1")==0){
		this->num_vertices = 1000;
	}
	else if(strcmp(filename,"powerlaw-10000-2.1")==0){
		this->num_vertices = 10000;
	}
	else if(strcmp(filename,"powerlaw-50000-2.1")==0){
		this->num_vertices = 50000;
	}else if(strcmp(filename,"powerlaw-100000-2.1")==0){
		this->num_vertices = 100000;
	}else if(strcmp(filename,"powerlaw-1000000-2.1")==0){
		this->num_vertices = 1000000;
	}else{
		printf("please select from 10000, 50000 and 100000 vertices.\n");
		exit(1);
	}
	this->load_from_file(filename,this->num_vertices,0);
}
//end of load_synthetic_file








void distributed_graph::load_synthetic_file_1k_series(const char *filename) {
    //load live journal graph from files
    if(strcmp(filename,"powerlaw-1000-2.1")==0){
        this->num_vertices = 1000;
    }
    else if(strcmp(filename,"powerlaw-10000-2.1")==0){
        this->num_vertices = 10000;
    }
    else if(strcmp(filename,"powerlaw-50000-2.1")==0){
        this->num_vertices = 50000;
    }else if(strcmp(filename,"powerlaw-100000-2.1")==0){
        this->num_vertices = 100000;
    }else if(strcmp(filename,"powerlaw-1000000-2.1")==0){
        this->num_vertices = 1000000;
    }else{
        printf("Error:please select from 10000, 50000 and 100000 vertices to argv[3].\n");
        exit(1);
    }
    this->load_from_file_1k_series(filename,this->num_vertices,0);
}
//end of load_synthetic_file_1k_series








//synthetic graphs
void distributed_graph::load_synthetic_powerlaw(size_t nverts, bool in_degree,	double alpha, size_t truncate){
	this->num_vertices = nverts;
	this->g = new Graph();
	std::vector<double> prob(std::min(nverts, truncate), 0);
	std::cout << "constructing pdf" << std::endl;
	for (size_t i = 0; i < prob.size(); ++i)
		prob[i] = std::pow(double(i + 1), -alpha);
	std::cout << "constructing cdf" << std::endl;
	pdf2cdf(prob);
	std::cout << "Building graph" << std::endl;	
	size_t target_index = 0;
	size_t addedvtx = 0;
	// add vertices
	for (int vi = 0; vi < nverts; vi++){
		MyVertex* v = new MyVertex();
		v->vertex_id = vi;//starting from 0	
		//////////////////
		// v->data->input_size = 10;
		//////////////////		
		v->location_id = (Amazon_EC2_regions)((int)this->rnd_generator(0, numofdcs)); 
		// v->data->routed_nodes.push_back(v->vertex_id);
		add_vertex(*v, *this->g);
	}
	// A large prime number
	const size_t HASH_OFFSET = 2654435761;
	for (size_t source = 0; source < nverts; source++) {
		const size_t out_degree = multinomial_cdf(prob) + 1;
		//printf("# out edges is %d\n", out_degree);
		for (size_t i = 0; i < out_degree; ++i) {
			target_index = (target_index + HASH_OFFSET) % nverts;
			while (source == target_index) {
				target_index = (target_index + HASH_OFFSET) % nverts;
			}
			float weight = this->rnd_generator(0,1000);
			if (in_degree) {
				std::pair<Edge, bool> cur_edge = add_edge(target_index, source, weight, *this->g);
				random_access_edges.emplace(addedvtx,cur_edge.first);
			}
			else {
				std::pair<Edge, bool> cur_edge = add_edge(source, target_index, weight, *this->g);
				random_access_edges.emplace(addedvtx,cur_edge.first);
			}
		}
		++addedvtx;
		if (addedvtx % 1000000 == 0) {

			std::cout << addedvtx << " inserted\n";
		}
	}
	this->num_edges = addedvtx;
} // end of load random powerlaw
//end of load_synthetic_powerlaw








//real-world graphs
/* void distributed_graph::load_live_journal_graph(){
	this->g = new Graph();
	//load live journal graph from files
	//load vertices
	std::string location_file_name = "livejournal/id2region";
	std::ifstream location_file(location_file_name);
	std::string line;
	this->num_vertices = 3577167;
	std::getline(location_file, line);
	std::stringstream strm(line);
	int curr_id = 0;
	strm >> curr_id;
	for (int i = 0; i < num_vertices; i++){
		MyVertex* v = new MyVertex();
		v->vertex_id = i;
		//////////////////
		v->data->input_size = 10;
		//////////////////
		std::string location_name;
		while (curr_id < v->vertex_id+1){
			if(!location_file.eof()){
				std::getline(location_file, line);
				strm = std::stringstream(line);
				strm >> curr_id;
			}else{
				break;
			}
			
		}
		if (curr_id == v->vertex_id+1){
			strm.ignore(1);
			std::getline(strm, location_name);
		}
		else if (curr_id > v->vertex_id+1){
			// cannot find userid i in the file
			// use default location for this user			
			location_name = "North America";
		}else{
			location_name = "North America";
		}
		v->location_id = location_to_region(location_name, "livejournal");
		v->data->routed_nodes.push_back(v->vertex_id);
		add_vertex(*v, *g);
		if (i % 100000 == 0) {
			std::cout << i << " vertices inserted\n";
		}
	}
	location_file.close();
	//load edges
	std::string graph_file_name = "livejournal/graph";
	std::ifstream graph_file(graph_file_name);
	std::getline(graph_file, line);
	int ecount = 0;
	while (!graph_file.eof()){
		line_parser(graph_file_name, line, ecount);
		std::getline(graph_file, line);
		ecount ++;	
		if(ecount % 100000 == 0)
			std::cout << ecount << " edges inserted\n";		
	}
	graph_file.close();
	
	std::cout << "load live journal done.\n";
} */


void distributed_graph::load_twitter_graph(){
	std::string filename = "/tmp/twitter_rv.net";
	int v = 61578415; //41652230; //81306; 
	int e = 1468365182; //1768149;
	this->load_from_file(filename,v,e);
}

void distributed_graph::load_googleweb_graph(){
	std::string filename = "googleweb/web-Google.txt";
	int v = 916427; //875713;
	int e = 5105039;
	this->load_from_file(filename,v,e);
}

void distributed_graph::load_facebook_graph(){
	std::string filename = "facebook/facebook_combined.txt";
	int v = 4039; 
	int e = 88234;
	this->load_from_file(filename,v,e);
}

void distributed_graph::load_wiki_graph(){
	std::string filename = "wikivote/Wiki-Vote.txt";
	int v = 8298; //7115;
	int e = 103689;
	this->load_from_file(filename,v,e);
}

void distributed_graph::load_roadnet_graph(){
	std::string filename = "roadca/weight-roadNet-CA.txt";
	int v = 1971281; //1965206;
	int e = 5533214;
	this->load_from_file(filename,v,e);
}

void distributed_graph::load_live_journal_graph(){
	std::string filename = "livejournal/graph";
	int v = 3577167; 
	int e = 44913072;
	this->load_from_file(filename,v,e);
}

void distributed_graph::load_p2p_graph(){
	std::string filename = "gnutella/p2p-Gnutella09.txt";
	int v = 8114; //8846;
	int e = 26013; //31839;
	this->load_from_file(filename,v,e);
}
/**
* vertices: 0,1,2,3
* edges: 0->1;0->2;1->3;2->3
*/
void distributed_graph::test_graph(){
	this->num_vertices = 4;
	this->num_edges = 4;
	this->g = new Graph();
	for (int i = 0; i < num_vertices; i++){
		MyVertex* v = new MyVertex();
		v->vertex_id = i;//starting from 0		
		//////////////////
		// v->data->input_size = 10;
		//////////////////		
		v->location_id = (Amazon_EC2_regions)((int)this->rnd_generator(0, numofdcs));
		// v->data->routed_nodes.push_back(v->vertex_id);
		add_vertex(*v, *this->g);
	}
	
	std::pair<Edge, bool> cur_edge = add_edge(0, 1, 10, *g);
	random_access_edges.emplace(0,cur_edge.first);
	cur_edge = add_edge(0, 2, 25, *g);
	random_access_edges.emplace(1,cur_edge.first);
	cur_edge = add_edge(1, 3, 20, *g);
	random_access_edges.emplace(2,cur_edge.first);
	cur_edge = add_edge(2, 3, 15, *g);
	random_access_edges.emplace(3,cur_edge.first);
}
