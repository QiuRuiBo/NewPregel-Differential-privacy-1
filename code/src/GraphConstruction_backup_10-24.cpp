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
	//std::cout << "Jerry: " << source << " " << target <<std::endl;
	
	//if the edge is not self-loop,add this edge into the graph
	if (source != target){
		std::pair<Edge, bool> cur_edge = add_edge(source, target, *(this->g));
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
	while (!graph_file.eof()){                              //eof():Determines whether the end of the file has been reached
		line_parser(filename, line, ecount);          //add edge
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









void distributed_graph::load_from_file_1k_series(std::string filename, int vertex_num, int edge_num){
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
    for(int i=0; i<1000000; i++){
        std::getline(rndfile,rndline);
        std::stringstream strm(rndline);
        float loc;
        strm >> loc;
        rndlocs.push_back(loc);
    }
    //add vertex_num
    for(int i = 0; i < num_vertices; i++){
        MyVertex* v = new MyVertex();
        v->vertex_id = i;
        //////////////////
        // v->data->input_size = 10;
        //////////////////
        int loc = std::floor(rndlocs[i%1000000]*(float)numofdcs);
        //std::cout<< " Jerry: loc: " << loc << std::endl;
        v->location_id = (Amazon_EC2_regions)loc;
        // v->data->routed_nodes.push_back(v->vertex_id);
        add_vertex(*v, *this->g);
    }
    //add edge_num
    std::getline(graph_file, line);
    int ecount = 0;
    while (!graph_file.eof()){
        line_parser(filename, line, ecount);
        std::getline(graph_file, line);
        ecount++;
        if(ecount % 1000000 == 0)
            std::cout << ecount << " edge_num inserted\n";
    }
	/*
    // Print out original graph's 1k-series
    std::cout << "Original graph in degree distribution: " << std::endl;
    int *original_in = utils::getIn1kSeries(this->g);                   
    for (int i = 0;i < num_vertices + 1;i++) {
        if (original_in[i] != 0) {
            cout << i << "," << original_in[i] << std::endl;
        }
    }
    std::cout << std::endl << "Original graph out degree distribution: " << std::endl;
    int *original_out = utils::getIn1kSeries(this->g);
    for (int i = 0;i < num_vertices + 1;i++) {
        if (original_in[i] != 0) {
            cout << i << "," << original_out[i] << std::endl;
        }
    }
	*/

	//std::cout << std::endl << "PI distribution: " << std::endl;
	//int* original_out = utils::getInPIdistribution(this->g);
	//int* original_in = utils::getOutPIdistribution(this->g);
	//output the original graph's degree distribution in each vertex
	/*
	for (int i = 0; i < num_vertices;i++) {
		std::cout << "In:" << original_in[i] << std::endl;
	}
	
	for (int i = 0; i < num_vertices; i++) {
		std::cout<< "Out:" << original_out[i] << std::endl;
	}
	std::cout << "Total:" << num_vertices << "In vertices" << std::endl;
	*/

    /**
     * Handling 1k series
     */

	*original_g = *g; //save the original graph.
	Count_subgraphs countSubGraph;     //Check whether the generated graph is a connected graph
     // First, iterate through edge_num and remove all the edges that is not in the 0th dc, to obtain the subgraph stored in dc0
    int removedEdgeNum = 0;
    // This stores how much nodes were there for each degree
    int *inDegreeAssignation = new int[num_vertices];
    int *outDegreeAssignation = new int[num_vertices];
    int *saved_inDegreeAssignation = new int[num_vertices];
    int *saved_outDegreeAssignation = new int[num_vertices];
    int *synthetic_inDegreeAssignation = new int[num_vertices];
    int *synthetic_outDegreeAssignation = new int[num_vertices];
    std::map<int,std::vector<int> > inDegreeAssignation_vertex;
    std::map<int,std::vector<int> > outDegreeAssignation_vertex;
    // First int is vertex id, second is the degree.
    std::map<int,int> outDegreeMap;
    std::map<int,int> inDegreeMap;

    for (int i = 0;i < num_vertices;i++) {
        inDegreeAssignation[i] = 0;
        outDegreeAssignation[i] = 0;
        saved_inDegreeAssignation[i] = 0;
        saved_outDegreeAssignation[i] = 0;
        synthetic_inDegreeAssignation[i] = 0;
        synthetic_outDegreeAssignation[i] = 0;
    }

    for (int i = 0;i < num_vertices;i++) {
        if((*this->g)[i].location_id == 0) {
            Vertex v = vertex(i,*this->g);
            int in_index = in_degree(v,*this->g);
            int out_index = out_degree(v,*this->g);
            inDegreeAssignation[in_index]--;                   //subtract the degree because DC0 is retained.
            outDegreeAssignation[out_index]--;
            inDegreeAssignation_vertex[in_index].push_back(i);
            outDegreeAssignation_vertex[out_index].push_back(i);
            outDegreeMap[i] = out_index;
            inDegreeMap[i] = in_index; // Keep the in_degree and out_degree in dc0 the same
        }
    }

    vertex_iter vi,vi_end;
    for (std::tie(vi,vi_end) = vertices(*this->g);vi != vi_end;vi++) {
        inDegreeAssignation[in_degree(*vi,*this->g)]++;
        outDegreeAssignation[out_degree(*vi,*this->g)]++;
        saved_inDegreeAssignation[in_degree(*vi,*this->g)]++;
        saved_outDegreeAssignation[out_degree(*vi,*this->g)]++;
    }
    edge_iter ei, ei_end, next,source_next;
    std::tie(ei,ei_end) = edges(*this->g);
    for (next = ei;next != ei_end; ei = next) {
        ++next;
        if ((*this->g)[source(*ei,*this->g)].location_id == (*this->g)[target(*ei,*this->g)].location_id
            && (*this->g)[source(*ei,*this->g)].location_id == 0) { // if the edges are in the same dc
            // Only deal with not 0 case, so the edges is not in the primary dc (0th).
        } else {
//            std::cout << "Source :" << (*this->g)[source(*ei,*this->g)].location_id << " target: " << (*this->g)[target(*ei,*this->g)].location_id << std::endl;
            remove_edge(*ei, *this->g);
            removedEdgeNum++;
        }
    }
    std::cout << "Jerry : " << removedEdgeNum << " edges removed, " << ecount - removedEdgeNum << " edges left" << std::endl;
	std::vector<int> permutation = utils::generateRandomPermutation(vertex_num);      //permutation is a random vector which size from 1 to n.
	
	/*
	std::cout << "Permutation begin:" << std::endl;
	for(int j=0;j<permutation.size();j++)
	std::cout << permutation[j] << std::endl;
	std::cout << "Permutation end:" << std::endl;
    */
	// Assign the dk1 series
    for (auto i : permutation) {                   //i equal to the element in permutation vector in turn.
        i = i - 1;
        if((*this->g)[i].location_id != 0) {
            Vertex v = vertex(i,*this->g);
            int in_index = 0;
            int out_index = 0;
            while (inDegreeAssignation[in_index] == 0) {
                in_index++;
            }
            while (outDegreeAssignation[out_index] == 0) {
                out_index++;
            }
            // TODO: low effiency here.
            inDegreeAssignation[in_index]--;
            outDegreeAssignation[out_index]--;
            inDegreeAssignation_vertex[in_index].push_back(i);
            outDegreeAssignation_vertex[out_index].push_back(i);
            outDegreeMap[i] = out_index - out_degree(vertex(i,*this->g),*this->g);
            inDegreeMap[i] = in_index - in_degree(vertex(i,*this->g),*this->g);
        }
    }
//    // Assignation done, Wire the edges.
    int count = 0;
    int self_loop_count = 0;
	
    for (int i = 0;i < num_vertices;i++) {
        inDegreeMap[i] -= in_degree(vertex(i,*this->g),*this->g);
        outDegreeMap[i] -= out_degree(vertex(i,*this->g),*this->g);
    }

    int sum = 0;
    for (auto out_pair : outDegreeMap) {
        sum += out_pair.second;
    }
    std::cout << "Jerry: total out degree = " << sum << std::endl;
    count = 0;
    int failed_count = 0;
// New Wiring method
	srand(time(NULL));    //Set random number seed, make each generation of random sequences different
    auto out_it = outDegreeMap.begin();
    auto out_end = outDegreeMap.end();
	auto in_it = inDegreeMap.begin();
    auto in_end = inDegreeMap.end();
	
	/*
	//inorder to limit the character of the power-law graph,we cut the degree of the vertex which's degree is too big
	for (; in_it != inDegreeMap.end(); in_it++)
	{
		if (in_it->second > 100)
			in_it->second = in_it->second * CUTCoefficient_0_1;
		//else if(in_it->second > 80)
			//in_it->second = in_it->second * CUTCoefficient_0_3;
		//else if (in_it->second > 50)
			//in_it->second = in_it->second * CUTCoefficient_0_5;
		//else if (in_it->second > 30)
			//in_it->second = in_it->second * CUTCoefficient_0_7;
		//else if (in_it->second > 20)
			//in_it->second = in_it->second * CUTCoefficient_0_9;
	}
	in_it = inDegreeMap.begin();
	for (; out_it != outDegreeMap.end(); out_it++)
	{
		if (out_it->second > 100)
			out_it->second = out_it->second * CUTCoefficient_0_1;
		
		//else if (out_it->second > 80)
			//out_it->second = out_it->second * CUTCoefficient_0_3;
		//else if (out_it->second > 50)
			//out_it->second = out_it->second * CUTCoefficient_0_5;
		//else if (out_it->second > 30)
			//out_it->second = out_it->second * CUTCoefficient_0_7;
		//else if (out_it->second > 20)
			//out_it->second = out_it->second * CUTCoefficient_0_9;
			
	}
	out_it = outDegreeMap.begin();

	*/

	/*
	*Save the backup data to restore the original state of the diagram 
	*if the composition fails 
	*/
	map<int,int> outDegreeMap_backup, inDegreeMap_backup;
	outDegreeMap_backup.clear();
	for (int i = 0; i < outDegreeMap.size(); i++)
		outDegreeMap_backup[i] = outDegreeMap[i];
	inDegreeMap_backup.clear();
	for (int i = 0; i < inDegreeMap.size(); i++)
		inDegreeMap_backup[i] = inDegreeMap[i];
	temp_g->clear();
	*this->temp_g = *this->g;
	  

	int num_edges = ecount;
	int generate_1k_graph_count=0;
	int numOfSubGraph = 0;

	while ( numOfSubGraph != 1)     //if the synthetic graph is not connected,generate another synthetic graph
	{


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


				if (p > random)    //Run the following code with probability p
				{
					////The adding edge is already exit if return false.
					bool success = add_edge(out_it->first,in_it->first, *this->g).second;
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
						//break;
					}
					else
					{
						//std::cout << "The edge is already exit." << std::endl;
						out_it++;
						if (out_it == out_end)
							break;
					}
				}
				in_it++;
				//if (in_it == in_end)
					//in_it = inDegreeMap.begin();
			}

		}
		generate_1k_graph_count++;
		numOfSubGraph = countSubGraph.operator()(*this->g);
	}
	
	std::cout << "The number of subgraph of PI distribution graph is:" << countSubGraph.operator()(*this->g) << std::endl;
	
	//Rewiring edges make the number of triangles the same.
	Triangles* CT=new Triangles();
	int Original_TRIANGLE = 0,Synthetic_TRIANGLE=0;
	Graph1::Graph gg;
	gg._g = *this->original_g; 
	Original_TRIANGLE = CT->countTriangles(gg);
	gg._g = *this->g;
	Synthetic_TRIANGLE = CT->countTriangles(gg);
	std::cout << "The number of triangles of original graph is:" << Original_TRIANGLE <<std::endl;
	std::cout << "The number of triangles of synthetic graph is:" << Synthetic_TRIANGLE << std::endl;
	std::cout << "Now begin to rewire the edges" << std::endl;

	int rewireEdge = 0;
	//create at 9_24 by Ribo to test the time consumer line.
	int executeTimes = 0,numSubgraph=0;
	size_t randomVertex = utils::generateRandomNumberFrom0TOmax(num_vertices);
	findRandomVertexNeighbor  FRVN;
	size_t target_neighbor = 0, source_neighbor = 0, fourthVertex=0,secondVertex=0; 
	bool little = false; //bigger = false;
	float random,p;
	//while (Synthetic_TRIANGLE < Synthetic_TRIANGLE*1.2 && Synthetic_TRIANGLE< Original_TRIANGLE*0.86) {
	//while (Synthetic_TRIANGLE < Original_TRIANGLE || Synthetic_TRIANGLE > Original_TRIANGLE) {
		while (Synthetic_TRIANGLE < 0) {     //choose this line to represent the number of triangles not reconstructed

		if (Synthetic_TRIANGLE < Original_TRIANGLE)
			little = true;
		else
			little = false;

		//select a random vertex from the 1k-series graph.
		randomVertex = utils::generateRandomNumberFrom0TOmax(num_vertices);
		target_neighbor = randomVertex;

		
		//It will continue only if the selected vertex has a in-degree greater than 2
		if (in_degree(target_neighbor, *this->g) >= 2)
		{
			//If the neighbor of neighbor of a random vertex has an outdegree greater than 2,
			//rewire the edge.
			std::tie(ei, ei_end) = edges(*this->g);
			for (next = ei; next != ei_end; next++) {
				//find the target vertex,and find the vertex's random out edge
				if (source(*next, *this->g) == target_neighbor) {
					target_neighbor = FRVN.FindRandomVertexNeighbor("out", target_neighbor, *this->g);
					break;
				}
			}
			secondVertex = target_neighbor;
			target_neighbor = FRVN.FindRandomVertexNeighbor("out", target_neighbor, *this->g);
			//find the random vertex's source vertex
			for (source_next = ei; source_next != ei_end; source_next++) {
				//find the target vertex
				if (target(*source_next, *this->g) == randomVertex)
				{
					source_neighbor = FRVN.FindRandomVertexNeighbor("in", randomVertex, *this->g);
					break;
				}
			}

			//The probability P determines whether the node is adopted or not
		/*
		*choose a vertex with probability P.
		*/
		//random = rand() % (N + 1) / (float)(N + 1);
		//p = (float)((float)(FRVN.countVertexOutDegree(randomVertex, *this->g)) / (float)ecount);
		//if (random < p)
		{
			//if the vertex's outdegree is greater than 2,rewire
			if (little == true && out_degree(target_neighbor, *this->g) >= 2)
			{

				//find the random next vertex and delete the edge
				fourthVertex = FRVN.FindRandomVertexNeighbor("out", target_neighbor, *this->g);

				//remove two edge
				remove_edge(target_neighbor, fourthVertex, *this->g);
				remove_edge(source_neighbor, randomVertex, *this->g);

				//add two edge
				bool addSucess_target = add_edge(target_neighbor, randomVertex, *this->g).second;
				bool addSucess_source = add_edge(source_neighbor, fourthVertex, *this->g).second;

				//if the number of triangles are more than the original graph,then save the rewiring
				//else undo the rewiring.
				int temp;
				gg._g = *this->g;
				temp = CT->countTriangles(gg);
				//std::cout << "After Rewire:" << temp << std::endl;
				//std::cout << "Before Rewire:" << Synthetic_TRIANGLE << std::endl;
				numSubgraph = countSubGraph.operator()(*this->g);
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
					//remove two edge
					remove_edge(target_neighbor, randomVertex, *this->g);
					remove_edge(source_neighbor, fourthVertex, *this->g);

					//add two edge
					add_edge(target_neighbor, fourthVertex, *this->g);
					add_edge(source_neighbor, randomVertex, *this->g);
					executeTimes++;
				}
			}
		}
				/*
				*If the selected vertex is a vertex in the triangle
				*then broke the triangle structural
				*/
				//random = rand() % (N + 1) / (float)(N + 1);
		       // p = (float)((float)(FRVN.countVertexOutDegree(randomVertex, *this->g)) / (float)ecount);
				//if (random < p)
				{
					if (little == false)
						//deal with the situation of more
					{
						//std::cout << "Begin to reduce the edges in the 1k-graph." << std::endl;
						//find the random next vertex and delete the edge
						fourthVertex = FRVN.FindRandomVertexNeighbor("in", target_neighbor, *this->g);

						//remove two random edges
						remove_edge(randomVertex, secondVertex, *this->g);
						remove_edge(fourthVertex, target_neighbor, *this->g);

						//add two edges
						bool addSucess_target = add_edge(randomVertex, target_neighbor, *this->g).second;
						bool addSucess_source = add_edge(fourthVertex, secondVertex, *this->g).second;

						//if the number of triangles are more than the original graph,then save the rewiring
						//else undo the rewiring.
						int temp;
						gg._g = *this->g;
						temp = CT->countTriangles(gg);
						//std::cout << "After Rewire:" << temp << std::endl;
						//std::cout << "Before Rewire:" << Synthetic_TRIANGLE << std::endl;
						numSubgraph = countSubGraph.operator()(*this->g);
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
			
		}
	}
	std::cout << "The number of triangles of synthetic graph is:" << CT->countTriangles(gg) << std::endl;

    // Rewiring self loop
    // TODO: bug: only one self loop?
    std::tie(ei,ei_end) = edges(*this->g);
    std::vector<Edge> self_loop_edges;
    for (next = ei;next != ei_end;next++) {
        if (source(*next,*this->g) == target(*next,*this->g)) {
            self_loop_edges.push_back(*next);
        }
    }
    if (self_loop_edges.size() == 1) {
        std::cout << " ERROR: Only One self loop edge in the graph" << std::endl;
    }
    if (!self_loop_edges.empty()) {
        int index1 = 0;
        int index2 = 1;
        while (index2 < self_loop_edges.size()) {
            add_edge(source(self_loop_edges[index1],*this->g),source(self_loop_edges[index2],*this->g),*this->g);
            count++;
            index1++;
            index2++;
        }
        add_edge(source(self_loop_edges[index2 - 1],*this->g),source(self_loop_edges[0],*this->g),*this->g);
        count++;

        for (index1 = 0;index1 < self_loop_edges.size();index1++) {
            remove_edge(self_loop_edges[index1],*this->g);
            count--;
        }
    }

    // Print out synthetic graph 1k-series
    std::cout << "Synthetic graph in degree distribution: " << std::endl;
    int *synthetic_in = utils::getIn1kSeries(this->g);
    for (int i = 0;i < num_vertices + 1;i++) {
        if (synthetic_in[i] != 0) {
            cout << i << "," << synthetic_in[i] << std::endl;
        }
    }
    std::cout << std::endl << "Synthetic graph out degree distribution: " << std::endl;
    int *synthetic_out = utils::getIn1kSeries(this->g);
    for (int i = 0;i < num_vertices + 1;i++) {
        if (synthetic_in[i] != 0) {
            cout << i << "," << synthetic_out[i] << std::endl;
        }
    }
    std::cout << count << " edges added. Self loop num: " << self_loop_count << std::endl;
    std::ofstream save("synthetic_g.txt");
    std::tie(ei,ei_end) = edges(*this->g);
    for (next = ei;next != ei_end;next++) {
        save << source(*next,*this->g) << " " << target(*next,*this->g) << std::endl;
    }

    save.close();
    graph_file.close();
    this->num_edges = ecount;
    std::cout << "load graph from file and split into 1k-series done.\n";
}
//end of load_from_file_1k_series









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
