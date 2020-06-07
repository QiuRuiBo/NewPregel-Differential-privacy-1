#ifndef APPS_H
#define APPS_H

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <limits>
#include "Cloud_instance.h"
#include "Simulator.h"
using namespace boost;

#endif

extern int Max_Monte_Carlo;
const float barrier_time = 0.03f; //randomly set to 3, need to verify
//const float pagerank_delta = 0.01; //can be changed, set arbitrarily
//const int pruned_DCs = 5;

enum AppType{
	pagerank,
	sssp,
	subgraph,
	typesofapps
};

class BaseApp{
public:

	AppType mytype;
	distributed_graph* global_graph;
	
	std::vector<int> sources;
	// for sssp 
	bool directed = true;
	// gather from both direction
	bool BiDirected = true;
	
	Graph* pattern_graph;
	/**
	* the rate of reducing message size
	* assume all vertices have the same value
	*/
	float red_rate;
	float budget;
	int ITERATIONS;
	/** msg size of each vertex in MB
	* assume all vertices have the same msg size
	*/
	float msg_size;
	//virtual float gather_msg_size(int vindex, distributed_graph*)=0;
	//virtual float apply_msg_size(int vindex, distributed_graph*)=0;
	
	virtual std::vector<int> get_gather_nbrs(int vindex,distributed_graph*)=0;	
	virtual std::vector<int> get_scatter_nbrs(int vindex,distributed_graph*)=0;	
	
	/**
	* application implementation using GAS model 
	*/
	virtual void init(distributed_graph*) = 0;
	virtual VertexData* gather(int,int,distributed_graph*)=0;
	virtual VertexData* sum(VertexData*,VertexData*)=0;
	virtual void apply(int,VertexData*,distributed_graph*)=0;
	virtual float scatter(int,int,distributed_graph*)=0;
	virtual void output(distributed_graph*)=0;
	//only implemented for PageRank and the Pregel model
        virtual void Compute(int,distributed_graph*)=0;
};


class PageRank:public BaseApp{
public:
	// Global random reset probability
	double RESET_PROB = 0.15;

	double TOLERANCE = 1.0E-2;

	bool USE_DELTA_CACHE = false;
	
	std::vector<int> get_gather_nbrs(int vindex,distributed_graph* dag){
		/**
		* gather neighbor: in_edges
		*/
		std::vector<int> gather_nbrs = dag->get_in_nbrs(vindex);
		return gather_nbrs;
	}
	std::vector<int> get_scatter_nbrs(int vindex,distributed_graph* dag){
		/**
		* scatter neighbor: out_edges
		*/		
		std::vector<int> scatter_nbrs = dag->get_out_nbrs(vindex);
		return scatter_nbrs;
	}

	void init(distributed_graph* dag){
		/*initial, activate all*/
		for (int viter = 0; viter < dag->g->m_vertices.size(); viter++){
			(*dag->g)[viter].status = activated;
			(*dag->g)[viter].next_status = deactivated;
			dynamic_cast<PageRankVertexData*>((*dag->g)[viter].data)->rank = 0;
		}
	}
	/* Gather the weighted rank of the adjacent page   */
	VertexData* gather(int vid, int nbr_id, distributed_graph* dag){
		int num_out_nbrs = dag->get_out_nbrs(nbr_id).size();
		if(num_out_nbrs == 0){
			std::cout << "gather error. sth wrong with local graph edges." <<std::endl;
			exit(1);
		}
		PageRankVertexData* result=new PageRankVertexData();
		result->rank = dynamic_cast<PageRankVertexData*>((*dag->g)[nbr_id].data)->rank / (float)num_out_nbrs;
		return result;
	}
	/*Sum the results from multiple gathering neighbors*/
	VertexData* sum(VertexData* a, VertexData* b){
		PageRankVertexData* result = new PageRankVertexData();
		result->rank = dynamic_cast<PageRankVertexData*>(a)->rank + dynamic_cast<PageRankVertexData*>(b)->rank;
		return result;
	}
	/* Use the total rank of adjacent pages to update this page */
	void apply(int vid, VertexData* total,distributed_graph* dag){
		float newvalue = (1.0 - RESET_PROB) * dynamic_cast<PageRankVertexData*>(total)->rank + RESET_PROB;
		float last_change = (newvalue - dynamic_cast<PageRankVertexData*>((*dag->g)[vid].data)->rank)/(float)dag->get_out_nbrs(vid).size();
		dynamic_cast<PageRankVertexData*>((*dag->g)[vid].data)->last_change = last_change;
		dynamic_cast<PageRankVertexData*>((*dag->g)[vid].data)->rank = newvalue;
		if(ITERATIONS) (*dag->g)[vid].next_status = activated;
	}
	/* The scatter function just signal adjacent pages */
	float scatter(int vid, int nbr_id,distributed_graph* dag) {
		// if(USE_DELTA_CACHE) {
		  // context.post_delta(edge.target(), last_change);
		// }
		if(USE_DELTA_CACHE){
			if(dynamic_cast<PageRankVertexData*>((*dag->g)[vid].data)->last_change > TOLERANCE || dynamic_cast<PageRankVertexData*>((*dag->g)[vid].data)->last_change < -TOLERANCE) {
				(*dag->g)[nbr_id].next_status = activated;
			} 				
		}else{
			(*dag->g)[nbr_id].next_status = activated;
		}	
		return dynamic_cast<PageRankVertexData*>((*dag->g)[vid].data)->last_change;
	  }
	  
	void output(distributed_graph* dag){
		for(int v=0; v<dag->g->m_vertices.size(); v++){
			std::cout << (*dag->g)[v].vertex_id << ": rank " <<dynamic_cast<PageRankVertexData*>((*dag->g)[v].data)->rank << std::endl;
		}			
	}
	/* float gather_msg_size(int vindex, distributed_graph* dag){
		float msg_size = 0.000008;		
		return msg_size;
	}
	float apply_msg_size(int vindex, distributed_graph* dag){
		float msg_size = 0.000008;		
		return msg_size;
	} */

	/* The message passing model in Pregel */
	void Compute(int vindex, distributed_graph* dag){
		(*dag->g)[vindex].status = deactivated;
	    int numofmsg = (*dag->g)[vindex].messages.size();
	    PageRankVertexData* sum = new PageRankVertexData();
        sum->rank = 0;
        for (int m=0; m < numofmsg; m++){
			//printf("当前是vertex %d,rank=%f\n", vindex, dynamic_cast<PageRankVertexData*>((*dag->g)[vindex].messages[m])->rank);
            sum->rank += dynamic_cast<PageRankVertexData*>((*dag->g)[vindex].messages[m])->rank;
        }
		float newvalue = RESET_PROB + (1-RESET_PROB) * sum->rank;
		float last_change = (newvalue - dynamic_cast<PageRankVertexData*>((*dag->g)[vindex].data)->rank)/(float)dag->get_out_nbrs(vindex).size();
		dynamic_cast<PageRankVertexData*>((*dag->g)[vindex].data)->last_change = last_change;
        dynamic_cast<PageRankVertexData*>((*dag->g)[vindex].data)->rank = newvalue;
		//clear message queue

		for(int m=0; m<numofmsg; m++){
			//delete (*dag->g)[vindex].messages[m];  
		}
		(*dag->g)[vindex].messages.clear();

		//printf("msg_size=%d\n", (*dag->g)[vindex].messages.size());

		//if(ITERATIONS || last_change > TOLERANCE || last_change < -TOLERANCE)
		if (last_change > TOLERANCE || last_change < -TOLERANCE)
		{
			//printf("activate vertex %d\n", vindex);
			(*dag->g)[vindex].status = activated;
		}
	}
	PageRank(){
		red_rate = 0.5; BiDirected = false;
		msg_size = 0.000008f;
	}
};

class SSSP :public BaseApp{
public:
	double TOLERANCE = 500;
	std::vector<int> get_gather_nbrs(int vindex,distributed_graph* dag){
		/**
		* gather neighbor: all neighbors if undirected
		*				   in edges if directed
		*/
		std::vector<int> gather_nbrs = dag->get_in_nbrs(vindex);
		if (!directed){
			std::vector<int> out_nbrs = dag->get_out_nbrs(vindex);
			for (int i = 0; i < out_nbrs.size(); i++)
				gather_nbrs.push_back(out_nbrs[i]);
		}
		return gather_nbrs;
	}
	std::vector<int> get_scatter_nbrs(int vindex,distributed_graph* dag){
		/**
		* scatter neighbor: all_edges if undirected; out edges if directed
		*/
		std::vector<int> scatter_nbrs = dag->get_out_nbrs(vindex);
		if (!directed){
			std::vector<int> in_nbrs = dag->get_in_nbrs(vindex);
			for (int i = 0; i < in_nbrs.size(); i++)
				scatter_nbrs.push_back(in_nbrs[i]);
		}		
		return scatter_nbrs;
	}

	void init(distributed_graph* dag){
		/*initial, activate all sources*/
		for (int viter = 0; viter < dag->g->m_vertices.size(); viter++){
			dynamic_cast<SSSPVertexData*>((*dag->g)[viter].data)->dist = std::numeric_limits<float>::max();
			for (int siter = 0; siter < sources.size(); siter++){
				/* global to local vertex id */
				if((*dag->g)[viter].vertex_id == sources[siter]){
					(*dag->g)[viter].status = activated;
					break;
				}
			}
		}
	}

	inline int get_other_vertex(int edge_id, int v_id, distributed_graph* dag){
		// edge_iter it, end;
		// std::tie(it, end) = edges(*dag->g);
		// for(int i=edge_id; i>0; i--)
			// ++it;		
		// return source(*it,*dag->g) == v_id ? target(*it,*dag->g) : source(*it,*dag->g);
		int src = source(dag->random_access_edges[edge_id], *dag->g);
		int tgt = target(dag->random_access_edges[edge_id], *dag->g);
		return src == v_id ? tgt : src;
	}
	
	/*Collect the distance to the neighbor*/
	VertexData* gather(int vid, int nbr_id, distributed_graph* dag){
		//std::cout << "gather not implemented for sssp" << std::endl;
		in_edge_iterator in_i, in_end;
		Edge e, edge_id;
		boost::tie(in_i, in_end) = in_edges(vid, *dag->g);
		for (; in_i != in_end; ++in_i) {
			e = *in_i;
			Vertex src = source(e, *dag->g);
			if(src == nbr_id){
				edge_id = e;
				break;
			}
		}
		SSSPVertexData* result = new SSSPVertexData();
		result->dist = (dynamic_cast<SSSPVertexData*>((*dag->g)[nbr_id].data)->dist + dag->WeightMap[edge_id]);
		return result;		
	}
	/*Sum the results from multiple gathering neighbors*/
	VertexData* sum(VertexData* a, VertexData* b){
		SSSPVertexData* result = new SSSPVertexData();
		result->dist = std::min(dynamic_cast<SSSPVertexData*>(a)->dist,dynamic_cast<SSSPVertexData*>(b)->dist);
		return result;
	}
	/*If the distance is smaller then update*/
	void apply(int vid, VertexData* acc,distributed_graph* dag){
		dynamic_cast<SSSPVertexData*>((*dag->g)[vid].data)->changed = false;
		if(dynamic_cast<SSSPVertexData*>((*dag->g)[vid].data)->dist > dynamic_cast<SSSPVertexData*>(acc)->dist) {
		  dynamic_cast<SSSPVertexData*>((*dag->g)[vid].data)->changed = true;
		  dynamic_cast<SSSPVertexData*>((*dag->g)[vid].data)->dist = dynamic_cast<SSSPVertexData*>(acc)->dist;
		}
	}
	/*The scatter function just signal adjacent pages */
	float scatter(int vid, int nbr_id,distributed_graph* dag){
		// int other = get_other_vertex(edge_id, vid, dag);
		// edge_iter e_it, e_end;
		// int count = 0;
		// for(std::tie(e_it, e_end) = edges(*dag->g); count<edge_id; ++e_it)
			// count++;
		out_edge_iterator out_i, out_end;
		Edge e,edge_id;
		boost::tie(out_i, out_end) = out_edges(vid, *dag->g);
		for (; out_i != out_end; ++out_i){
			e = *out_i;
			Vertex tgt = target(e, *dag->g);
			if(tgt == nbr_id){
				edge_id = e;
				break;
			}
		}

		float newd = dynamic_cast<SSSPVertexData*>((*dag->g)[vid].data)->dist + dag->WeightMap[edge_id];
		if(dynamic_cast<SSSPVertexData*>((*dag->g)[vid].data)->changed)
			(*dag->g)[nbr_id].next_status = activated;
		return newd;
	}
	void output(distributed_graph* dag){
		std::cout << "output not implemented for sssp." <<std::endl;
	}

	/* The message passing model in Pregel */
	void Compute(int vindex, distributed_graph* dag) {
		int numofmsg = (*dag->g)[vindex].messages.size();

		//std::cout << "sources = " << dynamic_cast<SSSPVertexData*>((*dag->g)[0].data)->dist;

		SSSPVertexData* min = new SSSPVertexData();
		dynamic_cast<SSSPVertexData*>((*dag->g)[vindex].data)->changed = false;
		float newvalue = dynamic_cast<SSSPVertexData*>((*dag->g)[vindex].data)->dist;
		for (int m = 0; m < numofmsg; m++) {
			if (newvalue > dynamic_cast<SSSPVertexData*>((*dag->g)[vindex].messages[m])->dist)
			{
				//std::cout << "这是更新操作" << std::endl;
				dynamic_cast<SSSPVertexData*>((*dag->g)[vindex].data)->changed = true;
				newvalue = dynamic_cast<SSSPVertexData*>((*dag->g)[vindex].messages[m])->dist;
			}
		}
		dynamic_cast<SSSPVertexData*>((*dag->g)[vindex].data)->last_change = abs(newvalue - dynamic_cast<SSSPVertexData*>((*dag->g)[vindex].data)->dist);
		float last_change = dynamic_cast<SSSPVertexData*>((*dag->g)[vindex].data)->last_change;
		dynamic_cast<SSSPVertexData*>((*dag->g)[vindex].data)->dist = newvalue;
		//clear message queue
		for (int m = 0; m < numofmsg; m++) {
			//delete (*dag->g)[vindex].messages[m];  
		}
		(*dag->g)[vindex].messages.clear();

		if (ITERATIONS || last_change > TOLERANCE || last_change < -TOLERANCE)
			(*dag->g)[vindex].status = activated;
	}
	SSSP(){ 
		red_rate=0.5f; msg_size = 0.000004f; 
		if(directed) BiDirected = false;
		else BiDirected = true;
	}
};


class SubgraphIsom: public BaseApp{
	
public:
	std::vector<int> get_gather_nbrs(int vindex,distributed_graph* dag){
		/**
		* gather neighbor: all neighbors 
		*/
		std::vector<int> gather_nbrs = dag->get_in_nbrs(vindex);
		std::vector<int> out_nbrs = dag->get_out_nbrs(vindex);
		for (int i = 0; i < out_nbrs.size(); i++)
			gather_nbrs.push_back(out_nbrs[i]);		
		return gather_nbrs;
	}
	std::vector<int> get_scatter_nbrs(int vindex,distributed_graph* dag){
		/**
		* scatter neighbor: all neighbors
		*/
		std::vector<int> scatter_nbrs = dag->get_out_nbrs(vindex);
		std::vector<int> in_nbrs = dag->get_in_nbrs(vindex);
		for (int i = 0; i < in_nbrs.size(); i++)
			scatter_nbrs.push_back(in_nbrs[i]);		
		return scatter_nbrs;
	}

	void init(distributed_graph* dag){
		/*initial, activate all*/
		//#pragma omp parallel for
		for (int viter = 0; viter < dag->g->m_vertices.size(); viter++){
			(*dag->g)[viter].status = activated;
			(*dag->g)[viter].next_status = deactivated;
			int gvid = (*dag->g)[viter].vertex_id;			
			for(int i=0; i<pattern_graph->m_vertices.size(); i++){
				//no label, so add all
				SubgraphVertexData::Message init_msg;
				init_msg.pairs.push_back(std::pair<int,int> (i,gvid));	
				init_msg.forwarding_trace.push_back(gvid);
				dynamic_cast<SubgraphVertexData*>((*dag->g)[viter].data)->matches.push_back(init_msg);				
			}			
		}
	}

	/*return the messages of neighbors*/
	VertexData* gather(int vid, int nbr_id, distributed_graph* dag){			
		SubgraphVertexData* result = new SubgraphVertexData();
		for(int i=0; i< dynamic_cast<SubgraphVertexData*>((*dag->g)[nbr_id].data)->matches.size(); i++)
			result->matches.push_back(dynamic_cast<SubgraphVertexData*>((*dag->g)[nbr_id].data)->matches[i]);
		return result;
	}
	
	/*Sum the results from multiple gathering neighbors*/
	VertexData* sum(VertexData* a, VertexData* b){
		SubgraphVertexData* result = new SubgraphVertexData();
		int asize = dynamic_cast<SubgraphVertexData*>(a)->matches.size();
		int bsize = dynamic_cast<SubgraphVertexData*>(b)->matches.size();
		for(int i=0; i<asize; i++)
			result->matches.push_back(dynamic_cast<SubgraphVertexData*>(a)->matches[i]);
		for(int i=0; i<bsize; i++)
			result->matches.push_back(dynamic_cast<SubgraphVertexData*>(b)->matches[i]);
		return result;
	}
	
	/** On receiving messages, vertex vg uses the following operations on each message according to different situations.
	* Forward: if vg is not in the forwarding trace, 
			   add it to the forwarding trace, and send the message to all neighbors.
	* Update (Spawn): If vg has not matched any vertex in the pattern and it can match to
		some neighbor of the vertex in the pattern, which the message source has matched, then
		for each possible matching vertex vp, vertex function generates a new message based on
		the old one by adding the matching pair and replacing all IDs in the forwarding trace
		with the ID of the current vertex vg. In this case, a message may spawn multiple new
		messages. If all the pattern vertices have been matched already after this operation in a
		new message, then an isomorphic subgraph has been found. Otherwise the messages are sent to all neighbors.
	* Drop: In other situations, the message is dropped, because it is hopeless or redundant to
		find further matching vertex.
	*/
	void apply(int vid, VertexData* acc, distributed_graph* dag){
		dynamic_cast<SubgraphVertexData*>((*dag->g)[vid].data)->matches.clear();
		//use global vid 
		int gvid = (*dag->g)[vid].vertex_id;
		if(false){
		//if(dag->init_indicator){
			/* VertexData::Message init_msg;
			for(int i=0; i<pattern_graph->m_vertices.size(); i++){
				//no label, so add all
				init_msg.pairs.push_back(std::pair<int,int> (i,gvid));		
			}
			init_msg.forwarding_trace.push_back(gvid);
			(*dag->g)[vid].data->matches.push_back(init_msg); */
			dag->init_indicator = false;
		}else{
			int num_msg = dynamic_cast<SubgraphVertexData*>(acc)->matches.size();
			for(int i=0; i<num_msg; i++){
				int j = 0;
				for(; j< dynamic_cast<SubgraphVertexData*>(acc)->matches[i].pairs.size(); j++){
					if(dynamic_cast<SubgraphVertexData*>(acc)->matches[i].pairs[j].second == gvid){
						if(std::find(dynamic_cast<SubgraphVertexData*>(acc)->matches[i].forwarding_trace.begin(),
							dynamic_cast<SubgraphVertexData*>(acc)->matches[i].forwarding_trace.end(),gvid) == dynamic_cast<SubgraphVertexData*>(acc)->matches[i].forwarding_trace.end()){
							//forward operation
							dynamic_cast<SubgraphVertexData*>(acc)->matches[i].forwarding_trace.push_back(gvid);
							dynamic_cast<SubgraphVertexData*>((*dag->g)[vid].data)->matches.push_back(dynamic_cast<SubgraphVertexData*>(acc)->matches[i]);
						}
						break;						
					}//end if				
				}//end for each match pair
				if(j == dynamic_cast<SubgraphVertexData*>(acc)->matches[i].pairs.size()){
					//update operation
					int last_trace = dynamic_cast<SubgraphVertexData*>(acc)->matches[i].forwarding_trace.back();
					int pattern_vertex = -1;
					for(j=0; j<dynamic_cast<SubgraphVertexData*>(acc)->matches[i].pairs.size(); j++){
						if(dynamic_cast<SubgraphVertexData*>(acc)->matches[i].pairs[j].second == last_trace){
							pattern_vertex = dynamic_cast<SubgraphVertexData*>(acc)->matches[i].pairs[j].first;
							break;							
						}
					}
					adja_iterator nbrIt, nbrEnd;
					boost::tie(nbrIt, nbrEnd) = adjacent_vertices( pattern_vertex, *pattern_graph );
					for (; nbrIt != nbrEnd; ++nbrIt){
						//check the matching criteria
						int k = 0;
						for(;k<dynamic_cast<SubgraphVertexData*>(acc)->matches[i].pairs.size(); k++){
							if(dynamic_cast<SubgraphVertexData*>(acc)->matches[i].pairs[k].first == *nbrIt)
								break;
						}
						if(k == dynamic_cast<SubgraphVertexData*>(acc)->matches[i].pairs.size()){
							//this nbr is not matched yet		
							bool nomatch = false;
							adja_iterator nnIt, nnEnd;									
							for(k = 0; k < dynamic_cast<SubgraphVertexData*>(acc)->matches[i].pairs.size(); k++){
								//if there is an edge between *nbrIt and the k-th matched pattern
								boost::tie(nnIt, nnEnd) = adjacent_vertices( *nbrIt, *pattern_graph );
								for(; nnIt != nnEnd; ++nnIt){
									if(*nnIt == dynamic_cast<SubgraphVertexData*>(acc)->matches[i].pairs[k].first){
										//check if the same edge exist in the graph
										adja_iterator localIt, localEnd;
										boost::tie(localIt, localEnd) = adjacent_vertices( dynamic_cast<SubgraphVertexData*>(acc)->matches[i].pairs[k].second, *global_graph->g );
										for(; localIt != localEnd; ++ localIt)
											if(*localIt == gvid)
												break;
										if(localIt == localEnd){
											nomatch = true;
											break;
										}
									}
								}
								if(nomatch)
									break;		
							}//for all pairs
							if(!nomatch){
								SubgraphVertexData::Message mprime;
								mprime.pairs = dynamic_cast<SubgraphVertexData*>(acc)->matches[i].pairs;
								mprime.pairs.push_back(std::pair<int,int>(*nbrIt,gvid));
								mprime.forwarding_trace.push_back(gvid);
								if(mprime.pairs.size() == pattern_graph->m_vertices.size()){
									//found a subgraph
									/* std::cout << "found a subgraph, with vertices: \n" ;
									for(int kk=0; kk<mprime.pairs.size(); kk++){
										std::cout << mprime.pairs[kk].first << ",\t" << mprime.pairs[kk].second << std::endl;
									} */
								}else{
									dynamic_cast<SubgraphVertexData*>(acc)->matches.push_back(mprime);
								}
							}
						}
					}//end for pv
				}
			}//end for each msg
		}		
	}
	
	/*The scatter function signal all neighbors */
	float scatter(int vid, int nbr_id, distributed_graph* dag){
		if(dynamic_cast<SubgraphVertexData*>((*dag->g)[vid].data)->matches.size() > 0){
			(*dag->g)[vid].next_status = activated;
			(*dag->g)[nbr_id].next_status = activated;	
		}				
		return 0.0;
	}
	void output(distributed_graph* dag){}
	/* The message passing model in Pregel */
    void Compute(int vindex, distributed_graph* dag){}

	SubgraphIsom(){ 
		red_rate=1.0; msg_size = 1.0; 
		BiDirected = true;
	} //msg size not constant
};
