#include "Simulator.h"
#include <omp.h>
#include <algorithm>

extern int N_THREADS;
extern float noise_budget;
extern bool privacy;
void GraphEngine::Simulate(){
	
	/**
	* Clear memory.
	*/
	Threads.clear();
	for(int i=0; i< num_threads; i++){
		Threads.push_back(new Thread());
	}
	/**
	* Build the local graphs for each thread
	* according to the replication plan
	* each vertex has a <local, global> ID
	*/
	float t_staging = 0.0;
	/**           
        * put all high-degree vertices in the same DC
        */
	int total_vertices = num_vertices(*myapp->global_graph->g);
	/*for(int vi=0; vi < total_vertices; vi++){
		std::vector<int> in_nbrs = myapp->global_graph->get_in_nbrs(vi);
		if(in_nbrs.size() >= total_vertices*0.003)
			(*myapp->global_graph->g)[vi].location_id = (Amazon_EC2_regions)0;
	}*/
	//debug
	omp_set_num_threads(N_THREADS);
	
	#pragma omp parallel for
        for(int i=0; i<num_threads; i++){
                if(num_threads != DCs.size()){
                        std::cout << "the number of threads has to be the same as the number of DCs" << std::endl;
                        exit(1);
                }
                printf("constructing local graph for thread %d\n",i);
                Threads[i]->DC_loc_id = i;
                Threads[i]->l_dag = new distributed_graph();
                Threads[i]->l_dag->g = new Graph();

                //add vertices
                int local_id = 0;
                for(int vi=0; vi < total_vertices; vi++){
                        if((*myapp->global_graph->g)[vi].location_id == i){
                                MyVertex* v = new MyVertex((*myapp->global_graph->g)[vi]);
                                v->local_vid = local_id;
                                v->is_master = true;
                                add_vertex(*v,*Threads[i]->l_dag->g);
                                //construct vid2lvid map
                                Threads[i]->vid2lvid.emplace(v->vertex_id,v->local_vid);
                                local_id ++;
                        }
                }
                Threads[i]->l_dag->num_vertices = num_vertices(*Threads[i]->l_dag->g);
                printf("Thread %d\n #vertices: %d\n",i,Threads[i]->l_dag->num_vertices);
        }	

	for(int tr=0; tr<num_threads; tr++){
                /**
                * Signal vertices
                */
                //myapp->init(Threads[tr]->l_dag);
                for(int v=0; v<Threads[tr]->l_dag->num_vertices; v++){
                        int vid = (*Threads[tr]->l_dag->g)[v].vertex_id;
                        dynamic_cast<PageRankVertexData*>((*myapp->global_graph->g)[vid].data)->rank = 0;
			//identify gateway nodes
			int distant_nbr = 0;
			std::vector<int> out_nbrs = myapp->global_graph->get_out_nbrs(vid);
	                for(int nid=0; nid < out_nbrs.size(); nid++){
                        	int other_dc = (*myapp->global_graph->g)[out_nbrs[nid]].location_id;
                                if(other_dc != tr){
                                	distant_nbr ++;
                                }
                        }
			(*myapp->global_graph->g)[vid].distant_nbr = distant_nbr;
                }
        }


	/** Start the engine, each thread executes at the same time */
	float t_step = 0.0;
	int iter_counter = 0;
	float totalcost = 0.0;
	float wan_usage = 0.0;
	/**
	* Start the engine
	* currently only implemented the synchronous engine 
	*/
	if(type == synchronous){
		if(myapp->ITERATIONS != 0){
			/* converge within the fixed number of iterations */
			iter_counter = myapp->ITERATIONS;
			while(iter_counter > 0){
				
				/**
				* Execute Compute on all active vertices
				* Sync before sending messages
				*/
				std::cout << "-------------- Compute stage of iteration "<<myapp->ITERATIONS - iter_counter<<" -------------" << std::endl;
				#pragma omp parallel for
                                for(int tr=0; tr<num_threads; tr++){
					#pragma omp parallel for
                                        for(int v=0; v<Threads[tr]->l_dag->num_vertices; v++){
						int vid = (*Threads[tr]->l_dag->g)[v].vertex_id;
						myapp->Compute(vid,myapp->global_graph);
					}	
				}
				/* BARRIER */				

				/**
				* Send messages to neighbors
				*/
				std::cout << "-------------- Message passing stage  "<<myapp->ITERATIONS - iter_counter<<" -------------" << std::endl;
				/////add locks: push_back to the same node at the same time
				int numoflocks = num_threads;
                                std::vector<omp_lock_t> writelock(numoflocks);
                                for(int lk=0; lk<numoflocks; lk++){
                                        omp_init_lock(&writelock[lk]);
                                }
				/******disabled omp to avoid the locking problem, improve later*******/
				// TODO: improve the locking problem
			
				//#pragma omp parallel for
                                for(int tr=0; tr<num_threads; tr++){
					/**
                                        * Imaging the largest diff neighboring graphs (containing only vertices in local DC)
                                        * G: a node k <- all the other nodes but v, and k -> v, a random node n -> v
                                        * G': edge k -> v is removed from G
                                        */
					float deltaf = 0; //the same for all gateway nodes in the same DC
					if(privacy){
						//starting from rank=0.15 for all vertices
						//the hub node
                                                float rank_hub = 0.15 + 0.85*(0.15*(Threads[tr]->l_dag->num_vertices - 3)+0.15/2.0);
						//the gateway node
						float fg= 0.15 + 0.85*(rank_hub+0.15/2.0);
						float fg1 = (float)(0.15 + 0.85 * 0.15 / 2.0);
						if(iter_counter < myapp->ITERATIONS)
							deltaf = std::abs(fg - fg1);
						printf("Gf is: %f\n", deltaf);
					}	
					//#pragma omp parallel for
                                        for(int v=0; v<Threads[tr]->l_dag->num_vertices; v++){
                                                int vid = (*Threads[tr]->l_dag->g)[v].vertex_id;
                                                std::vector<int> out_nbrs = myapp->global_graph->get_out_nbrs(vid);
						float msg_rank = dynamic_cast<PageRankVertexData*>((*myapp->global_graph->g)[vid].data)->rank/(float)out_nbrs.size();
						float noise_rank = msg_rank;
						/** 
                                                * 1. Gf=max diff between nbr graphs
                                                * 2. noise=Lap(Gf/noise_budget)
                                                */
						if(privacy){
							if((*myapp->global_graph->g)[vid].distant_nbr > 0){
								float epsilon = noise_budget/std::pow(2.0,(myapp->ITERATIONS-iter_counter+1))/(float)(*myapp->global_graph->g)[vid].distant_nbr;
								//printf("epsilon in iteration %d is %f\n", myapp->ITERATIONS-iter_counter, epsilon);
								noise_rank = msg_rank + laplace_generator(0,deltaf/epsilon);
								//printf("generate noise value %f\n", noise_rank-msg_rank);
							}
						}
						for(int nid=0; nid < out_nbrs.size(); nid++){
							VertexData* l_accum;
                                                        if(myapp->mytype == pagerank){
                                                                l_accum = new PageRankVertexData();
								if((*myapp->global_graph->g)[out_nbrs[nid]].location_id != tr){
									dynamic_cast<PageRankVertexData*>(l_accum)->rank = noise_rank;
								}
								else {
									dynamic_cast<PageRankVertexData*>(l_accum)->rank = msg_rank;
								}
                                                        }else if(myapp->mytype == sssp){
                                                                l_accum = new SSSPVertexData();
                                                        }else if(myapp->mytype == subgraph){
                                                                l_accum = new SubgraphVertexData();
                                                        }
							(*myapp->global_graph->g)[out_nbrs[nid]].messages.push_back(l_accum);
						}
                                        }
				}

				for(int lk=0; lk<numoflocks; lk++){
                                        omp_destroy_lock(&writelock[lk]);
                                }
				iter_counter -- ;	

				//log to check the convergence
			        for(int tr=0; tr<num_threads; tr++){
					for(int v=0; v<Threads[tr]->l_dag->num_vertices; v++){
			                        int vid = (*Threads[tr]->l_dag->g)[v].vertex_id;
						if((*myapp->global_graph->g)[vid].distant_nbr > 0)
				                        printf("rank value of gateway %d is %.4f\n",vid,dynamic_cast<PageRankVertexData*>((*myapp->global_graph->g)[vid].data)->rank);		                					
					}
			        }

				float avg_rank = 0;
                                for(int vi=0; vi<total_vertices; vi++){
                                        avg_rank += dynamic_cast<PageRankVertexData*>((*myapp->global_graph->g)[vi].data)->rank;
                                }
                                avg_rank /= (float)total_vertices;
                                printf("iter %d: average rank %.4f\n",myapp->ITERATIONS-iter_counter,avg_rank);

			}
			
		}else{
			/* converge according to the tolerance */
			
			iter_counter ++;
		}
	}else if(type == asynchronous){
		std::cout <<"asynchronous engine is not implemented. " << std::endl;
	}

	//output result
	//for(int tr=0; tr<num_threads; tr++){
	//	myapp->output(Threads[tr]->l_dag);
	//}
	
	return;
}
