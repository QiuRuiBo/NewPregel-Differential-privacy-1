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
				* Data aggregation
				*/
				std::cout << "-------------- Data aggregation stage  "<<myapp->ITERATIONS - iter_counter<<" -------------" << std::endl;
				for(int tr=0; tr<num_threads; tr++){
					for(int ai=0; ai<DCs[tr]->aggregators.size(); ai++){
						delete DCs[tr]->aggregators[ai];
					}
					DCs[tr]->aggregators.clear();
				}
                    for(int tr=0; tr<num_threads; tr++){
                        for(int v=0; v<Threads[tr]->l_dag->num_vertices; v++){
						int vid = (*Threads[tr]->l_dag->g)[v].vertex_id;
                        std::vector<int> out_nbrs = myapp->global_graph->get_out_nbrs(vid);
                        float msg_rank = dynamic_cast<PageRankVertexData*>((*myapp->global_graph->g)[vid].data)->rank/(float)out_nbrs.size();

						for(int nid=0; nid<out_nbrs.size(); nid++){
							VertexData* l_accum;
                            if(myapp->mytype == pagerank){
                                l_accum = new PageRankVertexData();
								dynamic_cast<PageRankVertexData*>(l_accum)->rank = msg_rank;
								int other_dc = (*myapp->global_graph->g)[out_nbrs[nid]].location_id;
                                    if(other_dc != tr){
										//search for the aggregator, if does not exist, create one
										bool foundagg = false;
										int founditer = -1;
										for(int ai=0; ai< DCs[tr]->aggregators.size(); ai++){
											if(DCs[tr]->aggregators[ai]->DC_id == other_dc && DCs[tr]->aggregators[ai]->vertex_id == out_nbrs[nid]){
												foundagg = true;
												founditer = ai;
												break;
									}
									if(foundagg){
										//aggregate to founditer
										dynamic_cast<PageRankVertexData*>(DCs[tr]->aggregators[founditer]->aggregated_data)->rank += msg_rank;
									}else{
										//create one
										Aggregator* new_agg = new Aggregator();
										new_agg->DC_id = other_dc;
										new_agg->vertex_id = out_nbrs[nid];
										new_agg->aggregated_data = l_accum;
										DCs[tr]->aggregators.push_back(new_agg);
									//	printf("create one aggregator in DC %d\n",tr);
									}
								}
                                   else {
									//simply send the message
									(*myapp->global_graph->g)[out_nbrs[nid]].messages.push_back(l_accum);
                                                                }
                                                        }else if(myapp->mytype == sssp){
                                                                l_accum = new SSSPVertexData();
                                                        }else if(myapp->mytype == subgraph){
                                                                l_accum = new SubgraphVertexData();
                                                        }
						}		
					}	
				}
				/*Barrier*/

				//After data aggregation, add noises to the aggregated data and send out
				int agg_count = 0;
				for(int tr=0; tr<num_threads; tr++){
					agg_count += DCs[tr]->aggregators.size();
				}
				printf("%d aggregators in total.\n",agg_count);
				
				for(int tr=0; tr<num_threads; tr++){
					/**
                                	* Imaging the largest diff neighboring graphs (containing only vertices in local DC)
                                	* G: a node k <- all the other nodes but v, and k -> v, a random node n -> v
                                	* G': edge k -> v is removed from G
                                	*/
                                	float deltaf = 0; //the same for all aggregators in the same DC
                                	float rank_hub = 0.15 + 0.85*(0.15*(Threads[tr]->l_dag->num_vertices - 3)+0.15/2.0);
                                	//the gateway node
                                	float fg= 0.15 + 0.85*(rank_hub+0.15/2.0);
                                	float fg1 = 0.15 + 0.85*0.15/2.0;
                                	if(iter_counter < myapp->ITERATIONS)
                                        	deltaf = std::abs(fg - fg1); 
                                	printf("Gf is: %f\n", deltaf);
                                
					for(int ai=0; ai < DCs[tr]->aggregators.size(); ai++){
						DCs[tr]->aggregators[ai]->noise_budget = noise_budget / (float)agg_count/(float)myapp->ITERATIONS;
						if(privacy)//generate noise value
							dynamic_cast<PageRankVertexData*>(DCs[tr]->aggregators[ai]->aggregated_data)->rank += laplace_generator(0,deltaf/DCs[tr]->aggregators[ai]->noise_budget);
						int other_vertex =DCs[tr]->aggregators[ai]->vertex_id; 
						//printf("aggregator in DC %d send message to vertex %d\n", tr,other_vertex);
						(*myapp->global_graph->g)[other_vertex].messages.push_back(DCs[tr]->aggregators[ai]->aggregated_data);
					}
				}

				iter_counter -- ;	
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
	
	for(int tr=0; tr<num_threads; tr++){
		printf("DC %d has %d aggregators\n",tr,DCs[tr]->aggregators.size());
	}
	//modified in 2019_7_5
	//return;
}
