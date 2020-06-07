#include "Simulator.h"
#include <omp.h>
/**
* the expected number of replicas for a vertex
* can be fixed by a user (we use 4 as an example)
* or calculated using methods in PowerGraph paper
*/
extern int N_THREADS;

// ***************************** Heterogeneity-aware Partitioning ***************************** //
bool mycompare(std::pair<int,float> a, std::pair<int,float> b){
	return a.second > b.second;
}

void Optimizer::StreamEdgeAssignment(){
	
	/**
	* The edge assignment rules:
	* 1) If R(v) and R(u) intersect, place edge (u,v) into one of
	  the intersected dcs r with the lowest g_r(i);
	  2) If R(v) and R(u) are not empty and do not intersect,
	  place (u,v) in a dc in R(v) if we have:
	  a_u(i)+g_m(i) < av_(i)+g_n(i); for any n in R(u) (8)
	  3) If only R(v) or R(u) is not empty, choose one dc
	  r from the non-empty set with the lowest g_r(i).
	  4) If both R(v) and R(u) are empty, place (u;v) in any
	  dc r with the lowest g_r(i).
	*/
	int vertexsize = dag->num_vertices;

	std::vector<int> orig_v = std::vector<int>(DCs.size(),0);
	for(int vi=0; vi<vertexsize; vi++){
		(*dag->g)[vi].replica_location.clear();
		(*dag->g)[vi].replica_location.push_back((*dag->g)[vi].location_id);
		orig_v[(*dag->g)[vi].location_id] ++;
	}
	for(int i=0; i<DCs.size(); i++){
		printf("datacenter %d: %d vertices\n",i,orig_v[i]);
	}
	int numedges = num_edges(*dag->g);
	edge_distribution = std::vector<std::vector<int> >(DCs.size());
	
	omp_set_num_threads(N_THREADS);
	int numoflocks = dag->num_vertices;
	std::vector<omp_lock_t> writelock(numoflocks);
	for(int i=0; i<numoflocks; i++){
		omp_init_lock(&writelock[i]);
	}
	int numoflocks1 = DCs.size();
	std::vector<omp_lock_t> writelock1(numoflocks1);
	for(int i=0; i<numoflocks1; i++){
		omp_init_lock(&writelock1[i]);
	}
#pragma omp parallel for
	for(int ecount=0; ecount<numedges; ecount++){
		Edge e = dag->random_access_edges[ecount];
		int src = source(e, *dag->g);
		int tgt = target(e, *dag->g);	
			
		int threadid = omp_get_thread_num();
		int lockid1, lockid2;
		lockid1 = src<tgt?std::fmod(src,(double)numoflocks):std::fmod(tgt,(double)numoflocks);
		lockid2 = src>tgt?std::fmod(src,(double)numoflocks):std::fmod(tgt,(double)numoflocks);
		// int lockid = std::fmod(pairing_function(src,tgt),(double)numoflocks);
		// printf("thread %d acquiring lock %d and %d: %d-->%d\n",threadid,lockid1,lockid2,src,tgt);
        omp_set_lock(&writelock[lockid1]);	
		omp_set_lock(&writelock[lockid2]);
		
		std::vector<int> notmoveto;
		std::pair<int,int> result = EdgeAssign(src,tgt,notmoveto,orig_v);
		int dc_to_assign = result.first;
		int tag = result.second;
		omp_set_lock(&writelock1[dc_to_assign]);
		edge_distribution[dc_to_assign].push_back(ecount);
		if(tag == 1 || tag == 2)
			orig_v[dc_to_assign] ++;
		if(tag == 3)
			orig_v[dc_to_assign] += 2;
		omp_unset_lock(&writelock1[dc_to_assign]);
		//printf("thread %d: edge %d of %d, place in dc %d\n",threadid,ecount,numedges,edge_distr[ecount]);
		
		if(ecount % 1000000 == 0){
			printf("thread %d addressed %d edges\n", threadid, ecount);
		}
		omp_unset_lock(&writelock[lockid2]);	
		omp_unset_lock(&writelock[lockid1]);	
	}	

	//update contained vertices for each dc
	std::vector<std::vector<int> > local_degree_counter;
	for(int i=0; i<vertexsize; i++){
		local_degree_counter.push_back(std::vector<int>(DCs.size(),0));
	}	
	#pragma omp parallel for
	for(int dc = 0; dc<DCs.size(); dc++){
		for(int ei=0; ei<edge_distribution[dc].size(); ei++){
			int e = edge_distribution[dc][ei];
			int src = source(dag->random_access_edges[e],*dag->g);
			int tgt = target(dag->random_access_edges[e],*dag->g);
			omp_set_lock(&writelock1[dc]);
			local_degree_counter[src][dc]++;
			local_degree_counter[tgt][dc]++;
			omp_unset_lock(&writelock1[dc]);
		}
	}
	for(int i=0; i<numoflocks; i++)
		omp_destroy_lock(&writelock[i]);
	for(int i=0; i<numoflocks1; i++)
		omp_destroy_lock(&writelock1[i]);	

	/* master location is selected as the replica with largest local degree */
#pragma omp parallel for
	for(int vi=0; vi<vertexsize; vi++){
		int master = -1;//(*dag->g)[vi].location_id;
		int largest_degree = 0;
		for(int dc=0; dc<DCs.size(); dc++){
			if(local_degree_counter[vi][dc] > largest_degree){
				largest_degree = local_degree_counter[vi][dc];
				master = dc;
			}
		}
		(*dag->g)[vi].master_location = master;
		// no edge to the v
		if(master == -1){
			(*dag->g)[vi].master_location = (*dag->g)[vi].location_id;
		}
	}
}


void Optimizer::InputDataMovement(){
	
	/** Step 1: vertex clustering **/
	
	/** Step 2: calculate the value for each cluster */
		
}

void Optimizer::PartitionPlacement(){
	
	/**
	* Step 1: identify the bottlenecked dcs in gather and apply stages
	*/
	int gather_btlneck = -1;
	int apply_btlneck = -1;
	int gather_link = -1;
	int apply_link = -1;
	int vertexsize = dag->num_vertices;
	in_nbr_distribution = std::vector<std::vector<int> >(vertexsize,std::vector<int>(DCs.size(),0));
	out_nbr_distribution = std::vector<std::vector<int> >(vertexsize,std::vector<int>(DCs.size(),0));

	omp_set_num_threads(N_THREADS);
	for(int dc=0; dc<DCs.size(); dc++){
		for(int ei=0; ei<edge_distribution[dc].size(); ei++){
			Edge e = dag->random_access_edges[edge_distribution[dc][ei]];
			int src = source(e, *dag->g);
			int tgt = target(e, *dag->g);;
			in_nbr_distribution[tgt][dc]++;
			out_nbr_distribution[src][dc]++;
		}
	}

	std::vector<int> tmp_master = std::vector<int>(vertexsize,-1);
#pragma omp parallel for
	for(int i=0; i<vertexsize; i++){
		tmp_master[i] = (*dag->g)[i].master_location;
	}
	std::vector<std::vector<int> > contained_vertices = std::vector<std::vector<int> >(DCs.size());
#pragma omp parallel for
	for(int dc = 0; dc < DCs.size(); dc++){
		std::vector<int> p_contained_vertices;
		for(int ei=0; ei<edge_distribution[dc].size(); ei++){
			int ecount = edge_distribution[dc][ei];
			int src = source(dag->random_access_edges[ecount], *dag->g);
			int tgt = target(dag->random_access_edges[ecount], *dag->g);	
			p_contained_vertices.push_back(src);
			p_contained_vertices.push_back(tgt);			
		}
		std::sort(p_contained_vertices.begin(),p_contained_vertices.end());
		std::vector<int>::iterator it = std::unique(p_contained_vertices.begin(),p_contained_vertices.end());
		p_contained_vertices.resize(std::distance(p_contained_vertices.begin(),it));
		contained_vertices[dc]=p_contained_vertices;
	}
	bool cont = false;
	do{
		float cur_wan_usage = 0.0;
		float gather_max_time = 0.0;
		float apply_max_time = 0.0;
		
		#pragma omp parallel for
		for(int dc=0; dc<DCs.size(); dc++){
			DCs[dc]->g_dnload_data = 0.0;
			DCs[dc]->g_upload_data = 0.0;
			DCs[dc]->a_dnload_data = 0.0;
			DCs[dc]->a_upload_data = 0.0;	
			DCs[dc]->data_size = 0.0;			
			for(int vi=0; vi<contained_vertices[dc].size(); vi++){
				//gather: master download, non-master upload
				//apply: master upload, non-master download
				int gv = contained_vertices[dc][vi];
				if((*dag->g)[gv].master_location == dc){
					if(!myapp->BiDirected){
						//gather neighbor: in_edges
						for(int di=0; di<DCs.size(); di++){
							if(dc != di){
								DCs[dc]->g_dnload_data += in_nbr_distribution[gv][di] * myapp->msg_size;// * myapp->red_rate;
							}
						}					
					}else {
						//gather neighbor: all neighbors
						for(int di=0; di<DCs.size(); di++){
							if(dc != di){
								DCs[dc]->g_dnload_data += (in_nbr_distribution[gv][di] + out_nbr_distribution[gv][di]) 
									* myapp->msg_size;// * myapp->red_rate;							
							}
						}
					}
					DCs[dc]->a_upload_data += myapp->msg_size * ((*dag->g)[gv].replica_location.size() - 1);
				}else{
					if(!myapp->BiDirected){
						//gather neighbor: in_edges					
						DCs[dc]->g_upload_data += in_nbr_distribution[gv][dc] * myapp->msg_size;// * myapp->red_rate;					
					}else{
						//gather neighbor: all neighbors
						DCs[dc]->g_upload_data += (in_nbr_distribution[gv][dc] + out_nbr_distribution[gv][dc]) * myapp->msg_size;// * myapp->red_rate;							
					}
					DCs[dc]->a_dnload_data += myapp->msg_size;
				}
				DCs[dc]->data_size += (*dag->g)[gv].data->input_size;
			}//end for each vi			
		}//end for each dc
		for(int dc=0; dc<DCs.size(); dc++){
			cur_wan_usage += DCs[dc]->g_upload_data + DCs[dc]->a_upload_data;
			float g_upload_time = DCs[dc]->g_upload_data / DCs[dc]->upload_band;
			float g_dnload_time = DCs[dc]->g_dnload_data / DCs[dc]->download_band;
			float a_upload_time = DCs[dc]->a_upload_data / DCs[dc]->upload_band;
			float a_dnload_time = DCs[dc]->a_dnload_data / DCs[dc]->download_band;
			if(g_upload_time > gather_max_time){
				gather_max_time = g_upload_time;
				gather_btlneck = dc;
				gather_link = 0;
			}
			if(g_dnload_time > gather_max_time){
				gather_max_time = g_dnload_time;
				gather_btlneck = dc;
				gather_link = 1;
			}		
			if(a_upload_time > apply_max_time){
				apply_max_time = a_upload_time;
				apply_btlneck = dc;
				apply_link = 0;
			}
			if(a_dnload_time > apply_max_time){
				apply_max_time = a_dnload_time;
				apply_btlneck = dc;
				apply_link = 1;
			}
		}
	
		//debug
		//std::cout<<"step 1 of parition placement is done." <<std::endl;
		/**
		* Step 2: move graph partitions out of bottlenecked dcs 
		*/
		std::cout << "bottlenecked datacenters for gather and apply: " << gather_btlneck << ",  "  << apply_btlneck << std::endl;	  
			  
		std::vector<int> alldcs = std::vector<int>(DCs.size());
		for(int dc=0; dc<DCs.size(); dc++){
			alldcs[dc] = dc;
		}
		int dc_to_switch = -1;
		int btlneck = -1;
		if(gather_btlneck == apply_btlneck){		
			std::pair<float,int> gain = EstimateGain(gather_btlneck,alldcs);
			std::cout << "switch to datacenter: " << gain.second << ", " << gain.first << std::endl;
			if(gain.second != -1){
				float new_wan_usage = cur_wan_usage + DCs[gain.second]->data_size + DCs[gather_btlneck]->data_size;
				if(new_wan_usage <= myapp->budget){
					dc_to_switch = gain.second;
					btlneck = apply_btlneck;
					//switch gather_btlneck and gain.second
					std::cout << "switching datacenter " << gather_btlneck << " with " << gain.second << std::endl;
				}
			}
		}else{
			std::pair<float,int> gain1 = EstimateGain(gather_btlneck,alldcs);
			std::pair<float,int> gain2 = EstimateGain(apply_btlneck,alldcs);
			std::cout << "switch to datacenters: " << gain1.second << ", " << gain1.first
				<< ";	" << gain2.second << ", " << gain2.first << std::endl;
			int larger_dc = gain1.first > gain2.first ? gain1.second : gain2.second;
			if(larger_dc != -1){
				btlneck = gain1.first > gain2.first ? gather_btlneck : apply_btlneck;
				float new_wan_usage = cur_wan_usage + DCs[larger_dc]->data_size + DCs[btlneck]->data_size;
				if(new_wan_usage <= myapp->budget){
					//switch btlneck and dc_to_switch
					dc_to_switch = larger_dc;
					std::cout << "switching datacenter " << btlneck << " with " << dc_to_switch << std::endl;
				}
			}
		}
		if(dc_to_switch != -1 && btlneck != -1){
			//update vertex replica location and in/out nbrs
		#pragma omp parallel for
			for(int vi=0; vi<contained_vertices[btlneck].size(); vi++){
				int gv = contained_vertices[btlneck][vi];
				in_nbr_distribution[gv][btlneck] = 0;
				out_nbr_distribution[gv][btlneck] = 0;
				for(int ri =0; ri < (*dag->g)[gv].replica_location.size(); ri++){
					if((*dag->g)[gv].replica_location[ri] == btlneck){
						(*dag->g)[gv].replica_location[ri] = dc_to_switch;
						break;
					}
				}
				if((*dag->g)[gv].master_location == btlneck)
					tmp_master[gv] = dc_to_switch;
			}
		#pragma omp parallel for
			for(int vi=0; vi<contained_vertices[dc_to_switch].size(); vi++){				
				int gv = contained_vertices[dc_to_switch][vi];
				in_nbr_distribution[gv][dc_to_switch] = 0;
				out_nbr_distribution[gv][dc_to_switch] = 0;
				for(int ri =0; ri < (*dag->g)[gv].replica_location.size(); ri++){
					if((*dag->g)[gv].replica_location[ri] == dc_to_switch){
						(*dag->g)[gv].replica_location[ri] = btlneck;
						break;
					}
				}
				if((*dag->g)[gv].master_location == dc_to_switch)
					tmp_master[gv] = btlneck;
			}
		#pragma omp parallel for
			for(int vi=0; vi<contained_vertices[btlneck].size(); vi++){
				int gv = contained_vertices[btlneck][vi];
				(*dag->g)[gv].master_location = tmp_master[gv];
			}
		#pragma omp parallel for
			for(int vi=0; vi<contained_vertices[dc_to_switch].size(); vi++){
				int gv = contained_vertices[dc_to_switch][vi];
				(*dag->g)[gv].master_location = tmp_master[gv];
			}		
			std::vector<int> tmp = edge_distribution[btlneck];
			edge_distribution[btlneck] = edge_distribution[dc_to_switch];
			edge_distribution[dc_to_switch] = tmp;
			for(int ei=0; ei<edge_distribution[btlneck].size(); ei++){
				int ecount = edge_distribution[btlneck][ei];
				int src = source(dag->random_access_edges[ecount], *dag->g);
				int tgt = target(dag->random_access_edges[ecount], *dag->g);
				in_nbr_distribution[tgt][btlneck]++;
				out_nbr_distribution[src][btlneck]++;
			}
			for(int ei=0; ei<edge_distribution[dc_to_switch].size(); ei++){
				int ecount = edge_distribution[dc_to_switch][ei];
				int src = source(dag->random_access_edges[ecount], *dag->g);
				int tgt = target(dag->random_access_edges[ecount], *dag->g);
				in_nbr_distribution[tgt][dc_to_switch]++;
				out_nbr_distribution[src][dc_to_switch]++;
			}
			float tmp_data_size = DCs[btlneck]->data_size;
			DCs[btlneck]->data_size = DCs[dc_to_switch]->data_size;
			DCs[dc_to_switch]->data_size = tmp_data_size;			
			tmp = contained_vertices[btlneck];
			contained_vertices[btlneck] = contained_vertices[dc_to_switch];
			contained_vertices[dc_to_switch] = tmp;
			// ContextSwitch(btlneck,dc_to_switch,DCs);
			cont = true;
		}else {
			cont = false;				
		}		
	}while(cont);

	for(int v=0; v<dag->num_vertices; v++){
		int master = (*dag->g)[v].master_location;
		if(master == -1){
			(*dag->g)[v].master_location = (*dag->g)[v].location_id;
		}
		if(std::find((*dag->g)[v].replica_location.begin(),(*dag->g)[v].replica_location.end(),master) == (*dag->g)[v].replica_location.end()){
			printf("vertex %d master location error in placement\n",v);
		}
		in_nbr_distribution[v].clear();
		out_nbr_distribution[v].clear();
	}
	in_nbr_distribution.clear();
	out_nbr_distribution.clear();
}

/** 
* Step 3: migrate edges from bottlenecked datacenters to other locations
* partition replacement cannot further improve performance
* given the bottlenecked dcs and links, move edges out of the bottlenecks
*/
void Optimizer::EdgeMigration(){
	
	omp_set_num_threads(N_THREADS);
		
	bool cont = true;
	int  total_migrated = 0;
	std::vector<int> num_migrated = std::vector<int>(N_THREADS,0);
	while(cont){
		float gather_max_time = 0.0;
		float apply_max_time = 0.0;
		float gather_sec_max = 0.0;
		float apply_sec_max = 0.0;
		int gather_btlneck = -1;
		int apply_btlneck = -1;
		int gather_link = -1;
		int apply_link = -1;
		
		/**
		* Step 1: Identify bottlenecks
		*/
		std::vector<std::vector<int> > contained_vertices = std::vector<std::vector<int> >(DCs.size());
		in_nbr_distribution = std::vector<std::vector<int> >(dag->num_vertices,std::vector<int>(DCs.size(),0));
		out_nbr_distribution = std::vector<std::vector<int> >(dag->num_vertices,std::vector<int>(DCs.size(),0));
	#pragma omp parallel for
		for(int dc = 0; dc < DCs.size(); dc++){
			std::vector<int> p_contained_vertices;
			for(int ei=0; ei<edge_distribution[dc].size(); ei++){
				int ecount = edge_distribution[dc][ei];
				int src = source(dag->random_access_edges[ecount], *dag->g);
				int tgt = target(dag->random_access_edges[ecount], *dag->g);	
				p_contained_vertices.push_back(src);
				p_contained_vertices.push_back(tgt);
				in_nbr_distribution[tgt][dc]++;
				out_nbr_distribution[src][dc]++;				
			}
			std::sort(p_contained_vertices.begin(),p_contained_vertices.end());
			std::vector<int>::iterator it = std::unique(p_contained_vertices.begin(),p_contained_vertices.end());
			p_contained_vertices.resize(std::distance(p_contained_vertices.begin(),it));
			contained_vertices[dc]=p_contained_vertices;
		}
		
		#pragma omp parallel for
		for(int dc=0; dc<DCs.size(); dc++){
			DCs[dc]->g_dnload_data = DCs[dc]->g_upload_data = DCs[dc]->a_dnload_data= DCs[dc]->a_upload_data =0.0;
			for(int vi=0; vi<contained_vertices[dc].size(); vi++){
				//gather: master download, non-master upload
				//apply: master upload, non-master download
				int gv = contained_vertices[dc][vi];
				if((*dag->g)[gv].master_location == dc){
					if(!myapp->BiDirected){
						//gather neighbor: in_edges
						for(int di=0; di<DCs.size(); di++){
							if(dc != di){
								DCs[dc]->g_dnload_data += in_nbr_distribution[gv][di] * myapp->msg_size;// * myapp->red_rate;
							}
						}					
					}else {
						//gather neighbor: all neighbors
						for(int di=0; di<DCs.size(); di++){
							if(dc != di){
								DCs[dc]->g_dnload_data += (in_nbr_distribution[gv][di] + out_nbr_distribution[gv][di]) 
									* myapp->msg_size;// * myapp->red_rate;							
							}
						}
					}
					DCs[dc]->a_upload_data += myapp->msg_size * ((*dag->g)[gv].replica_location.size() - 1);
				}else{
					if(!myapp->BiDirected){
						//gather neighbor: in_edges					
						DCs[dc]->g_upload_data += in_nbr_distribution[gv][dc] * myapp->msg_size;// * myapp->red_rate;					
					}else {
						//gather neighbor: all neighbors
						DCs[dc]->g_upload_data += (in_nbr_distribution[gv][dc] + out_nbr_distribution[gv][dc]) * myapp->msg_size;// * myapp->red_rate;							
					}
					DCs[dc]->a_dnload_data += myapp->msg_size;
				}
				DCs[dc]->data_size += (*dag->g)[gv].data->input_size;
			}//end for each vi			
		}//end for each dc
		for(int dc=0; dc<DCs.size(); dc++){
			float g_upload_time = DCs[dc]->g_upload_data / DCs[dc]->upload_band;
			float g_dnload_time = DCs[dc]->g_dnload_data / DCs[dc]->download_band;
			float a_upload_time = DCs[dc]->a_upload_data / DCs[dc]->upload_band;
			float a_dnload_time = DCs[dc]->a_dnload_data / DCs[dc]->download_band;
			if(g_upload_time > gather_max_time){
				gather_max_time = g_upload_time;
				gather_btlneck = dc;
				gather_link = 0;
			}
			if(g_dnload_time > gather_max_time){
				gather_max_time = g_dnload_time;
				gather_btlneck = dc;
				gather_link = 1;
			}		
			if(g_upload_time < gather_max_time && g_upload_time > gather_sec_max)
				gather_sec_max = g_upload_time;
			if(g_dnload_time < gather_max_time && g_dnload_time > gather_sec_max)
				gather_sec_max = g_dnload_time;
			if(a_upload_time > apply_max_time){
				apply_max_time = a_upload_time;
				apply_btlneck = dc;
				apply_link = 0;
			}
			if(a_dnload_time > apply_max_time){
				apply_max_time = a_dnload_time;
				apply_btlneck = dc;
				apply_link = 1;
			}
			if(a_upload_time < apply_max_time && a_upload_time > apply_sec_max)
				apply_sec_max = a_upload_time;
			if(a_dnload_time < apply_max_time && a_dnload_time > apply_sec_max)
				apply_sec_max = a_dnload_time;
		}	
		printf("****************************edge migration step****************************\n");
		printf("gather btlneck: dc %d link %d\n", gather_btlneck, gather_link);
		printf("apply btlneck: dc %d link %d\n", apply_btlneck, apply_link);
		printf("gather max and sec max: %f, %f\n",gather_max_time,gather_sec_max);
		printf("apply max and sec max: %f, %f\n",apply_max_time,apply_sec_max);
		
		if(gather_btlneck == -1 && apply_btlneck == -1)
			break;
	
		std::vector<std::pair<int,float> > masters;
		std::vector<std::pair<int,float> > mirrors;
		std::vector<std::pair<int,float> > priorities;
		int btldc = -1;
		int btllink = -1;
		//bound on the same dc and same link
		if(gather_btlneck == apply_btlneck && gather_link == apply_link){			
			btldc = gather_btlneck;
			btllink = gather_link;	
			
			for(int i=0; i<contained_vertices[btldc].size(); i++){
				int vid = contained_vertices[btldc][i];
				if((*dag->g)[vid].master_location == btldc){
					masters.push_back(std::pair<int,float> (vid,0.0));
				}else{
					mirrors.push_back(std::pair<int,float> (vid,0.0));
				}
			}			
							
			if(btllink == 0){
				// uplink
				// estimate the gain of the gather stage, mirror replicas
				for(int i=0; i< mirrors.size(); i++){
					if(!myapp->BiDirected){
						//gather neighbor: in_edges					
						mirrors[i].second = in_nbr_distribution[mirrors[i].first][btldc] * myapp->msg_size;// * myapp->red_rate;					
					}else{
						//gather neighbor: all neighbors
						mirrors[i].second = (in_nbr_distribution[mirrors[i].first][btldc] + out_nbr_distribution[mirrors[i].first][btldc]) * myapp->msg_size;// * myapp->red_rate;							
					}
					priorities.push_back(mirrors[i]);
				}				
				// estimate the gain of apply stage, master replicas
				for(int i=0; i<masters.size();i++){
					masters[i].second = myapp->msg_size * ((*dag->g)[masters[i].first].replica_location.size() - 1);
					priorities.push_back(masters[i]);
				}										
			}else if(btllink == 1){
				//downlink
				// estimate the gain of the gather stage, master replicas
				#pragma omp parallel for
				for(int i=0; i< masters.size(); i++){
					if(!myapp->BiDirected){
						//gather neighbor: in_edges
						for(int di=0; di<DCs.size(); di++){
							if(btldc != di){
								masters[i].second += in_nbr_distribution[masters[i].first][di] * myapp->msg_size;// * myapp->red_rate;
							}
						}					
					}else{
						//gather neighbor: all neighbors
						for(int di=0; di<DCs.size(); di++){
							if(btldc != di){
								masters[i].second += (in_nbr_distribution[masters[i].first][di] + out_nbr_distribution[masters[i].first][di]) 
									* myapp->msg_size;// * myapp->red_rate;							
							}
						}
					}	
					#pragma omp critical					
					priorities.push_back(masters[i]);
				}				
				// estimate the gain of apply stage, mirror replicas
				for(int i=0; i<mirrors.size();i++){
					mirrors[i].second = myapp->msg_size;
					priorities.push_back(mirrors[i]);
				}						
			}else{
				printf("no bottleneck link?\n");
				exit(1);
			}		
		}else{
			//bound on different links, select one link to resolve
			float gather_improve = gather_max_time - gather_sec_max;
			float apply_improve = apply_max_time - apply_sec_max;
			if(gather_improve > apply_improve){
				btldc = gather_btlneck;
				btllink = gather_link;
				for(int i=0; i<contained_vertices[btldc].size(); i++){
					int vid = contained_vertices[btldc][i];
					if((*dag->g)[vid].master_location == btldc){
						masters.push_back(std::pair<int,float> (vid,0.0));
					}else{
						mirrors.push_back(std::pair<int,float> (vid,0.0));
					}
				}
				if(gather_link == 0){
					//uplink
					#pragma omp parallel for
					for(int i=0; i< mirrors.size(); i++){
						if(!myapp->BiDirected){
							//gather neighbor: in_edges					
							mirrors[i].second = in_nbr_distribution[mirrors[i].first][btldc] * myapp->msg_size;// * myapp->red_rate;					
						}else {
							//gather neighbor: all neighbors
							mirrors[i].second = (in_nbr_distribution[mirrors[i].first][btldc] + out_nbr_distribution[mirrors[i].first][btldc]) * myapp->msg_size;// * myapp->red_rate;							
						}							
					}	
				}else if(gather_link == 1){
					//downlink
					#pragma omp parallel for
					for(int i=0; i< masters.size(); i++){
						if(!myapp->BiDirected){
							//gather neighbor: in_edges
							for(int di=0; di<DCs.size(); di++){
								if(btldc != di){
									masters[i].second += in_nbr_distribution[masters[i].first][di] * myapp->msg_size;// * myapp->red_rate;
								}
							}					
						}else{
							//gather neighbor: all neighbors
							for(int di=0; di<DCs.size(); di++){
								if(btldc != di){
									masters[i].second += (in_nbr_distribution[masters[i].first][di] + out_nbr_distribution[masters[i].first][di]) 
										* myapp->msg_size;// * myapp->red_rate;							
								}
							}
						}
					}	
				}else{
					printf("no bottleneck link?\n");
					exit(1);
				}	
			}
			else{
				btldc = apply_btlneck;
				btllink = apply_link;
				for(int i=0; i<contained_vertices[btldc].size(); i++){
					int vid = contained_vertices[btldc][i];
					if((*dag->g)[vid].master_location == btldc){
						masters.push_back(std::pair<int,float> (vid,0.0));
					}else{
						mirrors.push_back(std::pair<int,float> (vid,0.0));
					}
				}
				if(apply_link == 0){
					//uplink
					// estimate the gain of apply stage, master replicas
					for(int i=0; i<masters.size();i++){
						masters[i].second += myapp->msg_size * ((*dag->g)[masters[i].first].replica_location.size() - 1);
						priorities.push_back(masters[i]);
					}
				}else if(apply_link == 1){
					//downlink
					// estimate the gain of apply stage, mirror replicas
					for(int i=0; i<mirrors.size();i++){
						mirrors[i].second += myapp->msg_size;
						priorities.push_back(mirrors[i]);
					}
				}else{
					printf("no bottleneck link?\n");
					exit(1);
				}						
			}							
		}//endif bound on diff links
		//iteratively remove the head of priorities
		//until btldc is no longer the bottleneck
		std::sort(priorities.begin(),priorities.end(),mycompare);
		int candi_size = priorities.size();
		if(candi_size > 100){
			//select the top 10-percent
			candi_size = std::ceil(candi_size*0.1);
			candi_size = candi_size>100?100:candi_size;
			priorities.resize(candi_size);
		}
		// printf("prioritie queue to remove\n");
		// for(int pi=0; pi<priorities.size(); pi++){
			// printf("%d, %f\n",priorities[pi].first, priorities[pi].second);
		// }		
		int suc = 0;
		
		for(int pvi = 0; pvi<priorities.size(); pvi++){
			bool mig_suc = false;
			int v_to_rm = priorities[pvi].first;
			int orig_in_nbr = in_nbr_distribution[v_to_rm][btldc];
			int orig_out_nbr = out_nbr_distribution[v_to_rm][btldc];
			
			float cur_time = EstimateTransferTime();
			//printf("****************time before removing vertex %d from dc %d is %f*****************\n",v_to_rm,btldc,cur_time);
			
			if((*dag->g)[v_to_rm].master_location == btldc){
				if(!myapp->BiDirected){
					//gather neighbor: in_edges
					for(int di=0; di<DCs.size(); di++){
						if(btldc != di){
							DCs[btldc]->g_dnload_data -= in_nbr_distribution[v_to_rm][di] * myapp->msg_size;// * myapp->red_rate;
						}
					}					
				}else {
					//gather neighbor: all neighbors
					for(int di=0; di<DCs.size(); di++){
						if(btldc != di){
							DCs[btldc]->g_dnload_data -= (in_nbr_distribution[v_to_rm][di] + out_nbr_distribution[v_to_rm][di]) 
								* myapp->msg_size;// * myapp->red_rate;							
						}
					}
				}
				DCs[btldc]->a_upload_data -= myapp->msg_size * ((*dag->g)[v_to_rm].replica_location.size() - 1);
			}else{
				float delta_traffic = 0;
				if(!myapp->BiDirected){
					//gather neighbor: in_edges					
					delta_traffic = in_nbr_distribution[v_to_rm][btldc] * myapp->msg_size;// * myapp->red_rate;	
				}else{
					//gather neighbor: all neighbors
					delta_traffic = (in_nbr_distribution[v_to_rm][btldc] + out_nbr_distribution[v_to_rm][btldc]) * myapp->msg_size;// * myapp->red_rate;
				}
				DCs[btldc]->g_upload_data -= delta_traffic;
				DCs[(*dag->g)[v_to_rm].master_location]->g_dnload_data -= delta_traffic;					
				DCs[btldc]->a_dnload_data -= myapp->msg_size;
				DCs[(*dag->g)[v_to_rm].master_location]->a_upload_data -= myapp->msg_size;
			}		
			for(int i=0; i<(*dag->g)[v_to_rm].replica_location.size(); i++){
				if((*dag->g)[v_to_rm].replica_location[i] == btldc){
					(*dag->g)[v_to_rm].replica_location.erase((*dag->g)[v_to_rm].replica_location.begin()+i);
					break;
				}
			}
			//update master if neccesary
			std::vector<int> orig_masters;
			orig_masters.push_back((*dag->g)[v_to_rm].master_location);
			in_nbr_distribution[v_to_rm][btldc] = 0;
			out_nbr_distribution[v_to_rm][btldc] = 0;
			
			std::vector<int > connected_edges;			
			for(int i=0; i<edge_distribution[btldc].size(); i++){
				int ecount = edge_distribution[btldc][i];
				int src = source(dag->random_access_edges[ecount],*dag->g);
				int tgt = target(dag->random_access_edges[ecount],*dag->g);
				if(src == v_to_rm || tgt == v_to_rm){
					if(src == v_to_rm ){
						in_nbr_distribution[tgt][btldc] --;
						orig_masters.push_back((*dag->g)[tgt].master_location);
					}else{
						out_nbr_distribution[src][btldc] --;
						orig_masters.push_back((*dag->g)[src].master_location);
					} 					
					connected_edges.push_back(edge_distribution[btldc][i]);
					edge_distribution[btldc].erase(edge_distribution[btldc].begin()+i);
					i--;					
				}
			}
			
			for(int i=0; i<orig_masters.size(); i++){
				int gv = -1;
				if(i==0){
					gv = v_to_rm;
				}else{
					int lsrc = source(dag->random_access_edges[connected_edges[i-1]],*dag->g);
					int ltgt = target(dag->random_access_edges[connected_edges[i-1]],*dag->g);		
					gv = lsrc == v_to_rm ? ltgt : lsrc;
				}
				if(orig_masters[i] == btldc){
					int largest_degree = -1;
					for(int dc=0; dc<(*dag->g)[gv].replica_location.size(); dc++){
						//if(dc != btldc){
							int newmaster = (*dag->g)[gv].replica_location[dc];
							int local_degree = in_nbr_distribution[gv][newmaster] + out_nbr_distribution[gv][newmaster];
							if(local_degree > largest_degree){
								largest_degree = local_degree;
								(*dag->g)[gv].master_location = newmaster;
							}
						//}
					}
					//printf("master location of v %d changed from %d to %d\n", gv,orig_masters[i],(*dag->g)[gv].master_location);	
				}				
			}
			//end of change status
			/********************Migrate edges: group k edges and migrate them to the same dc at the same time***********************************/
			std::pair<int,std::vector<int> > results;	
			std::vector<float> esttimes;
			std::vector<int> notmoveto;
			notmoveto.push_back(btldc);	
			//migrate all edges at the same time
			results = EdgeMigrate(connected_edges,notmoveto,cur_time);	
			int dc_to_place = results.first;
			if(dc_to_place != -1){		
				mig_suc = true;	
				total_migrated += connected_edges.size();				
				// change the status
				for(int ei=0; ei<connected_edges.size(); ei++){
					int src = source(dag->random_access_edges[connected_edges[ei]],*dag->g);
					int tgt = target(dag->random_access_edges[connected_edges[ei]],*dag->g);
					int tag = results.second[ei];
					if(tag == 0){
						if((*dag->g)[src].master_location != dc_to_place && myapp->BiDirected){
							DCs[dc_to_place]->g_upload_data += myapp->msg_size;// * myapp->red_rate;
						}						
						if((*dag->g)[tgt].master_location != dc_to_place){
							DCs[dc_to_place]->g_upload_data += myapp->msg_size;// * myapp->red_rate;
						}	
					}else if(tag == 1){
						(*dag->g)[src].replica_location.push_back(dc_to_place);
						
						if((*dag->g)[tgt].master_location != dc_to_place)
							DCs[dc_to_place]->g_upload_data += myapp->msg_size;// * myapp->red_rate;
						//add a src mirror: increase g_upload_data and a_dnload_data
						if(myapp->BiDirected){
							DCs[dc_to_place]->g_upload_data += myapp->msg_size;									
							DCs[(*dag->g)[src].master_location]->g_dnload_data += myapp->msg_size;
						}
						DCs[dc_to_place]->a_dnload_data += myapp->msg_size;	
						DCs[(*dag->g)[src].master_location]->a_upload_data += myapp->msg_size;	
					}else if(tag == 2){
						(*dag->g)[tgt].replica_location.push_back(dc_to_place);
						
						if((*dag->g)[src].master_location != dc_to_place && myapp->BiDirected)
							DCs[dc_to_place]->g_upload_data += myapp->msg_size;// * myapp->red_rate;
						//add a tgt mirror: increase g_upload_data and a_dnload_data											
						DCs[dc_to_place]->g_upload_data += myapp->msg_size;	
						DCs[dc_to_place]->a_dnload_data += myapp->msg_size;	
						DCs[(*dag->g)[tgt].master_location]->g_dnload_data += myapp->msg_size;
						DCs[(*dag->g)[tgt].master_location]->a_upload_data += myapp->msg_size;		
					}else if(tag == 3){
						(*dag->g)[src].replica_location.push_back(dc_to_place);
						(*dag->g)[tgt].replica_location.push_back(dc_to_place);
						
						if(myapp->BiDirected){
							DCs[dc_to_place]->g_upload_data += myapp->msg_size;// * myapp->red_rate;					
							DCs[(*dag->g)[src].master_location]->g_dnload_data += myapp->msg_size;
						}
						DCs[dc_to_place]->a_dnload_data += myapp->msg_size;
						DCs[(*dag->g)[src].master_location]->a_upload_data += myapp->msg_size;
						DCs[dc_to_place]->g_upload_data += myapp->msg_size;// * myapp->red_rate;
						DCs[dc_to_place]->a_dnload_data += myapp->msg_size;			
						DCs[(*dag->g)[tgt].master_location]->g_dnload_data += myapp->msg_size;
						DCs[(*dag->g)[tgt].master_location]->a_upload_data += myapp->msg_size;
					}	
					in_nbr_distribution[tgt][dc_to_place] ++;
					out_nbr_distribution[src][dc_to_place] ++;
					edge_distribution[dc_to_place].push_back(connected_edges[ei]);
					// if(in_nbr_distribution[src][(*dag->g)[src].master_location] == 0 && out_nbr_distribution[src][(*dag->g)[src].master_location] == 0)
						// (*dag->g)[src].master_location = dc_to_place;
					// if(in_nbr_distribution[tgt][(*dag->g)[tgt].master_location] == 0 && out_nbr_distribution[tgt][(*dag->g)[tgt].master_location] == 0)
						// (*dag->g)[tgt].master_location = dc_to_place;					
				}//end for each result	
				for(int ei=0; ei<connected_edges.size(); ei++){
					int src = source(dag->random_access_edges[connected_edges[ei]],*dag->g);
					int tgt = target(dag->random_access_edges[connected_edges[ei]],*dag->g);
					for(int it=0; it<2; it++){
						int node = -1;
						if(it == 0) node = src;
						else node = tgt;
						int master = (*dag->g)[node].master_location;
						if(std::find((*dag->g)[node].replica_location.begin(),(*dag->g)[node].replica_location.end(),master)==(*dag->g)[node].replica_location.end()){
							int largest_degree = -1;
							for(int di=0; di<(*dag->g)[node].replica_location.size(); di++){							
								int local_degree = in_nbr_distribution[node][(*dag->g)[node].replica_location[di]] + out_nbr_distribution[node][(*dag->g)[node].replica_location[di]];
								if(local_degree > largest_degree){
									largest_degree = local_degree;
									(*dag->g)[node].master_location = (*dag->g)[node].replica_location[di];
								}
							}
						}
					}
				}				
			}else{
				//don't remove, recover context
				(*dag->g)[v_to_rm].replica_location.push_back(btldc);
				for(int i=0; i<orig_masters.size(); i++){
					int gv = -1;
					if(i==0){
						gv = v_to_rm;
					}else{
						int lsrc = source(dag->random_access_edges[connected_edges[i-1]],*dag->g);
						int ltgt = target(dag->random_access_edges[connected_edges[i-1]],*dag->g);		
						gv = lsrc == v_to_rm ? ltgt : lsrc;
					}
					if(orig_masters[i] == btldc){						
						(*dag->g)[gv].master_location = btldc;
					}	
				}
				in_nbr_distribution[v_to_rm][btldc] = orig_in_nbr;
				out_nbr_distribution[v_to_rm][btldc] = orig_out_nbr;
				
				for(int i=0; i<connected_edges.size(); i++){
					edge_distribution[btldc].push_back(connected_edges[i]);
					int lsrc = source(dag->random_access_edges[connected_edges[i]],*dag->g);
					int ltgt = target(dag->random_access_edges[connected_edges[i]],*dag->g);		
					if(lsrc == v_to_rm ){
						in_nbr_distribution[ltgt][btldc] ++;
					}else{
						out_nbr_distribution[lsrc][btldc] ++;
					} 
				}				
				
				if(orig_masters[0] == btldc){
					if(!myapp->BiDirected){
						//gather neighbor: in_edges
						for(int di=0; di<DCs.size(); di++){
							if(btldc != di){
								DCs[btldc]->g_dnload_data += in_nbr_distribution[v_to_rm][di] * myapp->msg_size;// * myapp->red_rate;
							}
						}					
					}else {
						//gather neighbor: all neighbors
						for(int di=0; di<DCs.size(); di++){
							if(btldc != di){
								DCs[btldc]->g_dnload_data += (in_nbr_distribution[v_to_rm][di] + out_nbr_distribution[v_to_rm][di]) 
									* myapp->msg_size;// * myapp->red_rate;							
							}
						}
					}
					DCs[btldc]->a_upload_data += myapp->msg_size * ((*dag->g)[v_to_rm].replica_location.size() - 1);
				}else{
					float delta_traffic = 0;
					if(!myapp->BiDirected){
						//gather neighbor: in_edges					
						delta_traffic = in_nbr_distribution[v_to_rm][btldc] * myapp->msg_size;// * myapp->red_rate;	
					}else{
						//gather neighbor: all neighbors
						delta_traffic = (in_nbr_distribution[v_to_rm][btldc] + out_nbr_distribution[v_to_rm][btldc]) * myapp->msg_size;// * myapp->red_rate;
					}
					DCs[btldc]->g_upload_data += delta_traffic;
					DCs[orig_masters[0]]->g_dnload_data += delta_traffic;					
					DCs[btldc]->a_dnload_data += myapp->msg_size;
					DCs[orig_masters[0]]->a_upload_data += myapp->msg_size;
				}
			}
			if(mig_suc) suc ++;
		}//endfor priorities	
		if(suc == 0) cont = false;
	}//end while cont

	// for(int i=0; i<numoflocks; i++)
		// omp_destroy_lock(&writelock[i]);
	// for(int i=0; i<numoflocks1; i++)
		// omp_destroy_lock(&writelock1[i]);	
	
	// for(int i=0; i<N_THREADS; i++){
		// total_migrated += num_migrated[i];
	// }
	
	/*//debug
	std::vector<std::vector<int> > tmp_in_nbr = std::vector<std::vector<int> >(dag->num_vertices,std::vector<int>(DCs.size(),0));
	std::vector<std::vector<int> > tmp_out_nbr = std::vector<std::vector<int> >(dag->num_vertices,std::vector<int>(DCs.size(),0));
	for(int dc=0; dc<DCs.size(); dc++){
		for(int ei=0; ei<edge_distribution[dc].size(); ei++){
			int src = edge_distribution[dc][ei].first;
			int tgt = edge_distribution[dc][ei].second;
			tmp_in_nbr[tgt][dc]++;
			tmp_out_nbr[src][dc]++;
		}
	}
	
	for(int v=0; v<dag->num_vertices; v++){
		for(int dc=0; dc<DCs.size(); dc++){
			if(tmp_in_nbr[v][dc]!=in_nbr_distribution[v][dc] ){
				printf("in nbr of %d in dc %d not correct: %d, %d\n", v, dc, tmp_in_nbr[v][dc],in_nbr_distribution[v][dc]);
			}
			if(tmp_out_nbr[v][dc]!=out_nbr_distribution[v][dc]){
				printf("out nbr of %d in dc %d not correct: %d, %d\n", v, dc, tmp_out_nbr[v][dc],out_nbr_distribution[v][dc]);				
			}
		}
	}	 */ 
	
	for(int v=0; v<dag->num_vertices; v++){
		int master = (*dag->g)[v].master_location;	
		if(master == -1){
			(*dag->g)[v].master_location = (*dag->g)[v].location_id;
		}
		if(std::find((*dag->g)[v].replica_location.begin(),(*dag->g)[v].replica_location.end(),master) == (*dag->g)[v].replica_location.end()){
			printf("vertex %d master location error in migration\n",v);
		}
	}
	printf("In total %d edges are migrated.\n",total_migrated);
}










