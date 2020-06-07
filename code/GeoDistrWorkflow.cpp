// GeoDistrWorkflow.cpp : Defines the entry point for the console application.
//


#include "src/Simulator.h"
#include <ctime>
#include <stdlib.h>  
#include <time.h>  

float OnDemandLag = 60; //seconds
//float budget;
//int Max_Iteration;//for graph algorithms to converge
int Max_Monte_Carlo; //for monte carlo simulation
int N_THREADS;
int numofdcs;
float noise_budget;
bool privacy;

int main(int argc, char* argv[])
{
	/*   //output the Boost version
    std::cout << "Using Boost "
              << BOOST_VERSION / 100000     << "."  // major version
              << BOOST_VERSION / 100 % 1000 << "."  // minor version
              << BOOST_VERSION % 100                // patch level
              << std::endl;
	*/
	std::vector<DataCenter*> DCs;
	if(strcmp(argv[9],"real") == 0){ //for Amazon EC2 experiments
		numofdcs = 8; // real cloud
		for (int i = 0; i < numofdcs; i++){
			DataCenter* DC = new DataCenter(Amazon_EC2_regions(i));
			DC->id = i;
			DC->location = Amazon_EC2_regions(i);
			if(i == 0 || i == 5)
				DC->download_band = 100;	//usually it's 500
			DCs.push_back(DC);		
		}
	} else {
		// LO: set the privacy rank randomly
		numofdcs = 20;//20; // simulation
		 
		// random array 
		//int pri_rank_arrs[] = { 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,20 };
		//int pri_rank_arrs[] = { 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20 };  //combination1, 50% without noise
		//int pri_rank_arrs[] = { 1,1,1,1,1,1,1,1,1,10,11,12,13,14,15,16,17,18,19,20 }; //combination2, 60% without noise
		//int pri_rank_arrs[] = { 1,1,1,1,1,1,1,1,1,1,1,1,1,14,15,16,17,18,19,20 }; //combination3, 70% without noise
		//int pri_rank_arrs[] = { 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,17,18,19,20 }; //combination4, 82% without noise
		//int pri_rank_arrs[] = { 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,20 }; //combination5, 95% without noise
		//int pri_rank_arrs[] = { 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}; //combination6, 100% without noise
		int pri_rank_arrs[] = { 1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5 };
		//int pri_rank_arrs[] = { 1,1,2};
		srand(time(NULL));
		for (int i = 0; i < numofdcs; i++)
		{
			int index = rand() % (numofdcs - i) + i;
			if (index != i)
			{
				int tmp = pri_rank_arrs[i];
				pri_rank_arrs[i] = pri_rank_arrs[index];
				pri_rank_arrs[index] = tmp;
			}
		}
		for(int i=0; i< numofdcs; i++){					
			DataCenter* DC = new DataCenter();
			DC->id = i;
			DC->location = Amazon_EC2_regions(i);
			DC->pri_rank = pri_rank_arrs[i];

			DCs.push_back(DC);

			/*
			//Insert to check the correctness of location0
			std::cout << "DC" << i << "'s location=" << DC->location << std::endl;
			*/
		}


		//Weight is calculated based on the privacy level of the DC
		for (int i = 0; i < numofdcs; i++)
		{
			int weight = 0;
			//DataCenter* DC = new DataCenter();
			//DC->id = i;
			//DC->location = Amazon_EC2_regions(i);
			for (int index = 0; index < numofdcs; index++)
			{
				if (index == i)
					continue;
				if (DCs[i]->pri_rank > DCs[index]->pri_rank)
					weight++;
			}
			DCs[i]->budget_weight = weight;
			weight = 0;
		}

		float sum_weight = 0;
		//Weight_percentage is calculated based on the privacy weight of the DC
		for (int i = 0; i < numofdcs; i++)
			sum_weight += DCs[i]->budget_weight;
		for (int i = 0; i < numofdcs; i++)
		{
			DCs[i]->budget_sum_weight = sum_weight;
			DCs[i]->budget_weight_percentage = DCs[i]->budget_weight / sum_weight;
		}


		if(strcmp(argv[9],"homo")==0){
			//all uploading and downloading bandwidth are the same
			for(int i=0; i<numofdcs; i++){
				DCs[i]->download_band = DCs[i]->upload_band;
			}
		}else if(strcmp(argv[9],"medium")==0){
			//do nothing, all dcs are the same
			//each dc with heterogenous upload/download bandwidth
		}else if(strcmp(argv[9],"high")==0){
			//select 5 of the dcs with lower bandwidth
			//the upload/download in the bottlenecks also heterogenous
			for(int i=10; i<15; i++){
				DCs[i]->download_band *= 0.5;
				DCs[i]->upload_band *= 0.5;
			}
		}	
	}
	
	int degree = atoi(argv[3]); //graph size for synthetic
	distributed_graph* dag = new distributed_graph();
	if (strcmp(argv[1], "livejournal") == 0){
		dag->load_live_journal_graph();
	}
	else if (strcmp(argv[1], "synthetic") == 0){
		dag->load_synthetic_powerlaw(degree, true);//pagerank cares about in edges
	}
	else if	(strcmp(argv[1], "file") == 0){
		string filename;
		if(degree == 1000){
			filename = "powerlaw-1000-2.1";
		}
		else if(degree == 10000){
			filename = "powerlaw-10000-2.1";
		}else if(degree == 50000){
			filename = "powerlaw-50000-2.1";
		}else if(degree == 100000){
			filename = "powerlaw-100000-2.1";
		}else if(degree == 1000000){
			filename = "powerlaw-1000000-2.1";
		}
		dag->load_synthetic_file(filename.c_str());
	}else if(strcmp(argv[1], "uniform") == 0){
		dag->load_from_file("powerlaw-1000-0",1000,0);
	}
	else if (strcmp(argv[1], "facebook") == 0){
		dag->load_facebook_graph();
	}else if(strcmp(argv[1], "googleweb") == 0){
		dag->load_googleweb_graph();
	}else if(strcmp(argv[1], "wiki") == 0){
		dag->load_wiki_graph();
	}else if(strcmp(argv[1],"roadca") == 0){
		dag->load_roadnet_graph();
	}else if(strcmp(argv[1],"twitter") == 0){
		dag->load_twitter_graph();
	}else if(strcmp(argv[1],"p2p") == 0){
		dag->load_p2p_graph();
	}
	else if (strcmp(argv[1], "test") == 0){
		dag->test_graph();
	}
	else if(strcmp(argv[1],"1k") == 0) {
        string filename;
        if(degree == 1000){
            filename = "powerlaw-1000-2.1";
        }
        else if(degree == 10000){
            filename = "powerlaw-10000-2.1";
        }else if(degree == 50000){
            filename = "powerlaw-50000-2.1";
        }else if(degree == 100000){
            std::cout << "Here" << std::endl;
            filename = "powerlaw-100000-2.1";
        }else if(degree == 1000000){
            filename = "powerlaw-1000000-2.1";
        }
	    dag->load_synthetic_file_1k_series(filename.c_str());
	}
	/* log down the original locations of vertices */
	/*std::vector<vector<int> > init_locations(numofdcs);
	for(int i=0; i<dag->num_vertices; i++){
		for(int j=0; j<numofdcs; j++){
			if((*dag->g)[i].location_id == j){
				init_locations[j].push_back(i);
				break;
			}
		} 
	}		
	for(int j=0; j<numofdcs; j++){
		printf("DC %d contains vertices: ",j);
		for(int i=0; i<init_locations[j].size(); i++){
			printf("%d, ", init_locations[j][i]);
		}
		printf("\n");
	}*/

	// dag->mutexlocks = std::vector<LockV>(dag->num_vertices,LockV());
	// dag->mutexlocks = std::vector<std::mutex>(dag->num_vertices);
	int Max_Iteration = atoi(argv[4]); //if==0, converge with threshold
	Max_Monte_Carlo = atoi(argv[5]);
	float budget = atof(argv[6]);
	N_THREADS = atoi(argv[7]);
	noise_budget = atof(argv[8]);
	if(strcmp(argv[10], "privacy") == 0)
		privacy = true;
	else privacy = false;

	BaseApp* myapp;
	if (strcmp(argv[2], "pagerank") == 0){
		myapp = new PageRank();
		myapp->mytype = pagerank;
		for(int v=0; v<dag->num_vertices; v++){
			(*dag->g)[v].data = new PageRankVertexData();
			(*dag->g)[v].accum = new PageRankVertexData();			
		}
	}
	else if (strcmp(argv[2], "sssp") == 0){
		myapp = new SSSP();
		myapp->mytype = sssp;		
		for(int v=0; v<dag->num_vertices; v++){
			(*dag->g)[v].data = new SSSPVertexData();
			(*dag->g)[v].accum = new SSSPVertexData();
			dynamic_cast<SSSPVertexData*>((*dag->g)[v].data)->routed_nodes.push_back(v);
		}
		
		/*the source node to find the path*/
		if (dag->graph_type == test){
		//if(dag->graph_type == synthetic) {
			myapp->sources.push_back(0);  //0 is set as the source vertex of sssp
			printf("sourceid=00000\n");
		}
		else{
			/*select the one with the largest degree as source*/
			int degree = 0, sourceid = 0;
			for (int i = 0; i < dag->g->m_vertices.size(); i++){
				int curr_degree = dag->get_in_nbrs(i).size() + dag->get_out_nbrs(i).size();
				if (curr_degree > degree){
					degree = curr_degree;
					sourceid = i;
				}
			}
			myapp->sources.push_back(sourceid);
			printf("sourceid=%d\n", sourceid);
		}
		for (int i = 0; i < myapp->sources.size(); i++){
			dynamic_cast<SSSPVertexData*>((*dag->g)[myapp->sources[i]].data)->dist = 0;
			std::cout << "sources_dist=" << dynamic_cast<SSSPVertexData*>((*dag->g)[myapp->sources[i]].data)->dist << std::endl;;
		} 
	}
	else if(strcmp(argv[2], "subgraph") == 0){
		myapp = new SubgraphIsom();
		myapp->mytype = subgraph;
		for(int v=0; v<dag->num_vertices; v++){
			(*dag->g)[v].data = new SubgraphVertexData();
			(*dag->g)[v].accum = new SubgraphVertexData();
		}
		/* the pattern to be matched, no selfcycle
		*  without loss of generosity, we consider structural match 
		*  criteria: if vg matches vp, 
		*			 the adjacent edges of vg must be the superset of vp’s adjacent edges
		*/
		myapp->pattern_graph = new Graph();
		int num_vertex = 3;		
		for (int i = 0; i < num_vertex; i++){
			MyVertex* v = new MyVertex();
			v->vertex_id = i;//starting from 0	
			add_vertex(*v, *myapp->pattern_graph);
		}		
		if(num_vertex == 2){
			//pipeline (pattern size = 2)
			add_edge(0, 1, 10, *myapp->pattern_graph);
		}else if(num_vertex == 3){
			//triangle (pattern size = 3)
			add_edge(0, 1, 10, *myapp->pattern_graph);
			add_edge(0, 2, 10, *myapp->pattern_graph);
		}else if(num_vertex == 4){
			//square (pattern size = 4)
			add_edge(0, 1, 10, *myapp->pattern_graph);
			add_edge(0, 2, 10, *myapp->pattern_graph);
			add_edge(1, 3, 10, *myapp->pattern_graph);
		}
	}
	
	myapp->budget = budget;
	myapp->ITERATIONS = Max_Iteration;
	myapp->global_graph = dag;
	
	//apply the replication plan to the vertices
	Optimizer* optimizer = new Optimizer();
	optimizer->myapp = myapp;
	optimizer->dag = dag;
	optimizer->DCs = DCs;

	EngineType type = synchronous;
	GraphEngine* engine = new GraphEngine(type,DCs.size());
	engine->myapp = myapp;
	engine->myopt = optimizer;  //ribo,新加
	engine->DCs = DCs;
	engine->Simulate(argv[3]);	 //count the pagerank
	return 0;
}

