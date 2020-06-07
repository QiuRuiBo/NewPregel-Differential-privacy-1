#include "Simulator.h"
#include "Utils.h"
#include <omp.h>
#include <string>
#include <algorithm> 
#include <iostream>
#include <fstream>   //Output the results to the document
#define INI_BUDGET 17
#define UPDATE_SEND 1 
   
//#define TIMES (float)1 * pow(10,-5)
//#define TIMES (float)0
#define TIMES 100*3
extern int N_THREADS;
extern float noise_budget;
extern bool privacy;
  
void GraphEngine::Simulate(char* sizeofgraph){

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

    //identify the most highly connected vertex and trace its rank value
    int highest = 0;
    int hub_index = -1;
    for(int vi=0; vi < total_vertices; vi++){
        std::vector<int> in_nbrs = myapp->global_graph->get_in_nbrs(vi);
        if(in_nbrs.size() > highest){
            highest = in_nbrs.size();
            hub_index = vi;
        }
    }
    printf("hub vertex is %d\n",hub_index);
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
			if (myapp->mytype == pagerank)
			{
				dynamic_cast<PageRankVertexData*>((*myapp->global_graph->g)[vid].data)->rank = 0;
			}
			else if (myapp->mytype == sssp)
			{
				dynamic_cast<SSSPVertexData*>((*myapp->global_graph->g)[vid].data)->dist = std::numeric_limits<float>::max();

			}
			else if (myapp->mytype == subgraph)
			{
				//subgraph
			}
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


	if (myapp->mytype == sssp)
	{
		for (int i = 0; i < myapp->sources.size(); i++) {
			dynamic_cast<SSSPVertexData*>((*myapp->global_graph->g)[myapp->sources[i]].data)->dist = 0;
			//std::cout << "sources=" << dynamic_cast<SSSPVertexData*>((*myapp->global_graph->g)[myapp->sources[i]].data)->dist << std::endl;;
		}
	}

    /** Start the engine, each thread executes at the same time */
    float t_step = 0.0;
    int iter_counter = 0;
    float totalcost = 0.0;
    float wan_usage = 0.0;
	std::vector<float> tr_send_step = std::vector<float>(num_threads, 0.0);
	std::vector<float> tr_rcv_step = std::vector<float>(num_threads, 0.0);

	// LO: define a new counter
	int newCounter = 0, times = 0,times_ = 0;
	bool communicate = false;
	float sum_noise = 0;
    /**
    * Start the engine
    * currently only implemented the synchronous engine
    */
    if(type == synchronous){   //synchronous=0
        if(myapp->ITERATIONS != 0){
            /* converge within the fixed number of iterations */
            iter_counter = myapp->ITERATIONS;
			vector<float> rank_;
			//int ri_test = 0,a;
			while (iter_counter > 0) {

				/**
				* Execute Compute on all active vertices
				* Sync before sending messages
				*/
				std::cout << "-------------- Compute stage of iteration " << myapp->ITERATIONS - iter_counter << " -------------" << std::endl;
#pragma omp parallel for
				for (int tr = 0; tr < num_threads; tr++) {
#pragma omp parallel for
					for (int v = 0; v < Threads[tr]->l_dag->num_vertices; v++) {
						int vid = (*Threads[tr]->l_dag->g)[v].vertex_id;
						myapp->Compute(vid, myapp->global_graph);
					}
				}
				/**
				* Data aggregation
				*/
				std::cout << "-------------- Data aggregation stage  " << myapp->ITERATIONS - iter_counter << " -------------" << std::endl;
				for (int tr = 0; tr < num_threads; tr++) {
					if (DCs[tr]->aggregators.size() > 0) {
						for (int ai = 0; ai < DCs[tr]->aggregators.size(); ai++) {
							delete DCs[tr]->aggregators[ai];
						}
						DCs[tr]->aggregators.clear();
					}
				}
				for (int tr = 0; tr < num_threads; tr++) {  //遍历所有DC
					for (int v = 0; v < Threads[tr]->l_dag->num_vertices; v++) {  //遍历DC的所有vertices
						int vid = (*Threads[tr]->l_dag->g)[v].vertex_id;
						std::vector<int> out_nbrs = myapp->global_graph->get_out_nbrs(vid);
						float msg_rank;
						if (myapp->mytype == pagerank)
							msg_rank = dynamic_cast<PageRankVertexData*>((*myapp->global_graph->g)[vid].data)->rank / (float)out_nbrs.size();
						else if (myapp->mytype == sssp)
						{
							msg_rank = dynamic_cast<SSSPVertexData*>((*myapp->global_graph->g)[vid].data)->dist;
						}
						else if (myapp->mytype == subgraph)
						{
							//subgraph
						}
						for (int nid = 0; nid < out_nbrs.size(); nid++) {  //遍历顶点vid的所有出边的顶点，如果
							VertexData* l_accum;
							float weight;
							//if(myapp->mytype == pagerank)
							{
								if (myapp->mytype == pagerank)
								{
									l_accum = new PageRankVertexData();
									dynamic_cast<PageRankVertexData*>(l_accum)->rank = msg_rank;
								}
								else if (myapp->mytype == sssp)
								{
									l_accum = new SSSPVertexData();
									std::pair<Edge, bool> ed = boost::edge(vid, out_nbrs[nid], (*myapp->global_graph->g));
									weight = boost::get(boost::edge_weight_t(), (*myapp->global_graph->g), ed.first);
									dynamic_cast<SSSPVertexData*>(l_accum)->dist = msg_rank + weight;
								}
								else if (myapp->mytype == subgraph)
								{
									//subgraph
								}
								 
								int other_dc = (*myapp->global_graph->g)[out_nbrs[nid]].location_id;
								if (other_dc != tr) {  //如果出边的顶点不是当前DC内部，则归类到其所属的aggregator
									//search for the aggregator, if does not exist, create one
									bool foundagg = false;
									int founditer = -1;

									//(*myapp->global_graph->g)[out_nbrs[nid]].messages.push_back(l_accum);  //即使是gateway顶点也直接push back消息，用于原图的PageRank计算

									for (int ai = 0; ai < DCs[tr]->aggregators.size(); ai++) {
										if (DCs[tr]->aggregators[ai]->DC_id == other_dc && DCs[tr]->aggregators[ai]->vertex_id == out_nbrs[nid] 
											&& DCs[tr]->aggregators[ai]->source_vertex_id == vid) {  //一条msg就创建一个agg
											foundagg = true;
											founditer = ai;
											break;
										}
									}  
									if (foundagg) { 
										//aggregate to founditer
										if (myapp->mytype == pagerank)
										{
											dynamic_cast<PageRankVertexData*>(DCs[tr]->aggregators[founditer]->aggregated_data)->rank += msg_rank;
											DCs[tr]->aggregators[founditer]->v_list.push_back(out_nbrs[nid]);
										}
										else if (myapp->mytype == sssp)
										{
											
											//如果agg中的消息为无穷大且新收集进来的消息不是无穷大，则收集器的值赋值为新的消息
											if (msg_rank != std::numeric_limits<float>::max()
												&& dynamic_cast<SSSPVertexData*>(DCs[tr]->aggregators[founditer]->aggregated_data)->dist == std::numeric_limits<float>::max())
											{
												DCs[tr]->aggregators[founditer]->v_list.push_back(out_nbrs[nid]);  //放回其所属的agg
												std::pair<Edge, bool> ed = boost::edge(vid, out_nbrs[nid], (*myapp->global_graph->g));
												weight = boost::get(boost::edge_weight_t(), (*myapp->global_graph->g), ed.first);
												//dynamic_cast<SSSPVertexData*>(l_accum)->dist = msg_rank + weight;
												//std::cout << vid << "->" << out_nbrs[nid] << " 's weight_被放弃agg顶点=" << weight << std::endl;
												dynamic_cast<SSSPVertexData*>(DCs[tr]->aggregators[founditer]->aggregated_data)->dist = msg_rank + weight;
											}
											else if (msg_rank != std::numeric_limits<float>::max()
												&& dynamic_cast<SSSPVertexData*>(DCs[tr]->aggregators[founditer]->aggregated_data)->dist != std::numeric_limits<float>::max())
											{
												DCs[tr]->aggregators[founditer]->v_list.push_back(out_nbrs[nid]);  //放回其所属的agg
												std::pair<Edge, bool> ed = boost::edge(vid, out_nbrs[nid], (*myapp->global_graph->g));
												weight = boost::get(boost::edge_weight_t(), (*myapp->global_graph->g), ed.first);
												//dynamic_cast<SSSPVertexData*>(l_accum)->dist = msg_rank + weight;
												//std::cout << vid << "->" << out_nbrs[nid] << " 's weight_被放弃agg顶点=" << weight << std::endl;
												//msg_rank += weight;
												dynamic_cast<SSSPVertexData*>(DCs[tr]->aggregators[founditer]->aggregated_data)->dist += msg_rank + weight;
											}
										}
										else if (myapp->mytype == subgraph)
										{

										}
									}
									else {
										//create one
										Aggregator* new_agg = new Aggregator();
										new_agg->DC_id = other_dc;
										new_agg->vertex_id = out_nbrs[nid];
										new_agg->source_vertex_id = vid;
										new_agg->v_list.push_back(out_nbrs[nid]);
										new_agg->aggregated_data = l_accum;
										DCs[tr]->aggregators.push_back(new_agg);
										//printf("create one aggregator in DC %d\n",tr);
									}
								}
								else {
									//simply send the message
									(*myapp->global_graph->g)[out_nbrs[nid]].messages.push_back(l_accum);
								}
							}
							/*
						else if(myapp->mytype == sssp){
								l_accum = new SSSPVertexData();
							}else if(myapp->mytype == subgraph){
								l_accum = new SubgraphVertexData();
							}
							*/
						}
					}
				}
				/*Barrier*/
				/*
				for(int tr=0;tr<num_threads;tr++)
					for (int ai = 0; ai < DCs[tr]->aggregators.size(); ai++)
					{
						std::cout << "DC " << tr << " has " << DCs[tr]->aggregators[ai]->v_list.size() << " agg. " << std::endl;
					}
				*/



				//After data aggregation, add noises to the aggregated data and send out
				int agg_count = 0;
				for (int tr = 0; tr < num_threads; tr++) {
					agg_count += DCs[tr]->aggregators.size();
				}
				printf("%d aggregators in total.\n", agg_count);

				// std::cout << "THE NEWCOUNTER IS " << newCounter << std::endl;
				// LO: when newCounter is even numbers, send messages
				for (int tr = 0; tr < num_threads; tr++) {
					/**
					 * Imaging the largest diff neighboring graphs (containing only vertices in local DC)
					 * G: a node k <- all the other nodes but v, and k -> v, a random node n -> v
					 * G': edge k -> v is removed from G
					 *Marked by Ribo:It is noted here that we consider the worst case of each iteration in terms of sensitivity.
					 *In fact, further discussion can be conducted here to determine more accurate sensitivity
					 */
					float deltaf = 0; //the same for all aggregators in the same DC
					/*The initial PageRank value for each vertex is 0.15*/
					float rank_hub = 0.15 + 0.85 * (0.15 * (Threads[tr]->l_dag->num_vertices - 3) + 0.15 / 2.0);
					//float rank_hub = 0.15 + 0.85 * (0.15 * (645 - 3) + 0.15 / 2.0);
					//the gateway node
					float fg = 0.15 + 0.85 * (rank_hub + 0.15 / 2.0);
					float fg1 = 0.15 + 0.85 * 0.15 / 2.0;
					if (iter_counter < myapp->ITERATIONS)
						deltaf = std::abs(fg - fg1);  //ribo

					printf("Gf is: %f\n", deltaf);


					if (newCounter % UPDATE_SEND == 0 || times_ > TIMES) {   //只在特定的迭代次数发送消息,确定msg_last序列消息已清空
					//if ((newCounter % UPDATE_SEND == 0 || communicate) ) {
						times = 1;
						for (int ai = 0; ai < DCs[tr]->aggregators.size(); ai++) {
							//if (newCounter % UPDATE_SEND == 0 || communicate) {
							if (newCounter % UPDATE_SEND == 0 || times_ > TIMES) {
								//printf("正在清空记录的消息队列.\n");
								for (int nid = 0; nid < DCs[tr]->aggregators[ai]->v_list.size(); nid++) {
									int other_vertex = DCs[tr]->aggregators[ai]->v_list[nid];
									for (int i = 0; i < (*myapp->global_graph->g)[other_vertex].messages_last.size(); i++)
										(*myapp->global_graph->g)[other_vertex].messages_last[i] = 0;
									(*myapp->global_graph->g)[other_vertex].messages_last.clear();
								}
							}
							/*
							*斐波那契数列方式分配budget
							*/
							//float sum = utils::feibonaqi_sum(myapp->ITERATIONS);
							//std::cout << "sum= " << sum << std::endl;
							//std::cout << "当前斐波那契数列项数为:" << myapp->ITERATIONS - iter_counter + 1<< std::endl;
							//float percent = utils::feibonaqi_n(0,0,(float)(myapp->ITERATIONS - iter_counter+1)) / sum;
							//DCs[tr]->aggregators[ai]->noise_budget = noise_budget * percent /DCs[tr]->budget_sum_weight;


							if (myapp->mytype == sssp) {
								DCs[tr]->aggregators[ai]->noise_budget = noise_budget / (float)agg_count / (float)myapp->ITERATIONS;
								//std::cout << "DC " << tr << "'s budget_percentage=" << DCs[tr]->budget_weight_percentage << std::endl;
								//std::cout << "DC " << tr << "'s budget_sum_weight=" << DCs[tr]->budget_sum_weight << std::endl;
								//std::cout << "DC " << tr << "'s budget_weight=" << DCs[tr]->budget_weight << std::endl;
							}
							
							
							/*
							*修正的指数分配budget
							*/
							if (myapp->mytype == pagerank) {
								DCs[tr]->aggregators[ai]->noise_budget = (float)((noise_budget / (float)agg_count / (float)pow(2, INI_BUDGET - (float)((float)(myapp->ITERATIONS - iter_counter + 1) / myapp->ITERATIONS) * INI_BUDGET)));
								float Correction_coefficient = (float)(1 - (float)pow(2, -INI_BUDGET / ((float)myapp->ITERATIONS)));
								float k = (float)pow(2, INI_BUDGET) / (float)((float)pow(2, INI_BUDGET) - 1);
								Correction_coefficient = Correction_coefficient * k;
								DCs[tr]->aggregators[ai]->noise_budget = DCs[tr]->aggregators[ai]->noise_budget * Correction_coefficient;
							} 

							//std::cout << "DC " << tr <<"'s budget=" << DCs[tr]->aggregators[ai]->noise_budget << std::endl;
							if (DCs[tr]->aggregators[ai]->noise_budget > TIMES)
								communicate = true;

							if (privacy)//generate noise value 
							{
								float noise;
								noise = myapp->global_graph->laplace_generator(0, deltaf / DCs[tr]->aggregators[ai]->noise_budget);
								// LO: add noise
								int otherId = DCs[tr]->aggregators[ai]->DC_id;
								//std::cout << "DC " << DCs[tr]->pri_rank << " send message to DC " << DCs[otherId]->pri_rank << std::endl;
								if (DCs[tr]->pri_rank <= DCs[otherId]->pri_rank)  //低等级DC向高等级DC发送消息时不需要加noise
									noise = 0;
								else {
									sum_noise += DCs[tr]->aggregators[ai]->noise_budget;
								}
								if (myapp->mytype == pagerank)
								{
									dynamic_cast<PageRankVertexData*>(DCs[tr]->aggregators[ai]->aggregated_data)->rank += noise;
								}
								else if (myapp->mytype == sssp)
								{
									dynamic_cast<SSSPVertexData*>(DCs[tr]->aggregators[ai]->aggregated_data)->dist += noise;
								}
								else if (myapp->mytype == subgraph)
								{
									//subgraph
								}
								//std::cout << "noise=" << noise << std::endl;
							}

							for (int nid = 0; nid < DCs[tr]->aggregators[ai]->v_list.size(); nid++) {
								int other_vertex = DCs[tr]->aggregators[ai]->v_list[nid];
								VertexData* newmsg;
								if (myapp->mytype == pagerank)
								{
									newmsg = new PageRankVertexData();
								}
								else if (myapp->mytype == sssp)
								{
									newmsg = new SSSPVertexData();

								}
								else if (myapp->mytype == subgraph)
								{

								}
								//msg队列记录的应该是PR(j)/L(j),PR(j)为网页j的PR值，L(j)是网页j的链出网页数

								if (myapp->mytype == pagerank)
								{
									dynamic_cast<PageRankVertexData*>(newmsg)->rank = dynamic_cast<PageRankVertexData*>(DCs[tr]->aggregators[ai]->aggregated_data)->rank / (float)DCs[tr]->aggregators[ai]->v_list.size();
								}
								else if (myapp->mytype == sssp)
								{
									dynamic_cast<SSSPVertexData*>(newmsg)->dist = dynamic_cast<SSSPVertexData*>(DCs[tr]->aggregators[ai]->aggregated_data)->dist / (float)DCs[tr]->aggregators[ai]->v_list.size();

								}
								else if (myapp->mytype == subgraph)
								{

								}
								(*myapp->global_graph->g)[other_vertex].messages.push_back(newmsg);
								(*myapp->global_graph->g)[other_vertex].messages_last.push_back(newmsg);  //messages_temp记录上次的message值

								wan_usage += myapp->msg_size;  //记录通信过程的wan usage
									//pregel模型一次通信cost只需计算发送方DC的上传费用以及接收方的下载价格
								totalcost += myapp->msg_size * DCs[Threads[tr]->DC_loc_id]->upload_price / 1000.0; //per GB,记录通信过程的上传花费
								totalcost += myapp->msg_size * DCs[Threads[tr]->DC_loc_id]->download_price / 1000.0; //per GB,记录通信过程的上传花费
								//耗时t的计算不同于GAS model,pregel model仅需发送消息一次（GAS两次）,时间仅由上传 or 下载最大的那个决定
								tr_send_step[tr] += myapp->msg_size / DCs[Threads[tr]->DC_loc_id]->upload_band;
								tr_rcv_step[tr] += myapp->msg_size / DCs[Threads[tr]->DC_loc_id]->download_band;
								
							}
						}
						times_++;

					}
					else {  //如果不是指定的迭代次数，则保留着上次的message值
						times = 0;
						for (int ai = 0; ai < DCs[tr]->aggregators.size(); ai++) {
							std::vector<int> hasRecover;
							for (int nid = 0; nid < DCs[tr]->aggregators[ai]->v_list.size(); nid++) {
								int other_vertex = DCs[tr]->aggregators[ai]->v_list[nid];

								std::vector<int>::iterator it;
								it = find(hasRecover.begin(), hasRecover.end(), other_vertex);
								if (!(it != hasRecover.end()))  //other_vertex不在已恢复队列中，说明还没恢复该消息队列的值
								{
									//std::cout << "正在恢复" << std::endl;
									hasRecover.push_back(other_vertex);
									for (int i = 0; i < (*myapp->global_graph->g)[other_vertex].messages_last.size(); i++) {
										(*myapp->global_graph->g)[other_vertex].messages.push_back((*myapp->global_graph->g)[other_vertex].messages_last[i]);
									}
								}
							}
							hasRecover.clear();  //清空已恢复记录队列
						}

					}
				}
				//一个DC处理完毕，计算本轮耗时。t=max_time(DC0,DC1,DC2)
				float time = 0.0;
				for (int dc = 0; dc < num_threads; dc++) {
					//printf("tr_rcv_step=%f\n", tr_rcv_step[dc]);
					time = time > tr_rcv_step[dc] ? time : tr_rcv_step[dc];
					tr_rcv_step[dc] = 0;  //取完顺便清空
					//printf("tr_send_step=%f\n", tr_send_step[dc]);
					time = time > tr_send_step[dc] ? time : tr_send_step[dc];
					tr_send_step[dc] = 0;  //取完顺便清空

				}
				//printf("time=%f\n", time);
				t_step += time;
				//std::cout << "一轮结束" << std::endl;

				// LO: newCounter++
				newCounter++;
				if (times)
					iter_counter--;
				if (myapp->mytype == pagerank)
				{
					if (hub_index != -1)
						printf("hub vertex rank value: %.4f\n", dynamic_cast<PageRankVertexData*>((*myapp->global_graph->g)[hub_index].data)->rank);
				float avg_rank = 0;
				for (int vi = 0; vi < total_vertices; vi++) {
					avg_rank += dynamic_cast<PageRankVertexData*>((*myapp->global_graph->g)[vi].data)->rank;
				}
				avg_rank /= (float)total_vertices;
				printf("iter %d: average rank %.4f\n", myapp->ITERATIONS - iter_counter, avg_rank);
				//if (newCounter % UPDATE_SEND == 0)  //只记录发送noise的迭代次数的average rank
				//if (newCounter % UPDATE_SEND == 0 || DCs[0]->aggregators[0]->noise_budget > TIMES)
				if (times)
					rank_.push_back(avg_rank);
				}
            }

			printf("rank值是：\n");
			for(int i=0;i<rank_.size();i++)
				std::cout << rank_[i] << std::endl;
        } 
		else 
		{
            /* converge according to the tolerance */
            iter_counter ++;
			newCounter--;
        }
		std::cout << "sum budget=" << sum_noise << std::endl;
    }//end of pagerank
	
	else if(type == asynchronous){
        std::cout <<"asynchronous engine is not implemented. " << std::endl;
    }

    for(int tr=0; tr<num_threads; tr++) {
        printf("DC %d has %d aggregators\n",tr,(int)DCs[tr]->aggregators.size());
    }

    
	for(int tr=0; tr<num_threads; tr++){
		printf("DC %d has %d aggregators\n",tr,(int)DCs[tr]->aggregators.size());
	}

	map<int, vector<float>> dc0;
	if (myapp->mytype == pagerank) {
		for (int i = 0; i < num_threads; i++) {
			// std::cout << "DC: " << i << std::endl;
			for (int v = 0; v < Threads[i]->l_dag->num_vertices; v++) {
				//output the pagerank value
				dc0[i].push_back(dynamic_cast<PageRankVertexData*>((*Threads[i]->l_dag->g)[v].data)->rank);
				//printf("%f\n", dynamic_cast<PageRankVertexData*>((*Threads[i]->l_dag->g)[v].data)->rank);
			}
		}
	}


	int size = atoi(sizeofgraph);
	std::cout << "size=" << size << std::endl;
	if (myapp->mytype == sssp) {

		std::vector<int> sourcesList;
		for (int i = 0; i < myapp->sources.size(); i++) {
			sourcesList.push_back(myapp->sources[i]);  //注意该ID为全局ID，后续使用先转为local ID
		}


		for (int i = 0; i < num_threads; i++) {
			for (int v = 0; v < Threads[i]->l_dag->num_vertices; v++) {
				bool is_sources = false;
				int vid = (*Threads[i]->l_dag->g)[v].vertex_id;
				for (int i = 0; i < myapp->sources.size(); i++)
					if (vid == myapp->sources[i])
					{
						is_sources = true;
						break;
					}
				if (is_sources == false)
				{
					dc0[i].push_back(dynamic_cast<SSSPVertexData*>((*Threads[i]->l_dag->g)[v].data)->dist);
					//std::cout << "从0出发到"<< v << "点的最短距离为：" <<dynamic_cast<SSSPVertexData*>((*myapp->global_graph->g)[v].data)->dist << std::endl;
					//std::cout << dynamic_cast<SSSPVertexData*>((*Threads[i]->l_dag->g)[v].data)->dist << std::endl;

				}
			}
		}
	} 


	for (int tr = 0; tr < num_threads; tr++)
	{
		std::cout << "#################################################" << std::endl
			<< tr << "'s upload_band=" << DCs[Threads[tr]->DC_loc_id]->upload_band << std::endl
			<< "#################################################" << std::endl
			<< tr << "'s download_band=" << DCs[Threads[tr]->DC_loc_id]->download_band << std::endl;

	}


	std::cout << "wan_usage=" << wan_usage << ", total_cost=" << totalcost << ", t_step=" << t_step << std::endl;
	std::cout << "msg_size=" << myapp->msg_size << std::endl;


	/*output the realtive error*/
	/*The 'datafile' should store
	**the PageRank of the vertices
	**in the DC0 in the original graph
	*/
	string datafile;
	if (myapp->mytype == pagerank)
		datafile = "pagerank_";
	else if (myapp->mytype == sssp)
		datafile = "sssp_";
	datafile += sizeofgraph;
	datafile += ".txt";
	//std::cout << "file name:" << datafile << std::endl;
	std::ifstream non_privacy(datafile);
	if (non_privacy.is_open() == false)
		std::cout << "file " << datafile << " is not exist." << std::endl;
	else {
		std::string rndline = "0";
		int index_to_pr = 0;
		float loc;
		float totalCount = 0;
		float errorSum = 0;

		float totalCount_global = 0;
		float errorSum_global = 0;

		for (int tr = 0; tr < num_threads; tr++)
		{
			while (1) {
				if (rndline.empty() == true || index_to_pr >= dc0[tr].size())
					//if (rndline.empty() == true)
						//if (rndline.empty() == true)
					break;
				std::getline(non_privacy, rndline);
				std::stringstream strm(rndline);  //stringstream:Convert a string to a number
				strm >> loc;
				//printf("loc=%f\n", loc);
				//if (dc0[index_to_pr] != std::numeric_limits<float>::max()) 
				{
					errorSum += abs(dc0[tr][index_to_pr] - loc) / loc;
					totalCount++;

					errorSum_global += abs(dc0[tr][index_to_pr] - loc) / loc;
					totalCount_global++;
				}
				index_to_pr++;
			}
			double average_relative_error = errorSum / (float)totalCount;
			std::cout << "DC " << tr << "'s num_vertices = " << totalCount << std::endl;
			std::cout << "DC " << tr << "'s Average relative error = " << average_relative_error << std::endl;
			totalCount = 0;
			errorSum = 0;
			index_to_pr = 0;
		}
		double average_relative_error_global = errorSum_global / (float)totalCount_global;
		std::cout << "Total num_vertices = " << totalCount_global << std::endl;
		std::cout << "Average relative error = " << average_relative_error_global << std::endl;

	}
}
