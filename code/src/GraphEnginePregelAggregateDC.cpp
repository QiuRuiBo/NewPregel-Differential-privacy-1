#include "Simulator.h"
#include "Utils.h"
#include <omp.h>
#include <string>
#include <algorithm> 
#include <iostream>
#include <fstream>   //Output the results to the document
#define INI_BUDGET 17  
#define UPDATE_SEND 1
#define AGG_NUM 1    //DC0发向DC1的agg数量 
#define ABANDON 1      
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
				//std::cout << "Threads=" << (*Threads[tr]->l_dag->g)[v].vertex_id << "Global=" << (*myapp->global_graph->g)[v].vertex_id;
				//int a;
				//scanf("%d", &a);
				/*
				Graph* pattern_graph;
				for (int i = 0; i < pattern_graph->myvertex.size(); i++) {
					//no label, so add all
					SubgraphVertexData::Message init_msg;
					init_msg.pairs.push_back(std::pair<int, int>(i, vid));
					init_msg.forwarding_trace.push_back(vid);
					dynamic_cast<SubgraphVertexData*>((*myapp->global_graph->g)[vid].data)->matches.push_back(init_msg);
				
				}
				*/
				
			}

			(*myapp->global_graph->g)[vid].status = activated;  //初始化所有顶点为activate												//identify gateway nodes
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

	//设置sssp算法的源顶点距离为0
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
	std::vector<int> abandoned_record;
	vector<float> budget_iter;
	vector<float> budget_dc;
	int abandoned_record_temp = 0;
	//float wan_usage = 0.0, cost = 0, performance_time = 0;
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
				int activated_num = 0;
				std::cout << "-------------- Compute stage of iteration " << myapp->ITERATIONS - iter_counter << " -------------" << std::endl;
				#pragma omp parallel for
				for (int tr = 0; tr < num_threads; tr++) {
					#pragma omp parallel for
					for (int v = 0; v < Threads[tr]->l_dag->num_vertices; v++) {
						int vid = (*Threads[tr]->l_dag->g)[v].vertex_id;
						myapp->Compute(vid, myapp->global_graph);
						if ((*myapp->global_graph->g)[vid].status == activated)
							activated_num++;
					}
				}
				printf("本轮activate顶点数量为：%d\n", activated_num);
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

				int num_vertices_all = 0;  //记录当前测试图总共的顶点个数
				for (int tr = 0; tr < num_threads; tr++) {  //遍历所有DC
					for (int v = 0; v < Threads[tr]->l_dag->num_vertices; v++) {  //遍历DC的所有vertices
						int vid = (*Threads[tr]->l_dag->g)[v].vertex_id;
						num_vertices_all++;
					}
				}

				for (int tr = 0; tr < num_threads; tr++) {  //遍历所有DC
					float max_rank = 0, min_rank = 0, dis_rank = 0;


					/*
					**找到当前DC所有顶点的最大、最小rank值
					*/
					for (int v = 0; v < Threads[tr]->l_dag->num_vertices; v++) {  //遍历DC的所有vertices
						int vid = (*Threads[tr]->l_dag->g)[v].vertex_id;
						std::vector<int> out_nbrs = myapp->global_graph->get_out_nbrs(vid);
						(*myapp->global_graph->g)[vid].is_urgently = false;  //初始化urgently程度为false
						float msg_rank;
						if (myapp->mytype == pagerank)
							msg_rank = dynamic_cast<PageRankVertexData*>((*myapp->global_graph->g)[vid].data)->rank / (float)out_nbrs.size();
						else if (myapp->mytype == sssp)
							msg_rank = dynamic_cast<SSSPVertexData*>((*myapp->global_graph->g)[vid].data)->dist;
						else if (myapp->mytype == subgraph)
						{
							//subgraph
						}
						if (msg_rank != std::numeric_limits<float>::max())
							max_rank = max_rank > msg_rank ? max_rank : msg_rank;
						min_rank = min_rank < msg_rank ? min_rank : msg_rank;
					}
					dis_rank = max_rank - min_rank;


					/*
					**创建好指定数量的agg
					*/
					for (int i = 0; i < AGG_NUM; i++) {
						//create
						for (int j = 0; j < num_threads; j++)
						{
							if (num_threads - j - 1 != tr) {
								Aggregator* new_agg = new Aggregator();
								new_agg->DC_id = num_threads - j - 1;
								new_agg->vertex_id = AGG_NUM - i - 1;
								VertexData* l_accum;
								if (myapp->mytype == pagerank)
									l_accum = new PageRankVertexData();
								else if (myapp->mytype == sssp)
									l_accum = new SSSPVertexData();
								else if (myapp->mytype == subgraph)
									l_accum = new SubgraphVertexData();  //subgraph

								new_agg->aggregated_data = l_accum;
								DCs[tr]->aggregators.push_back(new_agg);
							}
						}
					}

					int aban = 0, save = 0;
					for (int v = 0; v < Threads[tr]->l_dag->num_vertices; v++) {  //遍历DC的所有vertices
						int vid = (*Threads[tr]->l_dag->g)[v].vertex_id;
						std::vector<int> out_nbrs = myapp->global_graph->get_out_nbrs(vid);
						int agg_belong = 0;
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
						//std::cout << "min=" << min_rank << std::endl << "max=" << max_rank << std::endl << "步进=" << dis_rank / AGG_NUM << std::endl;
						//std::cout << "msg_rank=" << msg_rank << std::endl;
						//判断当前vertex属于哪个aggregator
						for (int i = 0; i < AGG_NUM; i++)
						{
							if (abs(msg_rank) >= abs(min_rank) + (float)i * abs(dis_rank) / (float)AGG_NUM && abs(msg_rank) < abs(min_rank) + (float)(i + 1) * abs(dis_rank) / (float)AGG_NUM) {  //不管正负，都应该抛弃那些绝对值小的
								agg_belong = i;
								break;
							}
						}
						//std::cout << "agg_belong=" << agg_belong << std::endl;
						//int a = 0;
						//scanf("%d", &a);

						srand((unsigned)time(NULL));
						for (int nid = 0; nid < out_nbrs.size(); nid++) {  //遍历顶点vid的所有出边的顶点，如果
							VertexData* l_accum;
							float weight;
							//if ((*myapp->global_graph->g)[out_nbrs[nid]].status == activated || iter_counter < 3)  //只有顶点处于activated状态才会受理
							{
								//if (myapp->mytype == pagerank || myapp->mytype == sssp)
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
										//std::cout << vid << "->" << out_nbrs[nid] <<" 's weight_内部顶点=" << weight <<std::endl;
										//dynamic_cast<SSSPVertexData*>(l_accum)->dist = msg_rank + weight;
										dynamic_cast<SSSPVertexData*>(l_accum)->dist = msg_rank + weight;
									}
									else if (myapp->mytype == subgraph)
									{
										//subgraph
									}

									int other_dc = (*myapp->global_graph->g)[out_nbrs[nid]].location_id;
									if (other_dc != tr) {  //search for the aggregator
										bool foundagg = false;
										//int founditer = -1;
										for (int ai = 0; ai < DCs[tr]->aggregators.size(); ai++) {
											//std::cout << "other dc=" << other_dc << std::endl;
											if (DCs[tr]->aggregators[ai]->DC_id == other_dc && DCs[tr]->aggregators[ai]->vertex_id == agg_belong) {
												//std::cout << "找到归属AGG，其DC_ID=" << DCs[tr]->aggregators[ai]->DC_id << ",VERTEX_ID=" << DCs[tr]->aggregators[ai]->vertex_id << std::endl;
												if (DCs[tr]->aggregators[ai]->vertex_id < ABANDON)  //将被抛弃的rank值最小的agg，从中随机取出几个保留下来
												//if (DCs[tr]->aggregators[ai]->vertex_id > 4)  //将被抛弃的rank值最大的agg，从中随机取出几个保留下来
												{
													foundagg = true;
													/** 
													*按照概率从被抛弃的agg中选择一些点重新保留
													*/
													//float p = 0.1;
													/**
													*这可能与我们的直观想法有点不一样，直观上应该是rank值越小的顶点对结果影响越小，
													*越不应该恢复出去，但是这是基于顶点的rank值是等于或者接近于其真实rank值的情况
													*实际应用中，试想在该抛弃顶点的模型中，rank值越小，意味着其一直被抛弃，越应该放出来
													*让它参与计算，修正其他顶点以及它自己本身的rank值，这也是选择rank值小的放出来效果优于
													*rank值大的效果的原因。
													*/
													//float p = abs(1 / msg_rank);  //rank值越小的顶点越容易被恢复
													float p;  //距离初始rank值0.15越近，越容易被恢复
													float random;
													if (myapp->mytype == pagerank)
													{
														random = rand() % (100000 + 1) / (float)(100000 + 1);
														//p = abs(1 / (0.15 - msg_rank));
														//p /= 5;
														
														float temp; 
														temp = dynamic_cast<PageRankVertexData*>((*myapp->global_graph->g)[vid].data)->last_change;
														p = 1 / abs(temp - 0.001);
														p /= 50;
														
													}
													else if (myapp->mytype == sssp)
													{
														float temp;
														
														//temp = msg_rank;
														//if (msg_rank == std::numeric_limits<float>::max())
															//temp = 0;
														//p = abs(1 / (temp - 0.001));
														//p *= 10;
														
														temp = dynamic_cast<SSSPVertexData*>((*myapp->global_graph->g)[vid].data)->last_change;
														p = 1 / abs(temp-0.001);
														p /= 50;
														
														random = rand() % (100000 + 1) / (float)(100000 + 1);
													}

													/*
													if (random < p)
													{
														std::vector<int> in_nbrs = myapp->global_graph->get_in_nbrs(vid);
														for (int in_nbr = 0; in_nbr < in_nbrs.size(); in_nbr++) {
															(*myapp->global_graph->g)[in_nbrs[in_nbr]].is_urgently = true;
														}
													}
													*/
													 
													//float random = rand() % (1) / (float)(1);
													//p = 10000000;
													//if (random < p || (*myapp->global_graph->g)[vid].is_urgently == true)  //命中的顶点分配到agg上
													if (random < p)
													{
														//std::cout << "找到归属AGG，其DC_ID=" << DCs[tr]->aggregators[ai]->DC_id << ",VERTEX_ID=" << DCs[tr]->aggregators[ai]->vertex_id << std::endl;
														int index = 0;
														int index_backup = 0;
														//while (DCs[tr]->aggregators[index++]->vertex_id != ABANDON - 1);  //从rank值小到大找到最近的agg
														index = agg_belong;
														//index = 5;
														//while (DCs[tr]->aggregators[index--]->vertex_id != 3);  //从rank值大到小找到最近的agg
														/*
														//找到最近的不为空的agg
														index_backup = index;
														while (DCs[tr]->aggregators[index]->v_list.size() == 0)
														{
															index++;
															if (index >= DCs[tr]->aggregators.size())
															{
																index = index_backup;
																break;
															}
														}
														*/
														if (myapp->mytype == pagerank)
														{
															DCs[tr]->aggregators[index]->v_list.push_back(out_nbrs[nid]);  //放回其所属的agg
															dynamic_cast<PageRankVertexData*>(DCs[tr]->aggregators[index]->aggregated_data)->rank += msg_rank;
														}
														else if (myapp->mytype == sssp)
														{
															//如果agg中的消息为无穷大且新收集进来的消息不是无穷大，则收集器的值赋值为新的消息
															if (msg_rank != std::numeric_limits<float>::max()
																&& dynamic_cast<SSSPVertexData*>(DCs[tr]->aggregators[index]->aggregated_data)->dist == std::numeric_limits<float>::max())
															{
																DCs[tr]->aggregators[index]->v_list.push_back(out_nbrs[nid]);  //放回其所属的agg
																std::pair<Edge, bool> ed = boost::edge(vid, out_nbrs[nid], (*myapp->global_graph->g));
																weight = boost::get(boost::edge_weight_t(), (*myapp->global_graph->g), ed.first);
																//dynamic_cast<SSSPVertexData*>(l_accum)->dist = msg_rank + weight;
																//std::cout << vid << "->" << out_nbrs[nid] << " 's weight_被放弃agg顶点=" << weight << std::endl;
																dynamic_cast<SSSPVertexData*>(DCs[tr]->aggregators[index]->aggregated_data)->dist = msg_rank + weight;
															}
															else if (msg_rank != std::numeric_limits<float>::max()
																&& dynamic_cast<SSSPVertexData*>(DCs[tr]->aggregators[index]->aggregated_data)->dist != std::numeric_limits<float>::max())
															{
																DCs[tr]->aggregators[index]->v_list.push_back(out_nbrs[nid]);  //放回其所属的agg
																std::pair<Edge, bool> ed = boost::edge(vid, out_nbrs[nid], (*myapp->global_graph->g));
																weight = boost::get(boost::edge_weight_t(), (*myapp->global_graph->g), ed.first);
																//dynamic_cast<SSSPVertexData*>(l_accum)->dist = msg_rank + weight;
																//std::cout << vid << "->" << out_nbrs[nid] << " 's weight_被放弃agg顶点=" << weight << std::endl;
																//msg_rank += weight;
																dynamic_cast<SSSPVertexData*>(DCs[tr]->aggregators[index]->aggregated_data)->dist += msg_rank + weight;
															}
														}
														else if (myapp->mytype == subgraph)
														{
															DCs[tr]->aggregators[index]->v_list.push_back(out_nbrs[nid]);  //放回其所属的agg
															//subgraph
														}


														save++;
													}
													else
														aban++;
												}
												else if (DCs[tr]->aggregators[ai]->vertex_id >= ABANDON) {  //没有被抛弃的agg，正常处理
													//printf("正常处理.\n");
													foundagg = true;

													//std::cout << "找到归属AGG，其DC_ID=" << DCs[tr]->aggregators[ai]->DC_id << ",VERTEX_ID=" << DCs[tr]->aggregators[ai]->vertex_id << std::endl;
													if (myapp->mytype == pagerank)
													{
														dynamic_cast<PageRankVertexData*>(DCs[tr]->aggregators[ai]->aggregated_data)->rank += msg_rank;
													}
													else if (myapp->mytype == sssp)
													{
														std::pair<Edge, bool> ed = boost::edge(vid, out_nbrs[nid], (*myapp->global_graph->g));
														weight = boost::get(boost::edge_weight_t(), (*myapp->global_graph->g), ed.first);
														//dynamic_cast<SSSPVertexData*>(l_accum)->dist = msg_rank + weight;
														//std::cout << vid << "->" << out_nbrs[nid] << " 's weight_未被放弃agg顶点=" << weight << std::endl;
														//msg_rank += weight;
														dynamic_cast<SSSPVertexData*>(DCs[tr]->aggregators[ai]->aggregated_data)->dist += msg_rank + weight;
													}
													else if (myapp->mytype == subgraph)
													{
														//subgraph
													}

													DCs[tr]->aggregators[ai]->v_list.push_back(out_nbrs[nid]);
													break;
												}
											}

										}
										if (foundagg == false)
											printf("找不到agg.\n");
									}
									else {
										//simply send the message
										(*myapp->global_graph->g)[out_nbrs[nid]].messages.push_back(l_accum);
									}
								}
							}
							/*
							else if(myapp->mytype == sssp){
								l_accum = new SSSPVertexData();
							}
							else if(myapp->mytype == subgraph){
								l_accum = new SubgraphVertexData();
							}
							*/
						} 
					}

					for (int i = 0; i < DCs[tr]->aggregators.size(); i++)
					{
						std::wcout << "agg" << i << "'s size=" << DCs[tr]->aggregators[i]->v_list.size() << ",DC_ID=" << DCs[tr]->aggregators[i]->DC_id << ",VERTEX_ID=" << DCs[tr]->aggregators[i]->vertex_id << std::endl;
					}
					std::wcout << "本DC共抛弃" << aban << "个顶点," << "共保留下" << save << "个顶点" << std::endl;
					std::cout << "-------------------DC " << tr << "---------------------" << std::endl;
					abandoned_record_temp += aban;

				}

				abandoned_record.push_back(abandoned_record_temp);
				abandoned_record_temp = 0;

				/*Barrier*/

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

						/*按照相反的顺序输出agg的vlist.size*/
						/*
						int vlist_size_total = 0;
						//float dc_rank = 0;
						std::map<int, float> vlist_distribution;
						std::vector<int> list_;  // 存储各个agg的编号
						for (int i = 0; i < DCs[tr]->aggregators.size(); i++)
						{
							if (DCs[tr]->aggregators[i]->v_list.size() != 0)
							{
								list_.push_back(i);  // 存储非空agg编号
								vlist_size_total += DCs[tr]->aggregators[i]->v_list.size();  //记录整个DC的所有rank值，用来对不同agg分配不同的budget
								vlist_distribution[i] = DCs[tr]->aggregators[i]->v_list.size();  //记录每个agg的总rank
							}
						}

						while (list_.size() > 1) {
							int max_index = 0, min_index = 0;
							for (int i = 0; i < list_.size(); i++) {
								if (vlist_distribution[list_[i]] >= vlist_distribution[list_[max_index]])
									max_index = i;
								if (vlist_distribution[list_[i]] < vlist_distribution[list_[min_index]])
									min_index = i;
							}
							//exchange the max and min rank
							float temp = vlist_distribution[list_[min_index]];
							vlist_distribution[list_[min_index]] = vlist_distribution[list_[max_index]];
							vlist_distribution[list_[max_index]] = temp;
							//delete the index in list
							std::vector<int>::iterator it_max = list_.begin() + max_index, it_min = list_.begin() + min_index;
							if (min_index < max_index) {
								list_.erase(it_max);
								list_.erase(it_min);
							}
							else if (min_index > max_index)
							{
								list_.erase(it_min);
								list_.erase(it_max);
							}
						}
						*/


						/**
						*只有不为空的agg才需要消耗budget，这里记录需要消耗budget的agg数量
						*按照rank值大小的倒序存储agg的rank值
						*/
						/*
						float dc_rank = 0;
						std::map<int, float> rank_distribution;
						std::vector<int> list;  // 存储各个agg的编号
						for (int i = 0; i < DCs[tr]->aggregators.size(); i++)
						{
							if (DCs[tr]->aggregators[i]->v_list.size() != 0)
							{
								list.push_back(i);  // 存储非空agg编号
								dc_rank += abs(dynamic_cast<PageRankVertexData*>(DCs[tr]->aggregators[i]->aggregated_data)->rank);  //记录整个DC的所有rank值，用来对不同agg分配不同的budget
								rank_distribution[i] = abs(dynamic_cast<PageRankVertexData*>(DCs[tr]->aggregators[i]->aggregated_data)->rank);  //记录每个agg的总rank
							}
						}


						while (list.size() > 1) {
							int max_index = 0, min_index = 0;
							for (int i = 0; i < list.size(); i++) {
								if (rank_distribution[list[i]] >= rank_distribution[list[max_index]])
									max_index = i;
								if (rank_distribution[list[i]] < rank_distribution[list[min_index]])
									min_index = i;
							} 
							//exchange the max and min rank
							float temp = rank_distribution[list[min_index]];
							rank_distribution[list[min_index]] = rank_distribution[list[max_index]];
							rank_distribution[list[max_index]] = temp;
							//delete the index in list
							std::vector<int>::iterator it_max = list.begin() + max_index , it_min = list.begin() + min_index;
							if (min_index < max_index) {
								list.erase(it_max);
								list.erase(it_min);
							}
							else if(min_index > max_index)
							{
								list.erase(it_min);
								list.erase(it_max);
							}
						}
						*/
						int total_all = 0;
						for (int dc = 0; dc < num_threads; dc++) {
							for (int i = 0; i < DCs[dc]->aggregators.size(); i++)
							{
								if (DCs[dc]->pri_rank > DCs[DCs[dc]->aggregators[i]->DC_id]->pri_rank)  //高等级DC向低等级DC发送消息才需要budget
									if (DCs[dc]->aggregators[i]->v_list.size() != 0)
									{
										total_all++;
									}
							}
						}
						std::wcout << "total_all=" << total_all << std::endl;


						for (int ai = 0; ai < DCs[tr]->aggregators.size(); ai++) {
							if (DCs[tr]->aggregators[ai]->v_list.size() != 0) {
								//if (newCounter % UPDATE_SEND == 0 || communicate) {
								if (newCounter % UPDATE_SEND == 0 || times_ > TIMES) {
									printf("正在清空记录的消息队列.\n");
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
								/*
								float sum = utils::feibonaqi_sum(myapp->ITERATIONS);
								float percent = utils::feibonaqi_n(0, 0, (float)(myapp->ITERATIONS - iter_counter + 1)) / sum;
								DCs[tr]->aggregators[ai]->noise_budget = noise_budget * percent;
								*/

								/*
								*线性分配
								*/
								//DCs[tr]->aggregators[ai]->noise_budget = noise_budget / (float)iter_counter;
								//std::cout << "DC " << tr << "'s budget_percentage=" << DCs[tr]->budget_weight_percentage << std::endl;
								//std::cout << "DC " << tr << "'s budget_sum_weight=" << DCs[tr]->budget_sum_weight << std::endl;
								//std::cout << "DC " << tr << "'s budget_weight=" << DCs[tr]->budget_weight << std::endl;

								/*
								*SSSP
								*平均分配，sssp算法由于每一个迭代都很重要，因此只需要平均分配
								*/
								if (myapp->mytype == sssp)
								{
									DCs[tr]->aggregators[ai]->noise_budget = noise_budget / (float)myapp->ITERATIONS;
								}
								
								/**
								*Pagerank
								*修正的指数分配budget,这里分配的是每个iter的budget，还需具体分配到每次iter的agg上
								*/
								else if (myapp->mytype == pagerank)
								{
									DCs[tr]->aggregators[ai]->noise_budget = (float)((noise_budget / (float)pow(2, INI_BUDGET - (float)((float)(myapp->ITERATIONS - iter_counter + 1) / myapp->ITERATIONS) * INI_BUDGET)));  //固定首项的指数分配
									float Correction_coefficient = (float)(1 - (float)pow(2, -INI_BUDGET / ((float)myapp->ITERATIONS)));
									float k = (float)pow(2, INI_BUDGET) / (float)((float)pow(2, INI_BUDGET) - 1);
									Correction_coefficient = Correction_coefficient * k;
									DCs[tr]->aggregators[ai]->noise_budget = DCs[tr]->aggregators[ai]->noise_budget * Correction_coefficient;
								}
								//budget_iter.push_back(DCs[tr]->aggregators[ai]->noise_budget); //记录本轮迭代所分配到的总budget

								else if (myapp->mytype == subgraph) {
									//TODO
								}


								/*拿到本轮iter的总budget，继续分配budget给具体的agg*/
								//float percent = abs(1 / dynamic_cast<PageRankVertexData*>(DCs[tr]->aggregators[ai]->aggregated_data)->rank / dc_rank);
								//std::cout << "total rank=" << dynamic_cast<PageRankVertexData*>(DCs[tr]->aggregators[ai]->aggregated_data)->rank << std::endl;

								//DCs[tr]->aggregators[ai]->noise_budget = (float)((noise_budget * DCs[tr]->aggregators[ai]->v_list.size() / vlist_size_total / (float)pow(2, INI_BUDGET - (float)((float)(myapp->ITERATIONS - iter_counter + 1) / myapp->ITERATIONS) * INI_BUDGET)));  //按照agg的总vlist值大小不平等对待各个agg
								//DCs[tr]->aggregators[ai]->noise_budget = (float)((noise_budget * abs(dynamic_cast<PageRankVertexData*>(DCs[tr]->aggregators[ai]->aggregated_data)->rank) /dc_rank / (float)pow(2, INI_BUDGET - (float)((float)(myapp->ITERATIONS - iter_counter + 1) / myapp->ITERATIONS) * INI_BUDGET)));  //按照agg的总rank值大小不平等对待各个agg
								DCs[tr]->aggregators[ai]->noise_budget = DCs[tr]->aggregators[ai]->noise_budget / (float)total_all;  //平等对待各个aggregator
								//DCs[tr]->aggregators[ai]->noise_budget = (float)((noise_budget /DCs[tr]->budget_sum_weight / (float)pow(2, INI_BUDGET - (float)((float)(myapp->ITERATIONS - iter_counter + 1) / myapp->ITERATIONS) * INI_BUDGET)));
								//DCs[tr]->aggregators[ai]->noise_budget = (float)((noise_budget / (float)agg_count / (float)pow(2, INI_BUDGET - (float)((float)(myapp->ITERATIONS - iter_counter + 1) / myapp->ITERATIONS) * INI_BUDGET)));  //所有agg（不管有没有msg）都平等分配到budget
								/*
								if (tr == 0)
									budget_dc.push_back(DCs[tr]->aggregators[ai]->noise_budget);
								//agg_num[tr]
								*/
								std::cout << "budget=" << DCs[tr]->aggregators[ai]->noise_budget << std::endl;


								//std::cout << "DC " << tr <<"'s budget=" << DCs[tr]->aggregators[ai]->noise_budget << std::endl;
								if (DCs[tr]->aggregators[ai]->noise_budget > TIMES)
									communicate = true;

								if (privacy)//generate noise value 
								{
									float noise;
									noise = myapp->global_graph->laplace_generator(0, deltaf / DCs[tr]->aggregators[ai]->noise_budget);
									// LO: add noise
									int otherId = DCs[tr]->aggregators[ai]->DC_id;
									std::cout << "DC " << DCs[tr]->pri_rank << " send message to DC " << DCs[otherId]->pri_rank << std::endl;
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
									std::cout << "noise=" << noise << std::endl;
								}

								for (int nid = 0; nid < DCs[tr]->aggregators[ai]->v_list.size(); nid++) {
									int other_vertex = DCs[tr]->aggregators[ai]->v_list[nid];

									//msg队列记录的应该是PR(j)/L(j),PR(j)为网页j的PR值，L(j)是网页j的链出网页数
									VertexData* newmsg;
									if (myapp->mytype == pagerank)
									{
										newmsg = new PageRankVertexData();
										dynamic_cast<PageRankVertexData*>(newmsg)->rank = dynamic_cast<PageRankVertexData*>(DCs[tr]->aggregators[ai]->aggregated_data)->rank / (float)DCs[tr]->aggregators[ai]->v_list.size();
									}
									else if (myapp->mytype == sssp)
									{
										newmsg = new SSSPVertexData();
										dynamic_cast<SSSPVertexData*>(newmsg)->dist = (dynamic_cast<SSSPVertexData*>(DCs[tr]->aggregators[ai]->aggregated_data)->dist / (float)DCs[tr]->aggregators[ai]->v_list.size());

									}
									else if (myapp->mytype == subgraph)
									{
										newmsg = new SubgraphVertexData();
										//subgraph
									}

									(*myapp->global_graph->g)[other_vertex].messages.push_back(newmsg);
									(*myapp->global_graph->g)[other_vertex].messages_last.push_back(newmsg);  //messages_temp记录上次的message值

									/*
									wan_usage += abs(dynamic_cast<PageRankVertexData*>(newmsg)->rank);  //记录通信过程的wan usage
									//pregel模型一次通信cost只需计算发送方DC的上传费用以及接收方的下载价格
									totalcost += abs(dynamic_cast<PageRankVertexData*>(newmsg)->rank) * DCs[Threads[tr]->DC_loc_id]->upload_price / 1000.0; //per GB,记录通信过程的上传花费
									totalcost += abs(dynamic_cast<PageRankVertexData*>(newmsg)->rank)* DCs[Threads[tr]->DC_loc_id]->download_price / 1000.0; //per GB,记录通信过程的上传花费
									//耗时t的计算不同于GAS model,pregel model仅需发送消息一次（GAS两次）,时间仅由上传 or 下载最小的那个决定
									tr_send_step[tr] += abs(dynamic_cast<PageRankVertexData*>(newmsg)->rank) / DCs[Threads[tr]->DC_loc_id]->upload_band;
									tr_rcv_step[tr] += abs(dynamic_cast<PageRankVertexData*>(newmsg)->rank) / DCs[Threads[tr]->DC_loc_id]->download_band;
									*/

									wan_usage += myapp->msg_size;  //记录通信过程的wan usage
									//pregel模型一次通信cost只需计算发送方DC的上传费用以及接收方的下载价格
									totalcost += myapp->msg_size * DCs[Threads[tr]->DC_loc_id]->upload_price / 1000.0; //per GB,记录通信过程的上传花费
									totalcost += myapp->msg_size * DCs[Threads[tr]->DC_loc_id]->download_price / 1000.0; //per GB,记录通信过程的上传花费
									//耗时t的计算不同于GAS model,pregel model仅需发送消息一次（GAS两次）,时间仅由上传 or 下载最大的那个决定
									tr_send_step[tr] += myapp->msg_size / DCs[Threads[tr]->DC_loc_id]->upload_band;
									tr_rcv_step[tr] += myapp->msg_size / DCs[Threads[tr]->DC_loc_id]->download_band;
								}
							}
						}
						times_++;

					}
					else {  //如果不是指定的迭代次数，则保留着上次的message值
						times = 0;
						//std::cout << "这是间隔的iter" << std::endl;
						for (int ai = 0; ai < DCs[tr]->aggregators.size(); ai++) {
							std::vector<int> hasRecover;
							for (int nid = 0; nid < DCs[tr]->aggregators[ai]->v_list.size(); nid++) {
								int other_vertex = DCs[tr]->aggregators[ai]->v_list[nid];
								std::vector<int>::iterator it;
								it = find(hasRecover.begin(), hasRecover.end(), other_vertex);
								if (!(it != hasRecover.end()))  //other_vertex不在已恢复队列中，说明还没恢复该消息队列的值
								{
									hasRecover.push_back(other_vertex);
									for (int i = 0; i < (*myapp->global_graph->g)[other_vertex].messages_last.size(); i++) {

										(*myapp->global_graph->g)[other_vertex].messages.push_back((*myapp->global_graph->g)[other_vertex].messages_last[i]);
										//local通信不需要wan_usage以及cost,但是可能需要耗时t

										//wan_usage += dynamic_cast<PageRankVertexData*>((*myapp->global_graph->g)[other_vertex].messages_last[i])->rank;  //记录通信过程的wan usage
										//pregel模型一次通信cost只需计算发送方DC的上传费用以及接收方的下载价格
										//totalcost += dynamic_cast<PageRankVertexData*>((*myapp->global_graph->g)[other_vertex].messages_last[i])->rank* DCs[Threads[tr]->DC_loc_id]->upload_price / 1000.0; //per GB,记录通信过程的上传花费
										//totalcost += dynamic_cast<PageRankVertexData*>((*myapp->global_graph->g)[other_vertex].messages_last[i])->rank* DCs[Threads[tr]->DC_loc_id]->download_price / 1000.0; //per GB,记录通信过程的上传花费

									}
								}
							}
							hasRecover.clear();  //清空已恢复记录队列
						}

					}
				}
				//std::cout << "一轮结束" << std::endl;
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
				if (myapp->mytype == pagerank)
				{
					printf("rank值是：\n");
					for (int i = 0; i < rank_.size(); i++)
						std::cout << rank_[i] << std::endl;
				}
			
        } 
		else 
		{
            /* converge according to the tolerance */
            iter_counter ++;
			newCounter--;
        }
		
    }//end of pagerank
	
	else if(type == asynchronous){
        std::cout <<"asynchronous engine is not implemented. " << std::endl;
    }

    for(int tr=0; tr<num_threads; tr++) {
        printf("DC %d has %d aggregators\n",tr,(int)DCs[tr]->aggregators.size());
    }

	int isolated_vertices = 0;
	map<int,vector<float>> dc0;
	if (myapp->mytype == pagerank) {
		for (int i = 0; i < num_threads; i++) {
			// std::cout << "DC: " << i << std::endl;
			for (int v = 0; v < Threads[i]->l_dag->num_vertices; v++) {

				//if (dynamic_cast<PageRankVertexData*>((*Threads[i]->l_dag->g)[v].data)->rank < 0)
					//dynamic_cast<PageRankVertexData*>((*Threads[i]->l_dag->g)[v].data)->rank = -dynamic_cast<PageRankVertexData*>((*Threads[i]->l_dag->g)[v].data)->rank;
				if (dynamic_cast<PageRankVertexData*>((*Threads[i]->l_dag->g)[v].data)->rank - 0.15 < 0.01)
					isolated_vertices++;

				//output the pagerank value
				dc0[i].push_back(dynamic_cast<PageRankVertexData*>((*Threads[i]->l_dag->g)[v].data)->rank);
				printf("%f\n", dynamic_cast<PageRankVertexData*>((*Threads[i]->l_dag->g)[v].data)->rank);
			}
		}
	} 


	int size = atoi(sizeofgraph);
	map<int,int> sources_dc;  //The first int points to the DC number, and the second points to the source vertex ID
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

				//if (dynamic_cast<SSSPVertexData*>((*Threads[i]->l_dag->g)[v].data)->dist < 0)
					//dynamic_cast<SSSPVertexData*>((*Threads[i]->l_dag->g)[v].data)->dist = -dynamic_cast<SSSPVertexData*>((*Threads[i]->l_dag->g)[v].data)->dist;
				for (int j = 0; j < myapp->sources.size(); j++)
					if (vid == myapp->sources[j])
					{
						is_sources = true;
						sources_dc.insert(make_pair(i, vid));
						//sources_dc[i] = vid;
						break;
					}
				if (is_sources == false)
				{
					dc0[i].push_back(dynamic_cast<SSSPVertexData*>((*Threads[i]->l_dag->g)[v].data)->dist);
					//std::cout << "从0出发到"<< v << "点的最短距离为：" <<dynamic_cast<SSSPVertexData*>((*myapp->global_graph->g)[v].data)->dist << std::endl;
					std::cout << dynamic_cast<SSSPVertexData*>((*Threads[i]->l_dag->g)[v].data)->dist << std::endl;
				}
			} 
		}
	}
	 

	std::wcout << "抛弃顶点的distribution是：" << std::endl;
	for (int i = 0; i < abandoned_record.size(); i++) {
		std::wcout << abandoned_record[i] << std::endl;
	}


	/*
	//输出budget在每轮迭代，每个dc中的分布
	printf("budget_iter分布:\n"); 
	float last_temp = 0; 
	for (int i = 0; i < budget_iter.size(); i++)
	{
		if(last_temp != budget_iter[i])
			printf("%f\n", budget_iter[i]);
		last_temp = budget_iter[i];
	}

	printf("budget_agg分布:\n");
	last_temp = 0;
	for (int i = 0; i < budget_dc.size(); i++) {
		if (last_temp != budget_dc[i])
			printf("%f\n", budget_dc[i]);
		last_temp = budget_dc[i];
	}
	*/



	for (int tr = 0; tr < num_threads; tr++)
	{
		std::cout << "#################################################" << std::endl
			<< tr << "'s upload_band=" << DCs[Threads[tr]->DC_loc_id]->upload_band << std::endl
			<< "#################################################" << std::endl
			<< tr << "'s download_band=" << DCs[Threads[tr]->DC_loc_id]->download_band << std::endl;

	}


	//std::cout << "wan_usage=" << wan_usage << ", total_cost=" << totalcost << ", t_step=" << t_step << std::endl;
	std::cout << "wan_usage=" << wan_usage << " " << t_step << std::endl;
	std::cout << "msg_size=" << myapp->msg_size << std::endl;
	std::cout << "sum budget=" << sum_noise << std::endl;
	
	map<int, int>::iterator iter;
	iter = sources_dc.begin();
	while (iter != sources_dc.end()) {
		std::cout << "source " << iter->second << " in DC " << iter->first << std::endl;
		iter++;
	}

	/*output the realtive error*/
	/*The 'datafile' should store 
	**the PageRank of the vertices
	**in the DC0 in the original graph
	*/
	string datafile; 
	if(myapp->mytype == pagerank)
		datafile = "pagerank_";
	else if(myapp->mytype == sssp)
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
			std::cout << "DC " << tr <<"'s num_vertices = " << totalCount << std::endl;
			std::cout << "DC " << tr << "'s Average relative error = " << average_relative_error << std::endl;
			totalCount = 0;
			errorSum = 0;
			index_to_pr = 0;
		}
		double average_relative_error_global = errorSum_global / (float)totalCount_global;
		std::cout  << "Total num_vertices = " << totalCount_global << std::endl;
		std::cout << "Average relative error = " << average_relative_error_global << std::endl;
		std::cout << "Number of isolated vertices = " << isolated_vertices << std::endl;
		
	}
}
