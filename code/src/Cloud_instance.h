#ifndef ENTITIES_H
#define ENTITIES_H

#include "Distributed_graph.h"

#endif

extern int NUM_DC;
extern float OnDemandLag;

enum Integer_VM{
	not_ready,
	ready,
	scheduled,
	running,
	failed
};

//consider only one type of VM
class VM{
public:
	float price;
	int DC_name;
	
	MyVertex* tk; //task running on the VM
	Integer_VM status;
	float turn_on; //turn on time
	float life_time;
};

//DataCenter definition
class DataCenter {
public:
	int id;
	int size;
	int pri_rank;   //骆佩君，新增DC隐私等级
	float budget_weight;  //ribo,Define the privacy weight of DC
	float budget_weight_percentage;
	float budget_sum_weight;  //Record the total number of times messages need to be sent
	Amazon_EC2_regions location;
	
	std::vector<VM*> vms; //the running vms in the current datacenter
	float data_size; //the total input data size of all vertices in this dc
	// std::vector<int> contained_vertices; //global id of vertex
	/**
	* cloud setting uses net_band
	* the network bandwidth between the current DC and other DCs, ordered by DC's ID
	*/
	// std::vector<float> net_band; 	
	// std::vector<float> net_latency; //the network latency between the current DC and other DCs
	float vm_price; //hourly price of instances, only consider m1.large instances
	/** hpc setting uses upload/download bandwidth*/
	float upload_band;
	float download_band;
	float upload_price;//outbound price to other regions, per GB
	float download_price;

	float g_upload_data = 0.0f;
	float g_dnload_data = 0.0f;
	float a_dnload_data = 0.0f;
	float a_upload_data = 0.0f;

	//privacy
	std::vector<Aggregator*> aggregators;

	DataCenter(){ 
		// net_band = std::vector<float>(num_regions, 3.0); 
		// net_band[location_id] = 100;
		// net_latency = std::vector<float>(num_regions, 0.001);
		// net_latency[location_id] = 0;
		upload_band = 100;
		download_band = 500;
		data_size = 0;	
		upload_price = 0.02f;
		download_price = 0;
	}
	DataCenter(Amazon_EC2_regions l){
		location = l;
		upload_band = 100;
		download_band = 500;
		if (l == useast){
			vm_price = 0.175f;
			upload_price = 0.02f;
		}
		else if (l == uswest_northcalifornia){
			vm_price = 0.19f;
			upload_price = 0.02f;
		}
		else if (l == uswest_oregon){
			vm_price = 0.175f;
			upload_price = 0.02f;
		}
		//else if (l == ap_seoul){
		//	vm_price = 1000;//do not support m1.large type
		//	net_price = 0.08;
		//}
		else if (l == ap_singapore){
			vm_price = 0.233f;
			upload_price = 0.09f;
		}
		else if (l == ap_sydney){
			vm_price = 0.233f;
			upload_price = 0.14f;
		}
		else if (l == ap_tokyo){
			vm_price = 0.243f;
			upload_price = 0.09f;
		}
		//else if (l == eu_frankfurt){
		//	vm_price = 1000;//do not support m1.large type
		//	net_price = 0.02;
		//}
		else if (l == eu_ireland){
			vm_price = 0.19f;
			upload_price = 0.02f;
		}
		else if (l == saeast){
			vm_price = 0.233f;
			upload_price = 0.16f;
		}
		// net_band = std::vector<float>(num_regions, 3.0);
		// net_band[location_id] = 100;
		// net_latency = std::vector<float>(num_regions, 0.001);
		// net_latency[location_id] = 0;
	}
};

