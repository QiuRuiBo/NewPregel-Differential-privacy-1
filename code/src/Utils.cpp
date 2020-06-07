//
// Created by Jerry Zhang on 2019-06-25.
//

#include "Utils.h"
#include <algorithm>
#include <string>
int getNum(std::vector<int>& v);
std::vector<int> utils::generateRandomPermutation(int size) {
    std::vector<int> v(size);
    std::vector<int> result;
    // Fill the vector with the values
    // 1, 2, 3, ..., n
    for (int i = 0; i < size; i++)
        v[i] = i + 1;

    // While vector has elements
    // get a random number from the vector and print it
    while (!v.empty()) {
        result.push_back(getNum(v));
    }
    return result;
}
// Function to return the next random number
int getNum(std::vector<int>& v)
{
    // Size of the vector
    int n = v.size();
    // Make sure the number is within
    // the index range
    int index = rand() % n;

    // Get random number from the vector
    int num = v[index];

    // Remove the number from the vector
    std::swap(v[index], v[n - 1]);
    v.pop_back();

    // Return the removed number
    return num-1;
}
int utils::generateRandomNumber(int max) {
	//std::cout << "The max is:"<< max << std::endl;
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 eng(rd()); // seed the generator
    std::uniform_int_distribution<> distr(0, max); // define the range
    return distr(eng);
}


//Modified by Ribo:change here to add pi distribution(divided into an in result(matrix) and an out result),
//which point to the probability of both vertices' connected. 
int* utils::getInPIdistribution(Graph* g) {
	vertex_iter v_it, v_end;
	std::tie(v_it, v_end) = vertices(*g);
	int num = num_vertices(*g) + 1;
	int* result = new int[num];
	for (int i = 0; i < num; i++) {
		result[i] = 0;
	}
	int indicate=0;
	for (; v_it != v_end; v_it++) {
		result[indicate++]=in_degree(*v_it, *g);             //in_degree() function will return in_edge_list().size();
		//std::cout << "In:" <<result[--indicate] <<std::endl;
	}
	return result;
}


//generate a random number from 0-(max vertex).
int utils::generateRandomNumberFrom0TOmax(int max) {
	//std::cout << "The max is:" << max << std::endl;
	//srand((unsigned)time(NULL));
	int i= (int)(rand() / (double)RAND_MAX * max);
	while(i>=max)
		i = (int)(rand() / (double)RAND_MAX * max);
	return i; //the value of RAND_MAX is 32767
}



int* utils::getOutPIdistribution(Graph* g) {
	vertex_iter v_it, v_end;
	std::tie(v_it, v_end) = vertices(*g);
	int num = num_vertices(*g) + 1;
	int* result = new int[num];
	for (int i = 0; i < num; i++) {
		result[i] = 0;
	}
	int indicate = 0;
	for (; v_it != v_end; v_it++) {
		result[indicate++]= out_degree(*v_it, *g);
		//std::cout << "Out:" << result[--indicate] << std::endl;
	}
	return result;
}


int* utils::getIn1kSeries(Graph *g) {
    vertex_iter v_it,v_end;
    std::tie(v_it,v_end) = vertices(*g);
    int num = num_vertices(*g) + 1;
    int *result = new int[num];

	//std::cout << "Result martix have " << num << " numbers." << std::endl;

    for (int i = 0;i < num;i++) {
        result[i] = 0;
    }
    for(;v_it != v_end;v_it++) {
        result[in_degree(*v_it,*g)]++;             //in_degree() function will return in_edge_list().size();
    }
    return result;
}
int* utils::getOut1kSeries(Graph *g) {
    vertex_iter v_it,v_end;
    std::tie(v_it,v_end) = vertices(*g);
    int num = num_vertices(*g) + 1;
    int *result = new int[num];
    for (int i = 0;i < num;i++) {
        result[i] = 0;
    }
    for(;v_it != v_end;v_it++) {
        result[out_degree(*v_it,*g)]++;
    }
    return result;
}

//ribo,Find the position of the corresponding elements in the array that satisfy a certain relationship
int utils::findLocation(int a[], int size, int which, Graph& g)
{
	int location = 0, index = -1;

	while (index != which)
	{
		if (g[location].location_id == 0)
			location++;
		else
		{
			location++;
			index++;
		}
		if (location > size)
		{
			printf("search NO %d failed.", which);
			printf("Out of index.\n");
			return -1;
		}
	}
	return location - 1;
}

vector<int> utils::adjacent_vertices(int u,Graph& g) {  //Looks for adjacent vertices in a directed graph
	vector<int> result;
	edge_iter ei, ei_end;
	std::tie(ei, ei_end) = boost::edges(g);
	for (auto next=ei; next!=ei_end; next++)
	{
		if (boost::source(*next, g) == u)
		{
			result.push_back(target(*next, g));
		}

		if (boost::target(*next, g) == u)
		{
			result.push_back(source(*next, g));
		}
	}
	return result;
}

vector<Edge> utils::selfloop_edges(Graph& g) {
	vector<Edge> self_loop_edges;
	edge_iter ei, ei_end;
	std::tie(ei, ei_end) = boost::edges(g);
	for (edge_iter next = ei; next != ei_end; next++) {
		if (source(*next, g) == target(*next, g)) {
			self_loop_edges.push_back(*next);
		}
	}
	return self_loop_edges;
}

float utils::feibonaqi_n(float n, float n_front, int k_num) {
	/**
	 * :n:第n项斐波那契数列值
	 * :n_front:第n-1项斐波那契数列值
	 * :k_num:当不传递n以及n_front值时，该项要求，表示返回第k_num项的值
	 * :return:第n+1项斐波那契数列值
	*/
	if (n && n_front)
		return n + n_front;

	else {
		float k1 = 1, k2 = 1;
		for (int i = 0; i < k_num - 2; i++) {
			float temp = k2;
			k2 = feibonaqi_n(k1, k2, 0);
			k1 = temp;
		}
		return k2;
	}
}

float utils::feibonaqi_sum(int n) {
	/**
 * :n:满足一定关系的n项斐波那契数列和
 * :result[]:用于保存结果的数组
*/
	float k1 = 1, k2 = 1, start = 1, sum = 2;
	for (int i = 0; i < n - 2; i++) {  //先计算从第一项开始的前n项和
		float temp = k2;
		k2 = feibonaqi_n(k1, k2, 0);
		sum = k2 + sum;
		k1 = temp;
	}
	return sum;
}