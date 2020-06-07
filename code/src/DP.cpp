#include "Distributed_graph.h"
using namespace std;


//Laplace distribution
double laplace_generator(double mu, double sigma){
	distributed_graph* g = new distributed_graph();
	double rnd=g->rnd_generator(0,1) - 0.5;
	double b = sigma / sqrt(2.0);
	float sign = -1;
	if(rnd == 0)
		sign = 0;
	else if(rnd > 0)
		sign = 1;
	rnd = mu - b * sign * log(1- 2* rnd*sign);
	//boost::laplace_distribution<> l( _location, _scale );
	//laplace l;	

	return rnd;

}


