#include <random>
#include <iostream>
#include <fstream>
#include <algorithm>
#include<string>
using namespace std;

void pdf2cdf(std::vector<double>& pdf) {
	double Z = 0;
	for (size_t i = 0; i < pdf.size(); ++i) Z += pdf[i];
	for (size_t i = 0; i < pdf.size(); ++i)
		pdf[i] = pdf[i] / Z + ((i>0) ? pdf[i - 1] : 0);
}

float rnd_generator(float lower, float upper){
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<> dis(lower, upper);
	return dis(gen);	
}

/**
* Generate a draw from a multinomial using a CDF.  This is
* slightly more efficient since normalization is not required
* and a binary search can be used.
*/
size_t multinomial_cdf(const std::vector<double>& cdf) {
	double rnd = rnd_generator(0,1);
	return std::upper_bound(cdf.begin(), cdf.end(),rnd) - cdf.begin();
}

//generate powerlaw graph with edge weights
int main1(int argc, char* argv[]){

	float alpha = atof(argv[2]);
	size_t nverts = (size_t)atoi(argv[1]);
	size_t truncate = (size_t)(-1);
	bool in_degree = true;

	std::vector<double> prob(std::min(nverts, truncate), 0);
//	std::cout << "constructing pdf" << std::endl;
	for (size_t i = 0; i < prob.size(); ++i)
		prob[i] = std::pow(double(i + 1), -alpha);
//	std::cout << "constructing cdf" << std::endl;
	pdf2cdf(prob);
//	std::cout << "Building graph" << std::endl;	
	size_t target_index = 0;
	size_t addedvtx = 0;
	
	ofstream myfile;
	string filename = "powerlaw-"+to_string(nverts)+"-"+to_string(alpha);
	myfile.open(filename);	
	// A large prime number
	const size_t HASH_OFFSET = 2654435761;
	for (size_t source = 0; source < nverts; source++) {
		const size_t out_degree = multinomial_cdf(prob) + 1;
		//printf("# out edges is %d\n", out_degree);
		for (size_t i = 0; i < out_degree; ++i) {
			target_index = (target_index + HASH_OFFSET) % nverts;
			while (source == target_index) {
				target_index = (target_index + HASH_OFFSET) % nverts;
			}
			float weight = rnd_generator(0,1000);
			if (in_degree) 
				myfile << target_index << "\t" << source << "\t" << weight << endl;
			else 
				myfile << source << "\t" << target_index << "\t" << weight << endl;
		}
		++addedvtx;
	}
	myfile.close();
	return 0;
}
