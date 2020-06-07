#ifndef DISTRIBUTED_GRAPH_H
#define DISTRIBUTED_GRAPH_H

#include <vector>
#include <random>
#include <mutex>
#include <unordered_map>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include "Geography.h"

using namespace boost;



enum Graph_type
{
    test,
    synthetic,
    livejournal,
    twitter,
    roadca
};

enum vertex_status{
    activated,
    deactivated
};

class VertexData{
public:
    /** size of raw input data in MB
    * for example, the tweets of a twitter user etc.
    */
    float input_size;
	virtual double size_of() { return 0; }
    VertexData(){
        input_size = 10;
    }
    VertexData& operator=(const VertexData& other)  {
        this->input_size = other.input_size;
        return *this;
    }
};

class PageRankVertexData:public VertexData{
public:
    /*for pagerank application*/
    float rank;
	//bool isAbandoned;
    float last_change;
    PageRankVertexData(){
        input_size = 10;
        rank = 0.0; last_change = 0.0;
    }
	//modified,float to double
    double size_of(){
        return 0.000008;
    }
    PageRankVertexData& operator=(const PageRankVertexData& other)  {
        VertexData::operator= (other);
        this->rank = other.rank;
        this->last_change = other.last_change;
        return *this;
    }
};

class SSSPVertexData:public VertexData{
public:
    /*for sssp application*/
    float dist;
    bool changed;
	float last_change;
    std::vector<int> routed_nodes;
    SSSPVertexData(){
        input_size = 10;
        dist = 0.0; changed = false;
    }
	//modified,float to double

    double size_of(){
        return 0.000008;
    }
    SSSPVertexData& operator=(const SSSPVertexData& other)  {
//		printf("sssp assignment called\n");
        VertexData::operator= (other);
        this->dist = other.dist;
        this->changed = other.changed;
        this->routed_nodes = other.routed_nodes;
        return *this;
    }
};

class SubgraphVertexData:public VertexData{
public:
    /*for subgraph isomorphism*/
    int v_pattern;

    class Message{
    public:
        /** Possible matching vertices pairs
        * the 1st int: vertex id in the pattern
        * the 2nd int: vertex id in the data graph
        */
        std::vector<std::pair<int,int> > pairs;
        std::vector<int> forwarding_trace;
    };
    std::vector<Message> matches;
    std::vector<std::vector<int> > matched_subgraphs;

    SubgraphVertexData(){
        input_size = 10;
        v_pattern = 0;
    }
	//modified.float to double
    double size_of(){
        return sizeof(Message)*matches.size();
        //return 1.0;
    }
    SubgraphVertexData& operator=(const SubgraphVertexData& other)  {
//		printf("subgraph assignment called\n");
        VertexData::operator= (other);
        this->v_pattern = other.v_pattern;
        this->matches = other.matches;
        this->matched_subgraphs = other.matched_subgraphs;
        return *this;
    }
};

class Aggregator{
public:
    int DC_id; //send message to DC_id
    int vertex_id;
	int source_vertex_id;
    VertexData* aggregated_data;
    float noise_budget;
    std::vector<int> v_list; //for DC-level aggregation
    Aggregator(){
        DC_id = vertex_id = -1;
        noise_budget = 0;
    }
};

class MyVertex{

public:
    /*identification*/
    int vertex_id;
    // char* username;
    /*init location*/
    Amazon_EC2_regions location_id;
    VertexData* data;
    VertexData* accum;
    // for Pregel model
    std::vector<VertexData*> messages;
	std::vector<VertexData*> messages_last;
	//std::vector<VertexData*> messages_temp;
    //for privacy
    int distant_nbr;

    vertex_status status;
    vertex_status next_status;

    /*replication*/
    std::vector<int> replica_location;
    // int expected_replicas;
    bool is_master;
	bool is_urgently;
    int master_location;
    int local_vid;

    /*vertex performance estimation*/
    float perf_estimation();

    MyVertex(){
        status = deactivated;
        next_status = deactivated;
    }

    /* std::vector<int> get_rep_location(){
        mtx.lock();
        std::vector<int> out = replica_location;
        mtx.unlock();
        return out;
    }

    void update_rep_location(std::vector<int> rep){
        mtx.lock();
        replica_location = rep;
        mtx.unlock();
    } */

};

class EdgeData{
public:
    int edge_id;
    int start_node_id;
    int end_node_id;
    float weight;
    Amazon_EC2_regions assigned_location;
};

typedef adjacency_list<boost::setS, boost::vecS, bidirectionalS, MyVertex, property<edge_weight_t, float> > Graph; //
typedef graph_traits<Graph>::vertex_descriptor Vertex;
typedef graph_traits<Graph>::vertex_iterator vertex_iter;
typedef graph_traits<Graph>::edge_iterator edge_iter;
typedef graph_traits<Graph>::out_edge_iterator out_edge_iterator;
typedef graph_traits<Graph>::in_edge_iterator in_edge_iterator;
typedef graph_traits<Graph>::edge_descriptor Edge;
typedef graph_traits<Graph>::adjacency_iterator adja_iterator;

class distributed_graph{
public:

    Graph_type graph_type;
    Graph* g;
	Graph* original_g;
	Graph* temp_g;
    property_map<Graph, edge_weight_t>::type WeightMap = get(edge_weight, *g);
    // implemented for random access to edges, which is not supported in boost::graph
    std::unordered_map<int,Edge> random_access_edges;

    int num_vertices;
    int num_edges;

    // std::vector<LockV> mutexlocks;
    // std::vector<std::mutex> mutexlocks;

    /*indicating at initilization stage
    * for SubgraphIsomorphism app only
    */
    bool init_indicator = true;

    std::vector<int> get_in_nbrs(int v){
        in_edge_iterator in_i, in_end;
        Edge e; 
        std::vector<int> in_nbrs;
        boost::tie(in_i, in_end) = in_edges(v, *g);
        for (; in_i != in_end; ++in_i) {//if only one parent
            e = *in_i;
            Vertex src = source(e, *g);
            in_nbrs.push_back(src);
        }
        return in_nbrs;
    }

    std::vector<int> get_out_nbrs(int v){
        out_edge_iterator out_i, out_end;
        Edge e;
        std::vector<int> out_nbrs;
        boost::tie(out_i, out_end) = out_edges(v, *g);
        for (; out_i != out_end; ++out_i){
            e = *out_i;
            Vertex tgt = target(e, *g);
            out_nbrs.push_back(tgt);
        }
        return out_nbrs;
    }

    void pdf2cdf(std::vector<double>& pdf) {
        double Z = 0;
        for (size_t i = 0; i < pdf.size(); ++i) Z += pdf[i];
        for (size_t i = 0; i < pdf.size(); ++i)
            pdf[i] = pdf[i] / Z + ((i>0) ? pdf[i - 1] : 0);
    }

    double rnd_generator(double lower, double upper){
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(lower, upper);
        return dis(gen);
    }

    //Laplace distribution
    double laplace_generator(double mu, double sigma){
        double rnd=rnd_generator(0,1) - 0.5;
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
    /**
    * Generate a draw from a multinomial using a CDF.  This is
    * slightly more efficient since normalization is not required
    * and a binary search can be used.
    */
    size_t multinomial_cdf(const std::vector<double>& cdf) {
        double rnd = rnd_generator(0,1);
        return std::upper_bound(cdf.begin(), cdf.end(),rnd) - cdf.begin();
    }

    /**
    *   the line parser returns true if the line is parsed successfully and
    *	calls graph.add_vertex(...) or graph.add_edge(...)
    */
    bool line_parser(const std::string& filename,const std::string& linename, int ecount);


    /**
    * \brief Constructs a synthetic power law graph.
    *
    * This function constructs a synthetic out-degree power law of "nverts"
    * vertices with a particular alpha parameter.
    * In other words, the probability that a vertex has out-degree \f$d\f$,
    * is given by:
    *
    * \f[ P(d) \propto d^{-\alpha} \f]
    *
    * By default, the out-degree distribution of each vertex
    * will have power-law distribution, but the in-degrees will be nearly
    * uniform. This can be reversed by setting the second argument "in_degree"
    * to true.
    *
    * \param nverts Number of vertices to generate
    * \param in_degree If set to true, the graph will have power-law in-degree.
    *                  Defaults to false.
    * \param alpha The alpha parameter in the power law distribution. Defaults
    *              to 2.1
    * \param truncate Limits the maximum degree of any vertex. (thus generating
    *                 a truncated power-law distribution). Necessary
    *                 for large number of vertices (hundreds of millions)
    *                 since this function allocates a PDF vector of
    *                 "nverts" to sample from.
    */
    void load_synthetic_powerlaw(size_t nverts, bool in_degree = false,
                                 double alpha = 2.1, size_t truncate = (size_t)(-1));

    /**
    * file format: src tgt weight
    */
    void load_from_file(std::string, int, int);

    void load_from_file_1k_series(std::string ,int, int);

    /** load synthetic graphs from file
    * vertices: 10,000; 50,000; 100,000
    */
    void load_synthetic_file(const char* filename);
    /**
     * load synthetic graphs from file and use 1k-series to regenerate new graphs
     * @param filename the name of the Graph
     * @author Jerry Zhang
     */
    void load_synthetic_file_1k_series(const char* filename);

    /**
    * \brief LiveJournal social network graph
    * citation: Mikl��s Kurucz, Andr��s A. Bencz��r, Attila Pereszl��nyi. Large-Scale Principal Component Analysis on LiveJournal Friends Network. In proc Workshop on Social Network Mining and Analysis Held in conjunction with The 13th ACM SIGKDD International Conference on Knowledge Discovery and Data Mining (KDD 2008) August 24-27, 2008, Las Vegas, NV
    * source: https://dms.sztaki.hu/en/letoltes/livejournal-data
    * #nodes: 3,577,166
    * #edges: 44,913,072 directed edges
    * #countries: 244
    */
    void load_live_journal_graph();

    /**
    * \brief Twitter social graph
    * citation: J. McAuley and J. Leskovec. Learning to Discover Social Circles in Ego Networks. NIPS, 2012.
    * source: https://snap.stanford.edu/data/egonets-Twitter.html
    * #nodes: 81,306
    * #edges: 1,768,149
    */
    void load_twitter_graph();

    /**
    * \brief Facebook social network
    * citation: J. McAuley and J. Leskovec. Learning to Discover Social Circles in Ego Networks. NIPS, 2012.
    * source: https://snap.stanford.edu/data/egonets-Facebook.html
    * #nodes: 4,039
    * #edges: 88,234
    * undirected
    */
    void load_facebook_graph();

    /**
    * \brief Google Web Graph
    * citation: J. Leskovec, K. Lang, A. Dasgupta, M. Mahoney. Community Structure in Large Networks: Natural Cluster Sizes and the Absence of Large Well-Defined Clusters. Internet Mathematics 6(1) 29--123, 2009.
    * source: https://snap.stanford.edu/data/web-Google.html
    * #nodes: 875,713
    * #edges: 5,105,039
    */
    void load_googleweb_graph();

    /**
    * \brief California Road network
    * citation: J. Leskovec, K. Lang, A. Dasgupta, M. Mahoney. Community Structure in Large Networks: Natural Cluster Sizes and the Absence of Large Well-Defined Clusters. Internet Mathematics 6(1) 29--123, 2009.
    * source: https://snap.stanford.edu/data/roadNet-CA.html
    * #nodes: 1,965,206
    * #edges: 2,766,607
    */
    void load_roadnet_graph();

    /**
    * \brief Wiki-Vote graph
    * Directed graph. Wikipedia voting on promotion to administratorship (till January 2008).
    * Directed edge A->B means user A voted on B becoming Wikipedia administrator.
    * Nodes: 7115 Edges: 103689
    */
    void load_wiki_graph();

    /**
    * \brief Gnutella graph
    * # Nodes: 8846 Edges: 31839
    */
    void load_p2p_graph();

    /*a simple and controllable graph*/
    void test_graph();

    ~distributed_graph(){delete g; g=NULL; random_access_edges.clear();}
};



//*********************** For differential privacy purpose *************************//
//add laplace noise
double laplace_generator(double lower, double upper);

#endif
//neighboring graph with distance=1
