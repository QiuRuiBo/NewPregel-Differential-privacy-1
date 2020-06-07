#include "Graph1.h"

namespace Graph1
{
	Vertex::Vertex(boost_vertex_t boost_vertex) :
		_vertex(boost_vertex) 
	{}
	
	std::ostream& operator<<(std::ostream& os, const Vertex& v)
	{
		return os << v._vertex;
	}
	
	Edge::Edge(boost_edge_t boost_edge):
		_edge(boost_edge)
	{}
	
	std::ostream& operator<<(std::ostream& os, const Edge& e)
	{
		return os << e._edge;
	}
	
	VertexIterator::VertexIterator(boost_vertex_iter_t boost_iterator):
		_iter(boost_iterator)
	{}
	
	OutEdgeIterator::OutEdgeIterator(boost_out_edge_iter_t boost_iterator):
		_iter(boost_iterator)
	{}
	
	InEdgeIterator::InEdgeIterator(boost_in_edge_iter_t boost_iterator):
		_iter(boost_iterator)
	{}
	
	Graph::Graph(int nb_vertices):
		_g(nb_vertices)
	{}
	
	int Graph::vertex_count() const
	{
		return boost::num_vertices(_g);
	}
	
	int Graph::edge_count() const
	{
		return boost::num_edges(_g);
	}
	
	bool Graph::is_directed()
	{
		#if GRAPH_TYPE == UNDIRECTED_GRAPH
		return false;
		#else
		return true;
		#endif
	}
	
	bool Graph::is_bidirectional()
	{
		#if GRAPH_TYPE == DIRECTED_GRAPH // Undirected graphs are considered bidirectional, so only directed is not
		return false;
		#else
		return true;
		#endif
	}
	
	void Graph::clear()
	{
		_g.clear();
	}
	
	Vertex Graph::add_vertex()
	{
		return Vertex(boost::add_vertex(_g));
	}
	
	void Graph::add_vertices(int nb_vertices)
	{
		for(int i = 0; i < nb_vertices; ++i)
			add_vertex();
	}
	
	void Graph::add_vertices(int nb_vertices, std::vector<Vertex>& vertex_ids)
	{
		for(int i = 0; i < nb_vertices; ++i)
			vertex_ids.push_back(add_vertex());
	}
	
	Vertex Graph::get_vertex(int id) const
	{
		return Vertex(boost::vertex(id, _g));
	}
	
	std::pair<VertexIterator, VertexIterator> Graph::get_vertices() const
	{
		std::pair<boost_vertex_iter_t, boost_vertex_iter_t> boost_iters = boost::vertices(_g);
		return std::make_pair(VertexIterator(boost_iters.first), VertexIterator(boost_iters.second));
	}
	
	void Graph::del_vertex(int id)
	{
		boost::remove_vertex(id, _g);
	}
	
	void Graph::del_vertices(const std::vector<Vertex>& vertices)
	{
		for (const Vertex& v : vertices)
		{
			del_vertex(v);
		}
	}
	
	Edge Graph::add_edge(int source, int target)
	{
		return Edge(boost::add_edge(source, target, _g).first);
	}
	
	Edge Graph::get_edge(int source, int target) const
	{
		std::pair<boost_edge_t, bool> boost_result = edge(source, target, _g);
		if (!boost_result.second)
			throw EdgeNotFoundException();
		return Edge(boost_result.first);
	}
	
	void Graph::del_edge(int source, int target)
	{
		boost::remove_edge(source, target, _g);
	}
	
	bool Graph::has_edge(int source, int target) const
	{
		return boost::edge(source, target, _g).second;
	}
	
	std::pair<OutEdgeIterator, OutEdgeIterator> Graph::get_out_edges(int id) const
	{
		std::pair<boost_out_edge_iter_t, boost_out_edge_iter_t> boost_iters = boost::out_edges(id, _g);
		return std::make_pair(OutEdgeIterator(boost_iters.first), OutEdgeIterator(boost_iters.second));
	}
	
	std::pair<InEdgeIterator, InEdgeIterator> Graph::get_in_edges(int id) const
	{
		#if GRAPH_TYPE == DIRECTED_GRAPH
			throw UnsupportedOperationException("Unsupported operation for directed graphs. Try compiling the Graph class with GRAPH_TYPE = BIDIRECTIONAL_GRAPH");
		#else
			std::pair<boost_in_edge_iter_t, boost_in_edge_iter_t> boost_iters = boost::in_edges(id, _g);
			return std::make_pair(InEdgeIterator(boost_iters.first), InEdgeIterator(boost_iters.second));
		#endif
	}
	
	int Graph::get_out_degree(int vertex_id) const
	{
		return boost::out_degree(vertex_id, _g);
	}
	
	int Graph::get_in_degree(int vertex_id) const
	{
		#if GRAPH_TYPE == DIRECTED_GRAPH
			throw UnsupportedOperationException("Unsupported operation for directed graphs. Try compiling the Graph class with GRAPH_TYPE = BIDIRECTIONAL_GRAPH");
		#else
			return boost::in_degree(vertex_id, _g);
		#endif
	}
	
	void Graph::clear_vertex(int vertex_id)
	{
		boost::clear_vertex(vertex_id, _g);
	}
	
	void Graph::clear_out_edges(int vertex_id)
	{
		#if GRAPH_TYPE == UNDIRECTED_GRAPH
			throw UnsupportedOperationException("Unsupported operation for undirected graphs. Use clear_vertex instead");
		#else
		boost::clear_out_edges(vertex_id, _g);
		#endif
	}
	
	void Graph::clear_in_edges(int vertex_id)
	{
		#if GRAPH_TYPE != BIDIRECTIONAL_GRAPH
			throw UnsupportedOperationException("Operation supported only on bidirectional graphs. Compile the Graph class with GRAPH_TYPE = BIDIRECTIONAL_GRAPH, or use clear_vertex instead");
		#else
			boost::clear_in_edges(vertex_id, _g);
		#endif
	}
	/*
	std::shared_ptr<VertexData>& Graph::vertex_data(int vertex_id)
	{
		return _g[vertex_id];
	}
	
	std::shared_ptr<EdgeData>& Graph::edge_data(const Edge& e)
	{
		return _g[e.get_boost_edge()];
	}
	*/
}