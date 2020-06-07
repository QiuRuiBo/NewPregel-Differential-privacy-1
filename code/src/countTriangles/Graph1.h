/**
 * Wrapper classes around boost graph library
 */
#ifndef GRAPH1_H
#define GRAPH1_H

#include <boost/graph/adjacency_list.hpp>
#include <boost/operators.hpp>
#include <iostream>
#include <vector>
#include <memory>
#include "GraphException.h"
#include "VertexData.h"
#include "EdgeData.h"

//#include <boost/graph/graph_traits.hpp>
#include "../Distributed_graph.h"

// Boost graph type to use. Note: Some functions may not be available for all types (e.g. get_in_edges)
#define UNDIRECTED_GRAPH 1
#define DIRECTED_GRAPH 2
#define BIDIRECTIONAL_GRAPH 3 // Directed graph, but a vertex can also access its in-edges (requires more space)
#ifndef GRAPH_TYPE
	#define GRAPH_TYPE BIDIRECTIONAL_GRAPH // Choose here the type of graph to be used
#endif
#if GRAPH_TYPE == UNDIRECTED_GRAPH
	#define BOOST_TYPE boost::undirectedS
#elif GRAPH_TYPE == DIRECTED_GRAPH
	#define BOOST_TYPE boost::directedS
#elif GRAPH_TYPE == BIDIRECTIONAL_GRAPH
	#define BOOST_TYPE boost::bidirectionalS
#endif

namespace Graph1
{
	// The type of graph that will be used
	//typedef adjacency_list<boost::setS, boost::vecS, bidirectionalS, MyVertex, property<edge_weight_t, float> > Graph; //
	typedef boost::adjacency_list<boost::setS, boost::vecS, bidirectionalS, MyVertex, property<edge_weight_t, float> > boost_graph_t;
	typedef boost::graph_traits<boost_graph_t>::vertex_descriptor boost_vertex_t;
	typedef boost::graph_traits<boost_graph_t>::edge_descriptor boost_edge_t;
	typedef boost::graph_traits<boost_graph_t>::vertex_iterator boost_vertex_iter_t;
	typedef boost::graph_traits<boost_graph_t>::out_edge_iterator boost_out_edge_iter_t;
	typedef boost::graph_traits<boost_graph_t>::in_edge_iterator boost_in_edge_iter_t;

	/**
	 * Wrapper class around boost vertex descriptor
	 */

	class Vertex : boost::totally_ordered<Vertex>
	{
	public:
		Vertex() = default;
		explicit Vertex(boost_vertex_t boost_vertex);
		virtual ~Vertex() = default;
		
		/**
		 * @return The vertex id
		 */
		inline int id() const { return _vertex; }
		
		inline bool operator==(const Vertex& v) const { return this->id() == v.id(); };
		inline bool operator<(const Vertex& v) const { return this->id() < v.id(); };
		friend std::ostream& operator<<(std::ostream& os, const Vertex& v);
		
	private:
		boost_vertex_t _vertex;
	};
	
	/**
	 * Wrapper class around boost edge descriptor
	 */
	class Edge
	{
	public:
		Edge() = default;
		explicit Edge(boost_edge_t boost_edge);
		virtual ~Edge() = default;
		
		/**
		 * @return The source vertex of this edge
		 */
		inline Vertex source() const { return Vertex(_edge.m_source); };
		/**
		 * @return The target vertex of this edge
		 */
		inline Vertex target() const { return Vertex(_edge.m_target); };
		
		/**
		 * @return The underlying boost edge descriptor
		 */
		inline const boost_edge_t& get_boost_edge() const { return _edge; };
		
		friend std::ostream& operator<<(std::ostream& os, const Edge& e);
		
	private:
		boost_edge_t _edge;
	};
	
	/**
	 * Wrapper around boost vertex iterator
	 */
	class VertexIterator
	{
	public:
		VertexIterator() = default;
		explicit VertexIterator(boost_vertex_iter_t boost_iterator);
		virtual ~VertexIterator() = default;
		
		/**
		 * Compare with another iterator
		 */
		inline bool operator==(const VertexIterator& v) const { return this->_iter == v._iter; };
		/**
		 * Compare with another iterator
		 */
		inline bool operator!=(const VertexIterator& v) const {return !operator==(v); };
		/**
		 * Move to the next item
		 */		
		inline VertexIterator& operator++() { ++_iter; return *this; };
		/**
		 * @return The vertex pointed by the iterator
		 */
		inline Vertex operator*() const { return Vertex(*_iter); };
	
	private:
		boost_vertex_iter_t _iter;
	};
	
	/**
	 * Wrapper around boost out_edge iterator
	 */
	class OutEdgeIterator
	{
	public:
		OutEdgeIterator() = default;
		explicit OutEdgeIterator(boost_out_edge_iter_t boost_iterator);
		virtual ~OutEdgeIterator() = default;
		
		/**
		 * Compare with another iterator
		 */
		inline bool operator==(const OutEdgeIterator& it) const { return this->_iter == it._iter; };
		/**
		 * Compare with another iterator
		 */
		inline bool operator!=(const OutEdgeIterator& it) const {return !operator==(it); };
		/**
		 * Move to the next item
		 */
		inline OutEdgeIterator& operator++() { ++_iter; return *this; };
		/**
		 * @return The edge pointed by the iterator
		 */
		inline Edge operator*() const { return Edge(*_iter); };
		
	private:
		boost_out_edge_iter_t _iter;
	};
	
	/**
	 * Wrapper around boost in_edge iterator.
	 * Not available for simple directed graph 
	 * (graph needs to be bidirectional for a vertex to be able to access in-edges)
	 */
	class InEdgeIterator
	{
	public:
		InEdgeIterator() = default;
		explicit InEdgeIterator(boost_in_edge_iter_t boost_iterator);
		virtual ~InEdgeIterator() = default;
		
		/**
		 * Compare with another iterator
		 */
		inline bool operator==(const InEdgeIterator& it) const { return this->_iter == it._iter; };
		/**
		 * Compare with another iterator
		 */
		inline bool operator!=(const InEdgeIterator& it) const {return !operator==(it); };
		/**
		 * Move to the next item
		 */
		inline InEdgeIterator& operator++() { ++_iter; return *this; };
		/**
		 * @return The edge pointed by the iterator
		 */
		inline Edge operator*() const { return Edge(*_iter); };
		
	private:
		boost_in_edge_iter_t _iter;
	};
	
	/**
	 * Wrapper class around boost adjacency_list graph representation
	 */
	class Graph
	{
	public:
		/**
		 * Create a new graph
		 * @param nb_vertices The number of vertices that the graph will have initially
		 */
		Graph(int nb_vertices = 0);
		virtual ~Graph() = default;
		
		/**
		 * @return The number of vertices in the graph
		 */
		int vertex_count() const;
		/**
		 * @return The number of edges in the graph
		 */
		int edge_count() const;
		/**
		 * @return Whether the graphs we are using are directed or not
		 */
		static bool is_directed();
		/**
		 * @return Whether the graphs we are using are bidirectional (i.e. from a vertex we can access both in-edges and out-edges)
		 * Note: undirected graphs are bidirectional by definition
		 */
		static bool is_bidirectional();
		
		/**
		 * Remove all edges and vertices from the graph
		 */
		void clear();
		
		/**
		 * Add a new vertex to the graph
		 * @return The descriptor of the new vertex
		 */
		Vertex add_vertex();
		/**
		 * Add more than one vertex to the graph
		 * @param nb_vertices The number of vertices to add
		 */
		void add_vertices(int nb_vertices);
		/**
		 * Add more than one vertex to the graph and store the descriptors
		 * @param nb_vertices The number of vertices to add
		 * @param vertex_ids Vertex descriptor of the new vertices will be pushed into this vector.
		 * The vector is never cleared.
		 */
		void add_vertices(int nb_vertices, std::vector<Vertex>& vertex_ids);
		/**
		 * Get a vertex
		 * @param id The id of the vertex
		 */
		Vertex get_vertex(int id) const;
		/**
		 * Get all vertices in the graph
		 * @return Start and end iterators on the set of vertices
		 */
		std::pair<VertexIterator, VertexIterator> get_vertices() const;
		
		/**
		 * Delete a vertex from the graph
		 * @param id The id of the vertex to delete
		 */
		void del_vertex(int id);
		/**
		 * Delete a vertex from the graph
		 * @param v The vertex descriptor
		 */
		inline void del_vertex(const Vertex& v) { del_vertex(v.id()); };
		/**
		 * Delete a list of vertices from the graph
		 * @param vertices The vertices to delete
		 */
		void del_vertices(const std::vector<Vertex>& vertices);
		/**
		 * Add an edge to the graph
		 * @param source The source vertex of the edge
		 * @param target The target vertex of the edge
		 * @return The descriptor of the new edge
		 */
		Edge add_edge(int source, int target);
		/**
		 * Add an edge to the graph
		 * @param source The source vertex of the edge
		 * @param target The target vertex of the edge
		 * @return The descriptor of the new edge
		 */
		inline Edge add_edge(const Vertex& source, const Vertex& target) { return add_edge(source.id(), target.id()); };
		/**
		 * Get an edge from the graph
		 * @param source The source vertex of the edge
		 * @param target The target vertex of the edge
		 * @return The descriptor of the edge
		 * @throw EdgeNotFoundException if the edge doesn't exist
		 */
		Edge get_edge(int source, int target) const;
		/**
		 * Get an edge from the graph
		 * @param source The source vertex of the edge
		 * @param target The target vertex of the edge
		 * @return The descriptor of the edge
		 * @throw EdgeNotFoundException if the edge doesn't exist
		 */
		inline Edge get_edge(const Vertex& source, const Vertex& target) const { return get_edge(source.id(), target.id()); };
		/**
		 * Delete an edge from the graph
		 * @param source The source vertex of the edge
		 * @param target The target vertex of the edge
		 */
		void del_edge(int source, int target);
		/**
		 * Delete an edge from the graph
		 * @param source The source vertex of the edge
		 * @param target The target vertex of the edge
		 */
		inline void del_edge(const Vertex& source, const Vertex& target) { del_edge(source.id(), target.id()); };
		/**
		 * Delete an edge from the graph
		 * @param e The edge to delete
		 */
		inline void del_edge(const Edge& e) { del_edge(e.source(), e.target()); };
		/**
		 * Checks whether an edge exist
		 */
		bool has_edge(int source, int target) const;
		/**
		 * Checks whether an edge exist
		 */
		inline bool has_edge(const Vertex& source, const Vertex& target) const { return has_edge(source.id(), target.id()); };
		/**
		 * Checks whether an edge exist
		 */
		inline bool has_edge(const Edge& e) const { return has_edge(e.source().id(), e.target().id()); };
		
		/**
		 * Get the edges starting from a vertex
		 * @param id The vertex id
		 * @return Start and end iterators on the set of edges
		 */
		std::pair<OutEdgeIterator, OutEdgeIterator> get_out_edges(int id) const;
		/**
		 * Get the edges starting from a vertex
		 */
		inline std::pair<OutEdgeIterator, OutEdgeIterator> get_out_edges(const Vertex& v) const { return get_out_edges(v.id()); };
		
		/**
		 * Get the edges ending on a vertex.
		 * Note: the vertex will be the 'target' attribute of the edges
		 * @param id The vertex id.
		 * @return Start and end iterators on the set of edges
		 * 
		 * If the class is compiled with GRAPH_TYPE == DIRECTED_GRAPH, this throws a Graph::UnsupportedOperationException
		 */
		std::pair<InEdgeIterator, InEdgeIterator> get_in_edges(int id) const;
		/**
		 * Get the edges ending on a vertex.
		 * Note: the vertex will be the 'target' attribute of the edges
		 * 
		 * If the class is compiled with GRAPH_TYPE == DIRECTED_GRAPH, this throws a Graph::UnsupportedOperationException
		 */
		inline std::pair<InEdgeIterator, InEdgeIterator> get_in_edges(const Vertex& v) const { return get_in_edges(v.id()); };

		/**
		 * Return the number of edges leaving a vertex
		 * @param vertex_id The vertex id
		 */
		int get_out_degree(int vertex_id) const;
		/**
		 * Return the number of edges leaving a vertex
		 * @param v The vertex
		 */
		inline int get_out_degree(const Vertex& v) const { return get_out_degree(v.id()); };
		
		/**
		 * Return the number of edges entering a vertex
		 * @param vertex_id The vertex id
		 * 
		 * If the class is compiled with GRAPH_TYPE == DIRECTED_GRAPH, this throws a Graph::UnsupportedOperationException
		 */
		int get_in_degree(int vertex_id) const;
		/**
		 * Return the number of edges entering a vertex
		 * @param v The vertex
		 * 
		 * If the class is compiled with GRAPH_TYPE == DIRECTED_GRAPH, this throws a Graph::UnsupportedOperationException
		 */
		inline int get_in_degree(const Vertex& v) const { return get_in_degree(v.id()); };
		
		/**
		 * Remove all edges to and from a vertex. The vertex itself is not removed
		 */
		void clear_vertex(int vertex_id);
		/**
		 * Remove all edges to and from a vertex. The vertex itself is not removed
		 */
		inline void clear_vertex(const Vertex& v) { clear_vertex(v.id()); };
		/**
		 * Removes all out-edges of a vertex. The vertex itself is not removed.
		 * 
		 * Not defined for undirected graphs, use clear_vertex instead.
		 */
		void clear_out_edges(int vertex_id);
		/**
		 * Removes all out-edges of a vertex. The vertex itself is not removed.
		 * 
		 * Not defined for undirected graphs, use clear_vertex instead.
		 */
		inline void clear_out_edges(const Vertex& v) { clear_out_edges(v.id()); };
		/**
		 * Removes all in-edges of a vertex. The vertex itself is not removed.
		 * 
		 * Only defined for bidirectional graphs, use clear_vertex or clear_out_edges for other types of graph
		 */
		void clear_in_edges(int vertex_id);
		/**
		 * Removes all in-edges of a vertex. The vertex itself is not removed.
		 * 
		 * Only defined for bidirectional graphs, use clear_vertex or clear_out_edges for other types of graph
		 */
		inline void clear_in_edges(const Vertex& v) { clear_in_edges(v.id()); };

		/**
		 * Get the VertexData associated with a vertex.
		 * Smart pointers are used to ensure deallocation of data.
		 * @param vertex_id The id of the vertex for which to get the data
		 * @return A reference to the smart pointer used to manage the data.
		 * 
		 * Since the return type is a reference to a smart pointer, it is possible
		 * to do vertx_data(id).reset(new CustomDataDerivedFromVertexData)
		 */
		std::shared_ptr<VertexData>& vertex_data(int vertex_id);
		/**
		 * Alias to vertex_data(int vertex_id)
		 */
		inline std::shared_ptr<VertexData>& vertex_data(const Vertex& v) { return vertex_data(v.id()); };
		
		/**
		 * Get the EdgeData associated with an edge.
		 * Smart pointers are used to ensure deallocation of data.
		 * @param e The edge for which to get the data
		 * @return A reference to the smart pointer used to manage the data.
		 * 
		 * Since the return type is a reference to a smart pointer, it is possible
		 * to do edge_data(e).reset(new CustomDataDerivedFromEdgeData)
		 */
		std::shared_ptr<EdgeData>& edge_data(const Edge& e);
		
		/**
		 * Get a vertex from the graph
		 * @param id The id of the vertex
		 */
		inline Vertex operator[](int id) const { return get_vertex(id); };
		
		
	public:
		boost_graph_t _g;
	};
}

#endif /* GRAPH1_H */
