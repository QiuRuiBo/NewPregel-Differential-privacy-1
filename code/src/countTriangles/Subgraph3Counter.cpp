#include "Subgraph3Counter.h"

#include <unordered_set>

namespace Graph1
{
	
	Subgraph3Count Subgraph3Count::operator+(const Subgraph3Count& other) const
	{
		Subgraph3Count result;
		result._triangles = this->_triangles + other._triangles;
		result._wedges = this->_wedges + other._wedges;
		return result;
	}

	Subgraph3Count& Subgraph3Count::operator+=(const Subgraph3Count& other)
	{
		this->_triangles = this->_triangles + other._triangles;
		this->_wedges = this->_wedges + other._wedges;
		return *this;
	}
	
	std::ostream& operator<<(std::ostream& os, const Subgraph3Count& count)
	{
		os << "+- Subgraph count" << std::endl;
		os << "| Triangles: " << count._triangles << std::endl;
		os << "| Wedges: " << count._wedges << std::endl;
		os << "+-" << std::endl;
		return os;
	}
	
	Subgraph3Count Subgraph3Counter::count_subgraphs(const Graph& g)
	{
		Subgraph3Count count;
		#pragma omp parallel for schedule(dynamic) reduction(+:count)
		for (int cur_node = 0; cur_node < g.vertex_count(); ++cur_node)
		{
			int cur_node_degree = g.get_in_degree(cur_node) + g.get_out_degree(cur_node);
			std::unordered_set<int> neighbors;
			
			// Find neighbors of current node
			{
				OutEdgeIterator e_it, e_end;
				for(std::tie(e_it, e_end) = g.get_out_edges(cur_node); e_it != e_end; ++e_it)
					neighbors.insert((*e_it).target().id());
			}
			{
				InEdgeIterator e_it, e_end;
				for(std::tie(e_it, e_end) = g.get_in_edges(cur_node); e_it != e_end; ++e_it)
					neighbors.insert((*e_it).source().id());
			}
			
			// Wedges are equivalent to all possible pairings of two neighbors
			count.incrementWedges(neighbors.size() * (neighbors.size() -1) / 2); // Binomial coefficient for k = 2
			
			// Find triplets of nodes
			// For optimization we delegate the subgraph count to the node in the triplet with the lowest degree (or lowest id for same-degree nodes)
			for (int u : neighbors)
			{
				int u_degree = g.get_in_degree(u) + g.get_out_degree(u);
				if (cur_node_degree > u_degree)
					continue;
				if ((cur_node_degree == u_degree) && (u < cur_node))
					continue;
				
				for (int v : neighbors)
				{
					if (v <= u) // equal: we avoid v == u. Lower than: we optimize by only counting once the triplet (cur_node,u,v) and not counting (cur_node,v,u)
						continue;
					int v_degree = g.get_in_degree(v) + g.get_out_degree(v);
					if (cur_node_degree > v_degree)
						continue;
					if ((cur_node_degree == v_degree) && (v < cur_node))
					continue;
					
					if (g.has_edge(u, v) || g.has_edge(v, u))
						count.incrementTriangles();
				}
			}
		}
		return count;
	}
}