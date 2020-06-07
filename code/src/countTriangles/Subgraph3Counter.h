#ifndef SUBGRAPH3COUNTER_H
#define SUBGRAPH3COUNTER_H

#include "Graph1.h"

namespace Graph1
{
	/**
	 * Class used to store the result of subgraph count for subgraphs of size 3
	 * It is returned by the count_subgraphs method of Subgraph3Counter
	 */
	class Subgraph3Count
	{
	public:
		Subgraph3Count() = default;
		virtual ~Subgraph3Count() = default;

		inline void setTriangles(long int triangles) { _triangles = triangles; };
		inline long int getTriangles() const { return _triangles; };
		inline void incrementTriangles() { ++_triangles; };
		inline void setWedges(long int wedges) { _wedges = wedges; };
		inline long int getWedges() const { return _wedges; };
		inline void incrementWedges() { ++_wedges; };
		inline void incrementWedges(long int increment) { _wedges += increment; };

		/**
		 * Add two count together: the result will be the sum of the two counts (for every subgraph count)
		 */
		Subgraph3Count operator+(const Subgraph3Count& other) const;
		/**
		 * Add this count to another: the values of this count will be the sum of the two counts (for every subgraph count)
		 */
		Subgraph3Count& operator+=(const Subgraph3Count& other);
		
		/**
		 * Nicely print a count
		 */
		friend std::ostream& operator<<(std::ostream& os, const Subgraph3Count& count);

	private:
		long int _triangles = 0;
		long int _wedges = 0;
	};
	
	// We allow the class to be used in a reduction clause of an openmp statement
	#pragma omp declare reduction (+: Subgraph3Count : omp_out += omp_in)

	/**
	 * Count subgraphs of size 3 in a graph.
	 * This implementation needs directed, bidirectional graphs (it would actually be easy to adapt it to undirected graphs, but it's not needed now).
	 * 
	 * Algorithm implemented with ideas from: http://theory.stanford.edu/~tim/s14/l/l1.pdf CS167: Reading in Algorithms Counting Triangles by Tim Roughgarden
	 */
	class Subgraph3Counter
	{
	public:
		Subgraph3Counter() = default;
		virtual ~Subgraph3Counter() = default;
		
		/**
		 * Execute the count
		 * @param g The graph on which the counting is performed.
		 * @return A Subgraph3Count object containing the result of the count.
		 */
		Subgraph3Count count_subgraphs(const Graph& g);
	};
}
#endif /* SUBGRAPH3COUNTER_H */
