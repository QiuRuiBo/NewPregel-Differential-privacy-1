#ifndef EDGEDATA_H
#define EDGEDATA_H

namespace Graph1
{
	/**
	 * Abstract class used to store data in a graph.
	 * To store custom data, derive from this class.
	 */
	class EdgeData
	{
	public:
		EdgeData() = default;
		virtual ~EdgeData() = default;
	};
}

#endif /* EDGEDATA_H */

