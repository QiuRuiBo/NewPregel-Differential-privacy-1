#ifndef VERTEXDATA_H
#define VERTEXDATA_H

namespace Graph1
{
	/**
	 * Abstract class used to store data in a graph.
	 * To store custom data, derive from this class.
	 */
	class VertexData
	{
	public:
		VertexData() = default;
		virtual ~VertexData() = default;
	};
}

#endif /* VERTEXDATA_H */

