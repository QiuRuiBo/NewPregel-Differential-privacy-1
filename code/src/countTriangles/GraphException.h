#ifndef GRAPHEXCEPTION_H
#define GRAPHEXCEPTION_H

#include <exception>

/**
 * Exceptions that can be thrown when executing graph operations
 */
namespace Graph1
{
	class GraphException : public std::exception
	{
	public:
		GraphException(const char* what = "Unknown graph error"):
			_what(what)
		{};
		inline const char* what() const noexcept { return _what; };
	private:
		const char* _what;
	};
	
	class EdgeNotFoundException : public GraphException
	{
	public:
		EdgeNotFoundException(const char* what = "Edge not found"):
			GraphException(what)
		{};
	};
	
	class UnsupportedOperationException : public GraphException
	{
	public:
		UnsupportedOperationException(const char* what = "Unsupported operation for this type of graph"):
			GraphException(what)
		{};
	};
}
#endif /* GRAPHEXCEPTION_H */
