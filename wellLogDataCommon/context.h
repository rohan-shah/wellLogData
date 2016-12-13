#ifndef CONTEXT_HEADER_GUARD
#define CONTEXT_HEADER_GUARD
#include <boost/noncopyable.hpp>
namespace wellLogData
{
	class context : public boost::noncopyable
	{
	public:
	private:
		context();
		context(const context&);
	};
}
#endif