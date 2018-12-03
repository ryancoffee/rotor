#include <vEnsemble.hpp>
#include <iterator>
#include <functional>


vEnsemble::vEnsemble(size_t nvibsin)
: m_nvs(nvibsin)
{
        ev.resize(m_nvs,double(0));
        pv.resize(m_nvs,double(0));
}
vEnsemble::~vEnsemble(void)
{
}

bool vEnsemble::checkgrowth(void){
	const float buffer(0.1); // setting 10% of length as the buffer size, we will grow the vector by 2x of this 
	size_t bufferi;
	size_t oldsz = pv.size();
	bufferi = size_t(oldsz * (1.- buffer));
	if ( pv[bufferi] > m_pthresh ) {
		size_t newsize = (oldsz * (1. + 2.*buffer));
		pv.resize(newsize,float(0));
		ev.resize(newsize);
		return true;
	}
	return false;
}

