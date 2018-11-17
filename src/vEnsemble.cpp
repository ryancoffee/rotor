#include <iterator>
#include <functional>


vEnsemble::vEnsemble(size_t nvibsin);
: m_nvibs(nvibsin)
{
}
vEnsemble::~vEnsemble(void)
{
}

bool vEnsemble::checkgrowth(void){
	const float buffer(0.05); // setting 5% of length as the buffer size, we will grow the vector by 2x of this 
	size_t bufferi;
	size_t oldsz = pv.size();
	bufferi = size_t(oldsz * (1.- buffer));
	if ( pv[bufferi] > m_pthresh ) {
		size_t newsize = (oldsz * (1. + 2.*buffer));
		pv.resize(newsize,float(0));
		ev.resize(newsize);
		m_molPtr->fill(*this,oldsz);
		return true;
	}
	return false;
}

