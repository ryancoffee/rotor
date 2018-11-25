#include <Constants.hpp> // for kb and Eh units
#include <jEnsemble.hpp>
#include <algorithm>

using Constants::kb;
using Constants::Eh;


jEnsemble::jEnsemble(size_t njsin = 50,size_t vin = 0)
: m_njs(njsin)
, v(vin)
{
        ej.resize(m_njs,double(0));
        pj.resize(m_njs,double(0));
	aajm.resize(m_njs,double(0));
	bbjm.resize(m_njs,double(0));
	ccjm.resize(m_njs,double(0));
}       
jEnsemble::~jEnsemble(void)
{
}

bool jEnsemble::checkgrowth(void)
{
	const float buffer(0.1); // setting 10% of length as the buffer size, we will grow the vector by 2x of this 
	unsigned bufferi;
	unsigned oldsz = pj.size();
	bufferi = unsigned(oldsz * (1.- buffer));
	unsigned newsize;
	if ( pj[bufferi] > m_pthresh ) {
		newsize = unsigned(oldsz * (1. + 2.*buffer));
		pj.resize(newsize,double(0));
		ej.resize(newsize);
		aajm.resize(newsize,double(0));
		bbjm.resize(newsize,double(0));
		ccjm.resize(newsize,double(0));
		fill6js(oldsz);
		return true;
	}
	return false;
}

void jEnsemble::fill6js(const size_t oldsz)
{
	aajm.resize(pj.size(),double(0));
	bbjm.resize(pj.size(),double(0));
	ccjm.resize(pj.size(),double(0));

	// 	these coefficients come from Arfken and Webber 4th ed. p753, eq.12.189 applied twice.
	for(size_t j=std::min(m,int(oldsz));j<aajm.size();++j){
		aajm[j] = double((j-m+1)*(j+m+1)*(2*j-1) + (2*j+3)*(j+m)*(j-m))/double((2*j-1)*(2*j+1)*(2*j+3));
		if (j>1) { // this term is 0 when j=0 m=0, j=1 m=0, j=1 m=1, and j=1 m=-1
			bbjm[j] = std::sqrt(double((j-m)*(j+m)*(j-m-1)*(j+m-1)) / double((2*j+1)*(2*j-3)*std::pow(int(2*j-1),int(2))) );  
		}
		ccjm[j] = std::sqrt(double((j-m+1)*(j+m+1)*(j-m+2)*(j+m+2)) / double((2*j+1)*(2*j+5)*std::pow(int(2*j+3),int(2))));  

	}
}
