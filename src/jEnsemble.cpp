#include <Constants.hpp> // for kb and Eh units
#include <jEnsemble.hpp>
#include <algorithm>

using Constants::kb;
using Constants::Eh;


jEnsemble::jEnsemble(size_t njsin = 50,size_t vin = 0)
: m_njs(njsin)
, m_v(vin)
{
        ej.resize(m_njs,double(0));
        pj.resize(m_njs,double(0));
	aajm.resize(m_njs,double(0));
	bbjm.resize(m_njs,double(0));
	ccjm.resize(m_njs,double(0));
	if (molPtrin != nullptr){
		m_molPtr = molPtrin;
		m_molPtr->fill<double>(*this);
	}
}       
jEnsemble::~jEnsemble(void)
{
}

void jEnsemble::checkgrowth(void){
	const float buffer(0.05); // setting 5% of length as the buffer size, we will grow the vector by 2x of this 
	unsigned bufferi;
	unsigned oldsz = pjnew.size();
	bufferi = unsigned(oldsz * (1.- buffer));
	unsigned newsize;
	if ( pj[bufferi] > m_pthresh ) {
		newsize = unsigned(oldsz * (1. + 2.*buffer));
		pj.resize(newsize,double(0));
		ej.resize(newsize);
		aajm.resize(newsize,double(0));
		bbjm.resize(newsize,double(0));
		ccjm.resize(newsize,double(0));
		mol.fill(*this,oldsz);
		fill6j(oldsz);
	}
}

void jEnsemble::fill6js(const unsigned oldsz = 0)
{
	aajm.resize(pj.size(),double(0));
	bbjm.resize(pj.size(),double(0));
	ccjm.resize(pj.size(),double(0));

	// 	these coefficients come from Arfken and Webber 4th ed. p753, eq.12.189 applied twice.
	for(unsigned j=std::min(m,oldsize);j<aajm.size();++j){
		aajmPtr[j] = double((j-m+1)*(j+m+1)*(2*j-1) + (2*j+3)*(j+m)*(j-m))/double((2*j-1)*(2*j+1)*(2*j+3));
		if (j>1) { // this term is 0 when j=0 m=0, j=1 m=0, j=1 m=1, and j=1 m=-1
			bbjmPtr[j] = std::sqrt(double((j-m)*(j+m)*(j-m-1)*(j+m-1)) / double((2*jj+1)*(2*jj-3)*std::power(int(2*jj-1),int(2))) );  
		}
		ccjmPtr[jj] = std::sqrt(double((jj-m+1)*(jj+m+1)*(jj-m+2)*(jj+m+2)) / double((2*jj+1)*(2*jj+5)*std::power(int(2*jj+3),int(2))));  

	}
}
