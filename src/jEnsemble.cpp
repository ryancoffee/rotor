#include <Constants.hpp> // for kb and Eh units
#include <jEnsemble.hpp>
#include <Ensemble.hpp>
#include <algorithm>

using Constants::kb;
using Constants::Eh;


jEnsemble::jEnsemble(unsigned njsin = 50,unsigned vin = 0,Molecules * molPtrin = nullptr)
: m_njs(njsin)
, m_popthresh(10*std::numeric_limits<double>::min())
, m_v(vin)
, m_molPtr(molPtrin)
{
	m_molPtr = new Molecules(*mol);
        ej.resize(m_njs,double(0));
        pj.resize(m_njs,double(0));
	mol.fill(*this);
	aajm.resize(m_njs,double(0));
	bbjm.resize(m_njs,double(0));
	ccjm.resize(m_njs,double(0));
}       
jEnsemble::~jEnsemble(void)
{
	if (m_molPtr~=nullptr){
		delete m_molePtr;
		m_molPtr = nullptr;
	}
}

unsigned jEnsemble::limitpops(double thresh = std::nextafter(double(0),double(1)))
{ 
	std::vector<double>::iterator it = std::find_end(pv.begin(),pv.end(),[](double x,double thresh){return bool(x>thresh);});
	unsigned newsize = it - pj.begin() + 1;
	pj.resize(newsize);
	ej.resize(newsize);
	return pv.size();
}

unsigned jEnsemble::initdist(void){
	setenergies();
	unsigned j = 0;
	float kT = m_molPtr->getkT();
	pj.assign(ej.begin(),ej.end());
	std::transform(pv.begin(), pv.end(), pv.begin(), 
			[float & kT,unsigned &j,Molecules * m_molPtr](float x){
			val = (x/kt>float(0.1) ? std::exp(-x/kT) : 1+std::expm1(-x/kT) );
				val *= m_molPtr->multiplicity(j)(2*(j)+1);
				j++;
				return val;
				}
		      );
	safe_normalize(pj);
	return pj.size()
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
