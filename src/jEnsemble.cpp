#include <Constants.hpp> // for kb and Eh units
#include <jEnsemble.hpp>
#include <Ensemble.hpp>
#include <algorithm>

using Constants::kb;
using Constants::Eh;


template <typename DType>
jEnsemble::jEnsemble(unsigned njsin,unsigned vin, const DType kTinau,const MolID id = Molecules::nno)
: m_njs(njsin)
, v(vin)
{
	m_molecule = Molecules(id);
        ej.resize(njs,DType(0));
        pj.resize(njs,DType(0));
	aajm.resize(njs,DType(0));
	bbjm.resize(njs,DType(0));
	ccjm.resize(njs,DType(0));
}       
template <typename DType>
jEnsemble::jEnsemble(unsigned njsin,unsigned vin,Molecules & mol)
: m_njs(njsin)
, v(vin)
{
	m_molecule = mol;
        ej.resize(m_njs,DType(0));
        pj.resize(m_njs,DType(0));
	mol.fill(*this);
	aajm.resize(m_njs,DType(0));
	bbjm.resize(m_njs,DType(0));
	ccjm.resize(m_njs,DType(0));
}       

template <typename DType>
unsigned jEnsemble::limitpops(DType thresh = std::nextafter(DType(0),DType(1)))
{ 
	std::vector<DType>::iterator it = std::find_end(pv.begin(),pv.end(),[](DType x,DType thresh){return bool(x>thresh);});
	unsigned newsize = it - pj.begin() + 1;
	pj.resize(newsize);
	ej.resize(newsize);
	return pv.size();
}

	/*
	template <typename IteratorT, typename FunctionT>
FunctionT jEnsemble::fillpop(IteratorT first, 
		IteratorT last, 
		typename std::iterator_traits<IteratorT>::difference_type initial,
		FunctionT func)
{
	for (;first != last; ++first, ++initial)
		func(initial, *first);
	return func;
}
*/
template <typename DType>
unsigned jEnsemble::initdist(void){
	setenergies();
	unsigned j = 0;
	DType kT = m_molecule.getkT();
	pj.assign(ej.begin(),ej.end());
	std::transform(pv.begin(), pv.end(), pv.begin(), 
			[DType & kT,unsigned &j,Molecules & m_molecule](DType x){
			val = (x/kt>DType(0.1) ? std::exp(-x/kT) : 1+std::expm1(-x/kT) );
				val *= m_molecule.multiplicity(j)(2*(j)+1);
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
		pj.resize(newsize,DType(0));
		ej.resize(newsize);
		aajm.resize(newsize,DType(0));
		bbjm.resize(newsize,DType(0));
		ccjm.resize(newsize,DType(0));
		mol.fill(*this,oldsz);
		fill6j(oldsz);
	}
}

template <typename DType>
void jEnsemble::fill6js(const unsigned oldsz = 0)
{
	aajm.resize(pj.size(),DType(0));
	bbjm.resize(pj.size(),DType(0));
	ccjm.resize(pj.size(),DType(0));

	// 	these coefficients come from Arfken and Webber 4th ed. p753, eq.12.189 applied twice.
	for(unsigned j=std::min(m,oldsize);j<aajm.size();++j){
		aajmPtr[j] = DType((j-m+1)*(j+m+1)*(2*j-1) + (2*j+3)*(j+m)*(j-m))/DType((2*j-1)*(2*j+1)*(2*j+3));
		if (j>1) { // this term is 0 when j=0 m=0, j=1 m=0, j=1 m=1, and j=1 m=-1
			bbjmPtr[j] = std::sqrt(DType((j-m)*(j+m)*(j-m-1)*(j+m-1)) / DType((2*jj+1)*(2*jj-3)*std::power(int(2*jj-1),int(2))) );  
		}
		ccjmPtr[jj] = std::sqrt(DType((jj-m+1)*(jj+m+1)*(jj-m+2)*(jj+m+2)) / DType((2*jj+1)*(2*jj+5)*std::power(int(2*jj+3),int(2))));  

	}
}
