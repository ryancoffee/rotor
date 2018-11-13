#include <iterator>
#include <functional>


vEnsemble::vEnsemble(unsigned nvibsin,Molecules * molPtrin);
: m_nvibs(nvibsin)
{
	m_molPtr = molPtrin;
}
vEnsemble::~vEnsemble(void)
{

}

void vEnsemble::setkT(float kTinauin)
{
	kTinau = kTinauin;
}

void vEnsemble::checkgrowth(void){
	const float buffer(0.05); // setting 5% of length as the buffer size, we will grow the vector by 2x of this 
	unsigned bufferi;
	unsigned oldsz = pv.size();
	bufferi = unsigned(oldsz * (1.- buffer));
	unsigned newsize;
	if ( pv[bufferi] > m_pthresh ) {
		newsize = unsigned(oldsz * (1. + 2.*buffer));
		pv.resize(newsize,float(0));
		ev.resize(newsize);
		m_molPtr->fill(*this,oldsz);
	}
}
void vEnsemble::setenergies(void){
	m_molPtr->fill(*this);
}

void vEnsemble::setmolecule(Molecules * molPtrin)
{
	m_molPtr = molPtrin;
}

unsigned vEnsemble::limitpops(float thresh = std::nextafter(float(0),float(1)))
{
	//std::replace_if(pv.begin(),pv.end(),std::bind1st(std::greater<DType>(),thresh),std::nextafter(DType(0),DType(-1)));
	//bool abovethresh(DType x){return x>thresh;}
	//std::vector<DType>::iterator it = std::find_end(pv.begin(),pv.end(),abovethresh);
	std::vector<float>::iterator it = std::find_end(pv.begin(),pv.end(),[](float x,float thresh){return bool(x>thresh);});
	unsigned newsize = it - pv.begin() + 1;
	pv.resize(newsize);
	ev.resize(newsize);
	return pv.size();
}

unsigned vEnsemble::initdist(void){
	setenergies();
	pv.assign(ev.begin(),ev.end());
	float kT = m_molPtr->getkT();
	std::transform(pv.begin(), pv.end(), pv.begin(), [float & kT](float x){return std::exp(-x/kT);});
	safe_normalize(pv);
	limitpops(float(1e-2));
	safe_normalize(pv);
	std::clog << "relevent vib pops = " << pv << std::endl;
	return pv.size()
}


