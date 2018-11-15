#include <iterator>
#include <functional>


vEnsemble::vEnsemble(size_t nvibsin,Molecules * molPtrin);
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
void vEnsemble::setenergies(void){
	m_molPtr->fill(*this);
}

void vEnsemble::setmolecule(Molecules * molPtrin)
{
	m_molPtr = molPtrin;
}

size_t vEnsemble::limitpops(double thresh = std::nextafter(double(0),double(1)))
{
	//std::replace_if(pv.begin(),pv.end(),std::bind1st(std::greater<double>(),thresh),std::nextafter(double(0),double(-1)));
	//bool abovethresh(double x){return x>thresh;}
	//std::vector<double>::iterator it = std::find_end(pv.begin(),pv.end(),abovethresh);
	std::vector<double>::iterator it = std::find_end(pv.begin(),pv.end(),[](double x,double thresh){return bool(x>thresh);});
	size_t newsize = it - pv.begin() + 1;
	pv.resize(newsize);
	ev.resize(newsize);
	return pv.size();
}

unsigned vEnsemble::initdist(void){
	setenergies();
	pv.assign(ev.begin(),ev.end());
	double kT = m_molPtr->getkT();
	std::transform(pv.begin(), pv.end(), pv.begin(), [double & kT](double x){return std::exp(-x/kT);});
	safe_normalize(pv);
	limitpops(double(1e-2));
	safe_normalize(pv);
	std::clog << "relevent vib pops = " << pv << std::endl;
	return pv.size()
}


