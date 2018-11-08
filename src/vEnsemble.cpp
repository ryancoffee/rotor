#include <iterator>
#include <functional>


vEnsemble::vEnsemble(unsigned nvibsin,Molecules mol)
: m_nvibs(nvibsin)
{
	m_molecule = mol;
}
vEnsemble::vEnsemble(unsigned nvibsin,const MolID id = Molecules::nno)
: m_nvibs(nvibsin)
{
	m_molecule = Molecules(id);
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
		pv.resize(newsize,DType(0));
		ev.resize(newsize);
		m_molecule.fill(*this,oldsz);
	}
}
void vEnsemble::setenergies(void){
	m_molecule.fill(*this);
}

template <typename DType>
void vEnsemble::molecule(Molecules<DType> & mol)
{
	m_molecule = mol;
}
template <typename DType>
unsigned vEnsemble::molecule(MolID & id,DType kTin = DType(1.f/40)*kb()/Eh())
{
	m_molecule = Molecules(id,kTin);
}

template <typename DType>
unsigned vEnsemble::limitpops(DType thresh = std::nextafter(DType(0),DType(1)))
{
	//std::replace_if(pv.begin(),pv.end(),std::bind1st(std::greater<DType>(),thresh),std::nextafter(DType(0),DType(-1)));
	//bool abovethresh(DType x){return x>thresh;}
	//std::vector<DType>::iterator it = std::find_end(pv.begin(),pv.end(),abovethresh);
	std::vector<DType>::iterator it = std::find_end(pv.begin(),pv.end(),[](DType x,DType thresh){return bool(x>thresh);});
	unsigned newsize = it - pv.begin() + 1;
	pv.resize(newsize);
	ev.resize(newsize);
	return pv.size();
}

template <typename DType>
unsigned vEnsemble::initdist(void){
	setenergies();
	pv.assign(ev.begin(),ev.end());
	DType kT = m_molecule.getkT();
	std::transform(pv.begin(), pv.end(), pv.begin(), [DType & kT](DType x){return std::exp(-x/kT);});
	safe_normalize(pv);
	limitpops(DType(1e-2));
	safe_normalize(pv);
	std::clog << "relevent vib pops = " << pv << std::endl;
	return pv.size()
}


