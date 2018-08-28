#ifndef VENSEMBLE_H
#define VENSEMBLE_H

#include <Constants.hpp>
#include <Ensemble.hpp>
#include <algorithm>
#include <Molecules.hpp>

using DataOps::operator<<;
using Constants::kb;
using Constants::Eh;

template <typename DType>
class vEnsemble {
	friend bool Molecules;

public:

	vEnsemble(unsigned nvibsin,const MolID id = Molecules::nno);
	~vEnsemble(void);
	void setkT(DType kTinauin);
	void molecule(Molecules & mol);
	unsigned limitpops(DType thresh = std::nextafter(DType(0),DType(1)));
	unsigned initdist(void);
	float getpop(void);
	void normalize(void);
	inline void printdist(ofstream & ofPtr) {of << pv;}


private;

	Molecules m_mol;
	
	unsigned setenergies(std::string & molstring);

	unsigned v;
	unsigned maxv;

	float kTinauin;
	std::vector<DType> pv;
	std::vector<DType> pvnew;
	std::vector<DType> ev;

};

#endif
