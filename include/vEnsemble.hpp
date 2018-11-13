#ifndef VENSEMBLE_H
#define VENSEMBLE_H

#include <Constants.hpp>
#include <Ensemble.hpp>
#include <algorithm>
#include <Molecules.hpp>
#include <fstream>

using DataOps::operator<<;

class Molecules;

class vEnsemble {

public:

	vEnsemble(unsigned nvibsin,Molecules * molPtrin = nullptr);
	~vEnsemble(void);
	void setkT(float kTinauin);
	void setmolecule(Molecules * molPtrin);
	unsigned limitpops(float thresh = std::nextafter(float(0),float(1)));
	unsigned initdist(void);
	float getpop(void);
	void normalize(void);
	inline void printdist(std::ofstream & of) {of << pv;}


private:

	Molecules * m_molPtr;
	
	unsigned setenergies(std::string & molstring);

	unsigned v;
	unsigned maxv;

	float kTinauin;
	std::vector<double> pv;
	std::vector<double> pvnew;
	std::vector<double> ev;

};

#endif
