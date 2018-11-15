#ifndef VENSEMBLE_H
#define VENSEMBLE_H

#include <algorithm>
#include <fstream>

#include <DataOps.hpp>
#include <Constants.hpp>
#include <Molecules.hpp>

using DataOps::operator<<;

class Molecules;

class vEnsemble {

//friend class Molecules;

public:

	vEnsemble(size_t nvibsin,Molecules * molPtrin = nullptr);
	~vEnsemble(void);
	void setkT(float kTinauin);
	void setmolecule(Molecules * molPtrin);
	size_t limitpops(double thresh = std::nextafter(double(0),double(1)));
	size_t initdist(void);
	double getpop(void);
	void normalize(void);
	inline void printdist(std::ofstream & of) {of << pv;}
	bool checkgrowth(void);


private:

	Molecules * m_molPtr;
	
	size_t setenergies(std::string & molstring);

	size_t v;
	size_t maxv;

	double kTinauin;
	std::vector<double> pv;
	std::vector<double> pvnew;
	std::vector<double> ev;

};

#endif
