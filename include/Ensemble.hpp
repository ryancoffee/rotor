#ifndef ENSEMBLE_H
#define ENSEMBLE_H

#include <iostream>
#include <fstream>
#include <vector>
#include <Molecules>

// base for ensemble classes //

class Ensemble {

public:
	Ensemble(const unsigned nstatesin,const float kTinkelvinin);
	Ensemble(const unsigned nstatesin);
	virtual ~Ensemble(void){}

	void setkTinau(const float kTinauin);
	void setkT(const float kTin);

	virtual void initdist(void);
        float getpop(void);
        virtual void normalize(void);
        virtual void printdist(std::ofstream & of);
	virtual void setmolecule(Molecules * molPtrin);
	
	void checkgrowth(void); // this should be something that checks if the state .95 away from the end has nonzero population, then grow the vector length by 5%
	virtual double getstateenergy(const unsigned statenum);


private:
	virtual unsigned setenergies(std::string & molstring);

	Molecules * m_molPtr;
	float m_kTinau;
	unsigned m_nstates;
	std::vector<double> m_pops;
	std::vector<double> m_energies;
};

#endif
