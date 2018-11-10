#ifndef JENSEMBLE_H
#define JENSEMBLE_H

#include <vector>
#include <iostream>
#include <fstream>
#include <Constants.hpp>
#include <DataOps.hpp>
#include <Molecules.hpp>
#include <limits>

using DataOps::operator<<;

class Molecules;

class jEnsemble {

	friend class Molecules;

public:

	jEnsemble(unsigned njsin,unsigned vin, const float kTinau);
	jEnsemble(unsigned njsin = 50,unsigned vin = 0,Molecules * mol);
	~jEnsemble(void);

	void setm(const unsigned m_in){m=m_in;}
	void setv(const unsigned v_in){v=v_in;}
	void setmolecule(Molecules & molin){m_moleculePtr = new Molecules(molin);}

        inline void printdist(std::ofstream & of) { of << pj;}
	inline double setpopthresh(double in){m_popthresh = in;return m_popthreash;}


private:
	unsigned m_realj;
	unsigned v;
	int m;
	unsigned m_njs;
	Molecules * m_moleculePtr;

	std::vector<double> pj;
	std::vector<double> ej;
        std::vector<double> aajm;
	std::vector<double> bbjm;
	std::vector<double> ccjm;


protected:
	double m_popthresh;//(10*std::numeric_limits<double>::min()); // used in testing for population leaking up past the 10% of end of the vector	
	void checkgrowth(void);
	void fill6j(const unsigned oldsz = 0);
};

#endif  

