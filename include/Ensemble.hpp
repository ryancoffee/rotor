#ifndef ENSEMBLE_H
#define ENSEMBLE_H

#include <Constants.hpp> // for kb and Eh

#include <iostream>
#include <fstream>
#include <vector>

// base for ensemble classes //

class Ensemble {

public:
	Ensemble(const unsigned nstatesin);
	virtual ~Ensemble(void){}

	void setkTinau(const float kTinauin);
	void setkT(const float kTin);

	virtual void initdist(void);
        float getpop(void);
        virtual void normalize(void);
        virtual void printdist(std::ofstream & of);
	
	void checkgrowth(void); // this should be something that checks if the state .95 away from the end has nonzero population, then grow the vector length by 5%
	virtual double getstateenergy(const unsigned statenum);


private:
	virtual unsigned setenergies(std::string & molstring);

	float kTinau;
	const unsigned nstates;
	std::vector<float> pops;
	std::vector<double> energies;
};

#endif
