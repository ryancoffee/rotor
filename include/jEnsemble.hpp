#ifndef JENSEMBLE_H
#define JENSEMBLE_H

#include <vector>
#include <iostream>
#include <fstream>
#include <limits>

#include <Constants.hpp>
#include <DataOps.hpp>
#include <Molecules.hpp>

using DataOps::operator<<;


class jEnsemble {

friend class Molecules;

public:

	jEnsemble(size_t njsin,size_t vin);
	~jEnsemble(void);

	inline int setm(const int min){m=min; return m;}
	inline size_t setv(const size_t vin){v=vin; return v;}

	double getpop(void);
	void normalize(void);
	size_t limitpops(double thresh);

	size_t initdist(void);
        inline void printdist(std::ofstream & of) { of << pj;}


private:
	size_t m_realj;
	size_t v;
	int m;
	size_t m_njs;

	std::vector<double> pj;
	std::vector<double> ej;
        std::vector<double> aajm;
	std::vector<double> bbjm;
	std::vector<double> ccjm;


protected:
	bool checkgrowth(void);
	void fill6j(const size_t oldsz = 0);
};

#endif  

