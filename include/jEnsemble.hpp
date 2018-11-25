#ifndef JENSEMBLE_H
#define JENSEMBLE_H

#include <vector>
#include <iostream>
#include <fstream>
#include <limits>

#include <Constants.hpp>
#include <DataOps.hpp>
#include <Molecules.hpp>
#include <Propagators.hpp>

using DataOps::operator<<;


class jEnsemble {

public:

	jEnsemble(size_t njsin,size_t vin);
	~jEnsemble(void);

	inline int setm(const int min){m=min; return m;}
	inline size_t setv(const size_t vin){v=vin; return v;}
	inline double set_pthresh(const double th=0.05){m_pthresh = th;return m_pthresh;}

	double getpop(void);
	void normalize(void);
	bool checkgrowth(void);
	inline size_t size(void){return m_njs;}

        inline void printdist(std::ofstream & of) { of << pj;}


private:
	double m_pthresh;
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
	void fill6js(const size_t oldsz);
};

#endif  

