#ifndef VENSEMBLE_H
#define VENSEMBLE_H

#include <algorithm>
#include <fstream>
#include <limits>

#include <DataOps.hpp>
#include <Constants.hpp>
#include <Molecules.hpp>

using DataOps::operator<<;

class Molecules;
/*
bool Molecules::fill(vEnsemble & vens, const size_t fillstart);
size_t Molecules::initdist(vEnsemble & vens);
size_t Molecules::updatedist(vEnsemble & vens, const size_t start);
size_t Molecules::limitpops(vEnsemble & vens, const double thresh);
*/

class vEnsemble {
	/*
	friend bool Molecules::fill(vEnsemble & vens, const size_t fillstart);
	friend size_t Molecules::initdist(vEnsemble & vens);
	friend size_t Molecules::updatedist(vEnsemble & vens, const size_t start);
	friend size_t Molecules::limitpops(vEnsemble & vens, const double thresh);
	*/
	friend class Molecules;
public:

	vEnsemble(size_t nvibsin);
	~vEnsemble(void);

	size_t limitpops(double thresh = std::nextafter(double(0),double(1)));
	size_t initdist(void);
	double getpop(void);
	void normalize(void);
	bool checkgrowth(void);
	inline size_t size(void){return m_nvs;}
	inline void printdist(std::ofstream & of) {of << pv;}
	inline double set_pthresh(const double th=0.05){m_pthresh = th;return m_pthresh;}

private:
	size_t m_nvs;
	double m_pthresh;

	std::vector<double> pv;
	std::vector<double> ev;

};

#endif
