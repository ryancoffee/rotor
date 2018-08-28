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

template <typename DType>
class jEnsemble {

	friend class Molecules;

public:

        jEnsemble(unsigned njsin,const Molecules::MolID id = Molecules::nno);
	jEnsemble(unsigned njsin = 50, Molecules & mol);

	void setm(const unsigned m_in){m=m_in;}
	void setv(const unsigned v_in){v=v_in;}

        inline void printdist(std::ofstream & of) { of << pj;}


private:
	unsigned m_realj;
	unsigned v;
	int m;
	unsigned m_njs;
	Molecules m_molecule;

	std::vector<DType> pj;
	std::vector<DType> ej;
        std::vector<DType> aajm;
	std::vector<DType> bbjm;
	std::vector<DType> ccjm;


protected:
	const DType m_popthresh(10*std::numeric_limits<DType>::min()); // used in testing for population leaking up past the 10% of end of the vector	
	void checkgrowth(void);
	void fill6j(const unsigned oldsz = 0);
	/*
	template <typename IteratorT, typename FunctionT> FunctionT jEnsemble::fillpop(
			IteratorT first,
			IteratorT last,
			typename std::iterator_traits<IteratorT>::difference_type initial,
			FunctionT func)
	*/

};

#endif  

