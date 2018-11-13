#ifndef MOLECULES_H
#define MOLECULES_H

#include <Constants.hpp>
#include <jEnsemble.hpp>
#include <vEnsemble.hpp>


class Molecules {

	enum MolID {nno=0,nn,oo,ii,oco};

	public:
	Molecules(float kTin, MolID in);
	~Molecules(void);

template <typename DType>
	bool fill(vEnsemble & vens);
template <typename DType>
	bool fill(jEnsemble & jens);
template <typename DType>
	DType & jmultiplicity(const unsigned j);

	inline void setkT(float kTin) {m_kT = kTin * Constants::kb<float>()/Constants::Eh<float>();}

	MolID & getMolID(void);
	std::string & getMolIDstring(void);
	std::string & setMolID(MolID idin = nno);

	//	INLINES		//
	inline double delta_alpha(void) 
	{
		switch (m_id) {
			case nno: return double(6.98374727360818)*Constants::aupolarizability<double>();
			case nn: return double(.696)*Constants::aupolarizability<double>();
			case oo: return double(1.14)*Constants::aupolarizability<double>();
			case ii: return double(7)*Constants::aupolarizability<double>();
			case oco: return double(13.2)*Constants::aupolarizability<double>();
			default: return double(.5)*Constants::aupolarizability<double>();
		} 
	};

	private:
	MolID m_id;
	float m_kT;
	
};

#endif
