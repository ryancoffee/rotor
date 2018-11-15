#ifndef MOLECULES_H
#define MOLECULES_H

#include <vector>

#include <Constants.hpp>
#include <vEnsemble.hpp>


class Molecules {

	enum MolID {nno=0,nn,oo,ii,oco};

	public:
	Molecules(float kTin, MolID in);
	~Molecules(void);

template <typename DType>
	bool fill(vEnsemble & vens);
template <typename DType>
	bool fill(std::vector<DType> & jens);
//	bool fill(jEnsemble & jens);
template <typename DType>
	DType & jmultiplicity(const unsigned j);

	inline void setkT(float kTin) {m_kT = kTin * Constants::kb<float>()/Constants::Eh<float>();}

	MolID & getMolID(void);
	std::string & getMolIDstring(void);
	std::string & setMolID(MolID idin = nno);

	//	INLINES		//
template <typename DType>
	inline DType delta_alpha(void) 
	{
		switch (m_id) {
			case nno: return DType(6.98374727360818)*Constants::aupolarizability<DType>();
			case nn: return DType(.696)*Constants::aupolarizability<DType>();
			case oo: return DType(1.14)*Constants::aupolarizability<DType>();
			case ii: return DType(7)*Constants::aupolarizability<DType>();
			case oco: return DType(13.2)*Constants::aupolarizability<DType>();
			default: return DType(.5)*Constants::aupolarizability<DType>();
		} 
	};

	private:
	MolID m_id;
	float m_kT;
	
};

#endif
