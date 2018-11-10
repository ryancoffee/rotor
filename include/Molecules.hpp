#ifndef MOLECULES_H
#define MOLECULES_H

#include <Constants.hpp>
#include <jEnsemble.hpp>
#include <vEnsemble.hpp>

class Molecules {

	static const enum MolID {nno,oo,ii,nno,oco};

	public:
	Molecules(float kTin = float(50) * Constants::kb<float>()/Constants::Eh<float>(), MolID in = nno);
	~Molecules(void);

template <typename DType>
	bool fill(vEnsemble & vens);
template <typename DType>
	bool fill(jEnsemble & jens);
template <typename DType>
	bool fill(jEnsemble & jens,const unsigned oldsz);
	inline void setkT(float kTin) {m_kT = kTin * Constants::kb<float>()/Constants::Eh<float>();}

	std::string & getMolID(void);
	std::string & setMolID(MolID idin = nno);

	//	INLINES		//
template <typename DType>
	inline DType delta_alpha(void) { return delta_alpha(this->m_id); }
template <typename DType>
	inline DType delta_alpha(MolID id) {
		switch (id) {
			case nno: return DType(6.98374727360818)*Constants::aupolarizability<DType>();
			case nn: return DType(.696)*Constants::aupolarizability<DType>();
			case ii: return DType(7)*Constants::aupolarizability<DType>();
			default: return DType(.5)*Constants::aupolarizability<DType>();
		} 
	}

	private:
	MolID m_id;
	float m_kT;
	
};

#endif
