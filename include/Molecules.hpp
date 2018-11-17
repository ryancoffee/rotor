#ifndef MOLECULES_H
#define MOLECULES_H

#include <vector>

#include <Constants.hpp>
#include <vEnsemble.hpp>
#include <jEnsemble.hpp>


class Molecules {


	public:

	enum MolID {nno,nn,oo,ii,oco};

	Molecules(float kTin, MolID in);
	~Molecules(void);

template <typename DType>
	bool fill(vEnsemble & vens);
template <typename DType>
	bool fill(jEnsemble & jens);
template <typename DType>
	DType & jmultiplicity(const unsigned j);
template <typename DType>
	DType kT(void){return m_kT;}

	inline void setkT(float kTin) {m_kT = kTin * Constants::kb<float>()/Constants::Eh<float>();}
template <typename DType>
	size_t initdist(vEnsemble & vens);
template <typename DType>
	size_t initdist(jEnsemble & jens);
template <typename DType>
	size_t limitpops(vEnsemble & vens,DType thresh = std::nextafter(DType(0),DType(1)));
template <typename DType>
	size_t limitpops(jEnsemble & vens,DType thresh = std::nextafter(DType(0),DType(1)));

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
