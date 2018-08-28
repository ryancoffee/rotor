#if ndef MOLECULES_H
#define MOLECULES_H
#include <jEnsemble.hpp>
#include <vEnsemble.hpp>

template <typename DType>
class Molecules {

	static const enum MolID {nno,oo,ii,nno,oco};

	using namespace Constants;

	public:
	Molecules(DType kTin = DType(50) * kb<DType>()/Eh<DType>(), MolID in = nno);
	~Molecules(void);

	bool fill(vEnsemble & vens);
	bool fill(jEnsemble & jens);
	bool fill(jEnsemble & jens,const unsigned oldsz);
	inline void setkT(DType kTin) {m_kT = kTin;}
	std::string & getMolID(void);
	std::string & setMolID(MolID idin = nno);

	//	INLINES		//
	inline DType delta_alpha(void) { return delta_alpha(this->m_id); }
	inline DType delta_alpha(MolID id) {
		switch (id) {
			case nno: return DType(6.98374727360818)*aupolarizability<DType>();
			case nn: return DType(.696)*aupolarizability<DType>();
			case ii: return DType(7)*aupolarizability<DType>();
			default: return DType(.5)*aupolarizability<DType>();
		} 
	}

	private:
	MolID m_id;
	DType m_kT;
	
};

#endif
