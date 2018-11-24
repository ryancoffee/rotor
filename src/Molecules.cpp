#include <Molecules.hpp>
#include <algorithm>

Molecules::Molecules(float kT = float(50) * Constants::kb<float>()/Constants::Eh<float>(), MolID id = nno)
: m_id(id)
, m_kT(kT)
{
}

Molecules::MolID & Molecules::getMolID(void)
{
	return m_id;
}

std::string & Molecules::getMolIDstring(void)
{
	std::string out("MolID_");
	switch (this->m_id) {
		case nno:
			out += "nno";	
		case oco:
			out += "oco";	
		case nn:
			out += "nn";	
		case oo:
			out += "oo";	
		case ii:
			out += "ii";
		default:
			out += "none";
	}
	out += "_" + std::to_string(this->m_id);
	return out;
}


/*
		Filling vibrational energies
*/

template <typename DType>
bool Molecules::fill(vEnsemble & vens)
{
	std::vector<DType> ens;
	std::vector<DType> cc;
	switch (this->m_id){
		case nno :
			// From Bohlin et al., J. Raman Spect. v43 p604 (2012)
			// units are icm as usual
			ens = {0.0, 588.767870,588.767870 ,1168.13230 ,1177.74467 ,1177.74467 ,1284.90334 ,1749.06523 ,1749.06515 ,1766.91238 ,1766.91224 ,1880.26574 ,1880.26574 ,2223.75676 ,2322.57308 ,2331.12151 ,2331.12145 ,2356.25242 ,2461.99644 ,2474.79870 ,2474.79865 ,2563.33944};
		case oco :
			// from Rothman and Young, J. Quonf. Spectrosc. Radia. Transfer Vol. 25, pp. 505-524. 1981
			ens = {0.0,667.379,1285.4087,1335.129,1388.1847,1932.472,2003.244,2076.855,2349.1433, 2548.373,2585.032,2671.146,2671.716,2760.735,2797.140,3004.012,3181.450,3240.564,3339.340,3340.501,3442.256,3500.590,3612.842,3659.271,3714.783 };
		case nn :
			ens = {1175.5, 3505.2, 5806.5, 8079.2, 10323.3, 12538.8, 14725.4, 16883.1, 19011.8, 21111.5, 23182.0};
		case oo :
			// from JOURNAL OF MOLECULAR SPECTROSCOPY 154,372-382 ( 1992) G. ROUILLE "High-Resolution Stimulated Raman Spectroscopy of O2"
			ens = {1556.38991/2.0, 1556.38991, 1532.86724, 1509.5275};
			for (size_t v=1;v<ens.size();++v){ ens[v] += ens[v-1]; }

		case ii :
			cc = {214.5481, -0.616259, 7.507e-5, -1.263643e-4, 6.198129e-6, -2.0255975e-7, 3.9662824e-9, -4.6346554e-11, 2.9330755e-13, -7.61000e-16};
			ens.resize(cc.size(),DType(0));
			for (size_t v=0;v<ens.size();++v){
				for(size_t n=0;n<cc.size();n++){
					ens[v] += cc[n] * std::pow(DType(v + 0.5) , int(n + 1));
				}
			}
		default :
			std::cerr << "Failed to choose molecule in Molecules::fill(vEnsemble &... method" << std::endl;
			return false;

	}

	if (vens.ev.size() > ens.size()){ vens.ev.resize(ens.size()); }
	std::copy_n(ens.begin,vens.ev.size(),std::back_inserter(vens.ev));
	vens.ev /= Constants::icmPau<DType>();
	return true;
}

template <typename DType>
DType & Molecules::jmultiplicity(const size_t j)
{
	switch (this->m_id){
		case nno:
		return 1;
		case oco:
		return (j%2==0? 0 : 1);
		case oo:
		return (j%2==0? 0 : 1);
		case ii:
		return (j%2==0? 7 : 5);
		case nn:
		return (j%2==0? 2 : 1);
		default:
		return 1;
		
	}
	// even odd intensity ratio = mul 2:1 for n2, 7:5 for i2, 0:1 for o2 I think
	// I co2 think this is same as O2 but maybe opposite if electronic is symmetric for CO2 // even odd intensity ratio = mul 2:1 for n2, 7:5 for i2, 
	// 0:1 for o2 I think, for N2O this should be an even 1:1 ratio
}


/*
		Filling rotational energies
*/

template <typename DType>
bool Molecules::fill(jEnsemble & jens)
{
	std::vector<DType> Bv;
	std::vector<DType> Dv;
	std::vector<DType> Hv;
	size_t v;

	switch (this->m_id){
		case nno:
			// From Bohlin et al., J. Raman Spect. v43 p604 (2012)
			// units are icm as usual
			Bv = {0.419011 , 0.419177 , 0.419969 , 0.419920 , 0.420125 , 0.420126 , 0.417255 , 0.419583 , 0.421079 , 0.420667 , 0.420671 , 0.417464 , 0.418372 , 0.415559 , 0.420618 , 0.420768 , 0.420772 , 0.421218 , 0.418147 , 0.418530 , 0.418531 , 0.415605 };
			Dv = {1.76e-7 , 1.78e-7 , 1.79e-7 , 2.49e-7 , 1.19e-7 , 1.18e-7 , 1.72e-7 , 2.11e-7 , 2.17e-7 , 1.61e-7 , 1.68e-7 , 1.74e-7 , 1.71e-7 , 1.75e-7 , 3.96e-7 , 0.17e-7 , 2.16e-7 , 2.72e-7 , 2.43e-7 , 1.20e-7 , 1.75e-7 , 1.63e-7 };
			Hv = {0.16e-13, -0.17e-13, -0.17e-13 , 29.55e-13, -29.50e-13, 0.95e-13, 1.46e-13, 12.22e-13, -3.59e-13, -9.91e-13, 30.15e-13, 1.07e-13, 2.17e-13, -0.13e-13, 140.28e-13, -142.92e-13, 4.83e-13, 1712.20e-13, 25.90e-13, -26.98e-13, 2.44e-13, 5.68e-13};

			assert(Bv.size() == Dv.size());
			assert(Bv.size() == Hv.size());
			v = std::min(Bv.size()-1,jens.v);
			for (size_t j=0;j<jens.ej.size();++j){
				jens.ej[j] = 
					( Bv[v]*j*(j+1)
					  - Dv[v]*std::pow(j,int(2))*std::pow(j+1,int(2)) 
					  + Hv[v]*std::pow(j,int(3))*std::pow(j+1,int(3))
					) / Constants::icmPau<DType>(); // in atomic units
			}

		case oco:
			// from Rothman and Young, J. Quonf. Spectrosc. Radia. Transfer Vol. 25, pp. 505-524. 1981
			Bv = {0.39021894,0.39064230,0.39048230,0.39167020,0.39018893,0.39073215,0.39238558,0.39041600,0.38714140,0.39110670,0.39193800,0.38954820,0.39308410,0.39153500,0.39059100,0.38759300,0.3910280,0.3926960,0.3900350,0.393908,0.3922100,0.3904610,0.38750493,0.38864000,0.38706227};
			Dv = {1.33373e-7,1.359e-7,1.57161e-7,1.389e-7,1.14952e-7,1.441e-7,1.403e-7,1.281e-7,1.33034e-7,1.7820e-7,1.390e-7,1.2630e-7,1.42e-7,1.44e-7,0.88e-7,1.349e-7,1.63e-7,1.51e-7,1.37e-7,1.44e-7,1.36e-7,1.14e-7,1.58150e-7,1.37445e-7,1.13570e-7};
			Hv = {0.16e-13,0.17e-13,2.33e-13,0.,1.91e-13,0.,0.,0.,0.17e-13,0.40e-13,0.,4.65e-13,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,2.74e-13,0.,1.10e-13};

			assert(Bv.size() == Dv.size());
			assert(Bv.size() == Hv.size());
			v = std::min(Bv.size()-1,jens.v);
			for (size_t j=0;j<jens.ej.size();j++){
				jens.ej[j] = 
					( Bv[v]*j*(j+1)
					- Dv[v]*std::pow(j,int(2))*std::pow(j+1,int(2)) 
					+ Hv[v]*std::pow(j,int(3))*std::pow(j+1,int(3))
					) / Constants::icmPau<DType>(); // in atomic units
			}

		case nn:
			// numbers come from Loftus and Kuprienie, J. Phys. Chem. ref. Data, Vol. 6, No. 1, 1977. p242
			Bv = {1.98957, 1.972, 1.9548, 1.9374, 1.9200, 1.9022, 1.8845, 1.8666, 1.8488, 1.8310, 1.8131, 1.7956, 1.7771, 1.7590, 1.7406, 1.7223};
			Dv = {1e-6*5.75};
			v = std::min(Bv.size()-1,jens.v);
			for (size_t j=0;j<jens.ej.size();++j){
				jens.ej[j] = 
					( Bv[v]*j*(j+1)
					  - Dv[0]*std::pow(j,int(2))*std::pow(j+1,int(2))
					) / Constants::icmPau<DType>(); // in atomic units
			} 
		case oo:
			// from JOURNAL OF MOLECULAR SPECTROSCOPY 154,372-382 ( 1992) G. ROUILLE "High-Resolution Stimulated Raman Spectroscopy of O2"
			Bv = {1.437676476, 1.42186454, 1.4061199, 1.39042};
			Dv = {4.84256e-6, 4.8418e-6, 4.8410e-6, 4.8402e-6};
			Hv = {2.8e-12};
			v = std::min(Dv.size()-1,jens.v);
			for (size_t j=0;j<jens.ej.size();++j){
				jens.ej[j] = 
					( Bv[v]*j*(j+1)
					- Dv[v]*std::pow(j,int(2))*std::pow(j+1,int(2)) 
					+ Hv[0]*std::pow(j,int(3))*std::pow(j+1,int(3))
					) / Constants::icmPau<DType>(); // in atomic units
			}
		case ii:
			Bv = {	3.7395e-2 
				- 1.2435e-4*(v+0.5) 
				+ 4.498e-7*std::pow(DType(v)+0.5,int(2)) 
				- 1.482e-8*std::pow(DType(v)+0.5,int(3)) 
				- 3.64e-11*std::pow(DType(v)+0.5,int(4)) 
			};
			Dv = {	4.54e-9 
				+ 1.7e-11*(DType(v)+0.5) 
				+ 7e-12*std::pow(DType(v)+0.5,int(2))
			};
			for (size_t j=0;j<jens.ej.size();++j){
				jens.ej[j] = ( Bv[0]*j*(j+1) - Dv[0]*std::pow(j,int(2))*std::pow(j+1,int(2)) ) / Constants::icmPau<DType>() ; // in atomic units
			}
		default:
			std::cerr << "Failed to choose molecule in Molecules::fill(jEnsemble &... method" << std::endl;
			return false;

	}
	return true;
}

template <typename DType>
size_t Molecules::initdist(jEnsemble & jens){
	size_t  j = 0;
	DType kT = m_kT;
	jens.pj.assign(jens.ej.begin(),jens.ej.end());
	std::transform(jens.pj.begin(), jens.pj.end(), jens.pj.begin(), 
			[&j,&kT](DType x){
				DType val = (x/kT>DType(0.1) ? std::exp(-x/kT) : 1+std::expm1(-x/kT) );
				val *= jmultiplicity<DType>(j)(2*(j)+1);
				j++;
				return val;
				}
		      );
	DataOps::safe_normalize(jens.pj);
	return jens.pj.size();
}

template <typename DType>
size_t Molecules::initdist(vEnsemble & vens){
	DType kT = m_kT;
	vens.pv.assign(vens.ev.begin(),vens.ev.end());
	std::transform(vens.pv.begin(), vens.pv.end(), vens.pv.begin(), 
				[kT](DType x){
					return (x/kT>DType(0.1)? std::exp(-x/kT) : 1+ std::expm1(-x/kT) );
				});
	DataOps::safe_normalize(vens.pv);
	vens.limitpops(DType(1e-2));
	DataOps::safe_normalize(vens.pv);
	std::clog << "relevent vib pops = " << vens.pv << std::endl;
	return vens.pv.size();
}

template <typename DType>
size_t Molecules::limitpops(jEnsemble & jens, DType thresh)
{ 
	auto it = std::find_end(jens.pj.begin(),jens.pj.end(),[](DType x,DType thresh){return bool(x>thresh);});
	//std::vector<DType>::iterator it = std::find_end(jens.pj.begin(),jens.pj.end(),[](DType x,DType thresh){return bool(x>thresh);});
	size_t newsize = it - jens.pj.begin() + 1;
	jens.pj.resize(newsize);
	jens.ej.resize(newsize);
	return jens.pj.size();
}
template <typename DType>
size_t Molecules::limitpops(vEnsemble & vens,DType thresh)
{
	auto it = std::find_end(vens.pv.begin(),vens.pv.end(),[](DType x,DType thresh){return bool(x>thresh);});
	//std::vector<DType>::iterator it = std::find_end(vens.pv.begin(),vens.pv.end(),[](DType x,DType thresh){return bool(x>thresh);});
	size_t newsize = it - vens.pv.begin() + 1;
	vens.pv.resize(newsize);
	vens.ev.resize(newsize);
	return vens.pv.size();
}

