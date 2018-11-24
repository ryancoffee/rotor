#ifndef PROPAGATORS_H
#define PROPAGATORS_H

// --- these classes are for the various coupled and free propagators for the rotor system --- //
//
// OK here is the rub.  let's convert to using the boost/numeric/odeint.hpp
// use the adaptive stepper if we are propagating througha  pulse
// use the fixed stepper if we are inside of a pulse
// simple phase accumulation when we are outside of a pulse
// And it looks like we can indeed use complex numbers easily
// Using the sandbox to play with the boost odeint stuff

// standard includes
#include <cmath>
#include <vector>
#include <complex>
#include <exception>

// boost includes
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/hermitian.hpp>
#include <boost/numeric/ublas/io.hpp>

// my headers
#include <Constants.hpp> // --- conversion constancts lie a0, Eh, icmPau, fsPau etc. --- //
#include <DataOps.hpp>
#include <PulseTime.hpp> // --- this handels the pulse envelope parameters like strength, duration and t0, and returns fales of FF and dFFdt --- //
#include <jEnsemble.hpp>

using std::complex_literals::operator""i;
// seems the 
namespace boost_blas = boost::numeric::ublas;

class jEnsemble;
class PulseTime;

class jKickPropagator
{
	friend class PulseTime;
	friend class jEnsemble;
	public:
		jKickPropagator(const size_t dimin);
		~jKickPropagator();

		inline bool apply(double &t,boost_blas::vector< std::complex<double> > &yin)
		{
			try {
			boost_blas::hermitian_adaptor<boost_blas::matrix<std::complex<double> >, boost_blas::lower> hal (Umat);
			auto temp = boost_blas::prod(hal,yin);
			yin = temp;
			t += kickstepsize;
			} catch (std::exception &e) {
				std::cerr << e.what() << "\n" << std::flush;
				return false;
			}
			return true;
		}
		inline bool apply(double &t,std::vector< std::complex<double> > &yin)
		{
			bool wasapplied = false;
			boost_blas::vector<std::complex<double> > y(yin.size(),*yin.data());
			wasapplied = apply(t,y);
			std::copy(y.begin(),y.end(),yin.begin());
			return wasapplied;
		}

	private:
		bool build();
		size_t dim;
		boost_blas::matrix<std::complex<double> > Umat;
		//boost_blas::vector_slice< boost_blas::vector <std::complex <double> > Uslice;

		/*
		// use boost
		gsl_vector_complex_view Ucol;
		gsl_vector_complex_view Upart;
		*/

		double kickstepsize;

		void printUmat();
		void printcornerUmat();
		void printcornerUmat_real();
		void printcornerUmat_imag();
		void sampleyPtr();
};





class jFreePropagator
{
	friend class PulseTime;
	friend class jEnsemble;

	public:
		jFreePropagator(const double dtinau);
		~jFreePropagator();

		bool apply(double &t, boost_blas::vector< std::complex<double> > &y);
		bool apply(double &t, const double &delta_tin, boost_blas::vector< std::complex< double> > &y);  
		bool apply(double &t, std::vector< std::complex<double> > &y);
		bool apply(double &t, const double &delta_tin, std::vector< std::complex< double> > &y);  
		bool build(jEnsemble & jens);
		bool build(jEnsemble & jens,const double & dt);

	private:
		void build();

		size_t dim;
		boost_blas::vector< std::complex<double> > Uvec;

		double m_dt;
};


#endif
