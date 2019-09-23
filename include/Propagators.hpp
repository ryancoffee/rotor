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
namespace boost_ublas = boost::numeric::ublas;
typedef typename boost_ublas::vector< std::complex <double> > cvec_t;
typedef typename boost_ublas::matrix< std::complex <double> > cmat_t;

class jEnsemble;
class PulseTime;

class jFreePropagator
{
	friend class PulseTime;
	friend class jEnsemble;

	public:
		jFreePropagator(jEnsemble & jens,const double dtinau);
		jFreePropagator(jEnsemble & jens,const double dtinau);
		~jFreePropagator();

		inline double apply(jEnsemble & jens,cvec_t & yin, double &t);
		inline double apply(jEnsemble & jens,cvec_t & yin, double &t, const double & delta_tin);
		bool build(jEnsemble & jens);
		bool build(jEnsemble & jens,const double & dt);

	private:
		size_t dim;
		cvec_t Uvec;

		double m_dt;
};


class jKickPropagator
{
	public:
		jKickPropagator(jEnsemble & jens,const PulseTime & pulse);
		jKickPropagator(jEnsemble & jens,const double & pulseduration)
		~jKickPropagator();

		inline bool apply(double &t, cvec_t &yin)
		{
			try {
				boost_ublas::hermitian_adaptor< cmat_t , boost_ublas::lower> hal (Umat);
				yin = boost_ublas::prod(hal,yin);
				t += kickstepsize;
			} catch (std::exception &e) {
				std::cerr << e.what() << "\n" << std::flush;
				return false;
			}
			return true;
		}

		bool build(jEnsemble & jens,PulseTime & pulse);

	private:
		size_t dim;
		cmat_t Umat;
		cvec_t m_state;

		double kickstepsize;

		std::vector< cvec_t > state_History;
		std::vector< double > t_History;

		void printUmat();
		void printcornerUmat();
		void printcornerUmat_real();
		void printcornerUmat_imag();
};



class Kicker 
{
	public:
	Kicker(const jEnsemble &jens,const PulseTime &pulse)
	void operator()(const cvec_t x, cvec_t xddt, const double t);
	private:
	size_t m_dim;

}



#endif
