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
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

// my headers
#include <Constants.hpp> // --- conversion constancts lie a0, Eh, icmPau, fsPau etc. --- //
#include <DataOps.hpp>
#include <PulseTime.hpp> // --- this handels the pulse envelope parameters like strength, duration and t0, and returns fales of FF and dFFdt --- //
#include <jEnsemble.hpp>
#include <vEnsemble.hpp>
/*
#include "Members.hpp" // --- This is the colleciton of member functions ( kickfunc and kickjac for now ) --- //
#include "FuncJac.hpp" // --- This defines kickfunc and kickjac

// my preprocessor defs

// gsl_blas defs
#define DAGGAR CblasConjTrans
#define TRANS CblasTrans
#define NOTRANS CblasNoTrans
#define RIGHT CblasRight
#define LEFT CblasLeft
#define UPPER CblasUpper
#define LOWER CblasLower
 */

using namespace std::complex_literals;

class jKickPropagator
{
	public:
		KickPropagator(const size_t dimin);
		~KickPropagator();

		inline bool apply(double &t,std::vector< std::complex<double> > &yin)
		{
			//    clog << "kicker.apply(t,dt,yPTr) is in question" << endl;
			bool wasapplied = false;
			gsl_blas_ztrmv(UPPER,NOTRANS,CblasNonUnit,UmatPtr,yinPtr);
			t += kickstepsize;
			//    clog << "\t\t\t... kicker.apply(t,dt,yPTr) is OK" << endl;
			return wasapplied;
		}

	private:
		bool build();
		size_t dim;
		boost::numeric::ublas::matrix<std::complex<double> > Umat;
		//gsl_matrix_complex *UmatPtr;
		//gsl_vector_complex *yPtr;

		// use boost
		vector_slice<vector<std::complex<double> > > Uslice;
		gsl_vector_complex_view Ucol;
		gsl_vector_complex_view Upart;

		double kickstepsize;

		void printUmat();
		void printcornerUmat();
		void printcornerUmat_real();
		void printcornerUmat_imag();
		void sampleyPtr();
};





class jFreePropagator
{
	public:
		FreePropagator(const double dtinau);
		~FreePropagator();

		double & apply(double &t, std::vector< std::complex<double> > &y);
		double & apply(double &t, const double &delta_tin, std::vector< std::complex< double> > &y);  

	private:
		void build();

		size_t dim;
		std::vector< std::complex<double> > Uvec:

			double dt;
};


#endif
