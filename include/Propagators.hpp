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
#include <boost/numeric/ublas/io.hpp>



// my headers
#include <Constants.hpp> // --- conversion constancts lie a0, Eh, icmPau, fsPau etc. --- //
#include <Rotor.hpp> // --- this is the molecule specifics like ej and pj vectors --- //
#include <DataOps>
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

class KickPropagator
{
public:
  KickPropagator(Rotor &rotorref);
  ~KickPropagator();

  //inline bool apply(double *t,gsl_vector_complex *yinPtr)
  inline bool apply(double &t,std::vector< std::complex<double> > &yin)
  {
    //    clog << "kicker.apply(t,dt,yPTr) is in question" << endl;
    bool wasapplied = false;
    gsl_blas_ztrmv(UPPER,NOTRANS,CblasNonUnit,UmatPtr,yinPtr);
    *t += kickstepsize;
    //    clog << "\t\t\t... kicker.apply(t,dt,yPTr) is OK" << endl;
    return wasapplied;
  }
  
private:
  bool build();
   Rotor &rotorRef;
  const size_t dim;
  void * voidPtr;
  gsl_matrix_complex *UmatPtr;
  gsl_vector_complex *yPtr;
  std::vector<std::complex<double>> y;
  
  gsl_vector_complex_view Ucol;
  gsl_vector_complex_view Upart;
  
  double kickstepsize;

  void printUmat();
  void printcornerUmat();
  void printcornerUmat_real();
  void printcornerUmat_imag();
  void sampleyPtr();
};

class FreePropagator
{
public:
  FreePropagator(Rotor &rotorref, const double dtinau);
  ~FreePropagator();

  inline void apply(double *t, gsl_vector_complex *y) 
  {
    //    clog << "fixed freestep.apply(t,yPtr) in question ... " << endl;
    for (size_t i = 0; i < stepvec->size; i++){
      gsl_vector_complex_set(y,
			     i,
			     gsl_complex_mul(gsl_vector_complex_get(y,i),
					     gsl_complex_polar(1.0,fmod(gsl_vector_get(stepvec,i),2*M_PI) ) 
					     )
			     );
      // add phase = -(Ej-Ejstart)dt to coeffs
    }
    *t += dt;
    //    clog << "\t\t\t... fixed freestep.apply(t,yPtr) is OK. " << endl;
  }

  void apply(double *t, const double delta_tin, gsl_vector_complex *y);  
  
private:
  void build();
  
  Rotor &rotorRef;
  const unsigned dim;
  gsl_vector *stepvec ;
  gsl_vector *variablestep;
  double dt;
};


#endif
