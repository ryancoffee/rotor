

// --- these classes are for the various coupled and free propagators for the rotor system --- //

// standard includes
#include <math.h>
#include <iostream>
#include <iomanip>
using namespace std;
// gsl includes
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_odeiv.h>

// my headers
#include "Constants.h" // --- conversion constancts lie a0, Eh, icmPau, fsPau etc. --- //
#include "Pulse.hpp" // --- this handels the pulse envelope parameters like strength, duration and t0, and returns fales of FF and dFFdt --- //
#include "Temp.hpp" // --- this give temperature object, and can represent it in various units --- //
#include "LookUp.hpp" // --- this is the 3j symbol lookup table --- //
#include "Signal.hpp" // --- this handles the time and data vectors, updates and file writing --- //
#include "Rotor.hpp" // --- this is the molecule specifics like ej and pj vectors --- //
#include "Inlines.hpp"  // --- dimension and index conversion functions ( and other inline functions ) --- //
#include "Members.hpp" // --- This is the colleciton of member functions ( kickfunc and kickjac for now ) --- //
#include "Params.hpp" // --- This is the collection of pointers to objects that gets passed to kickfunc and kickjac mainly --- //
#include "Propagators.hpp"

// my preprocessor defs
#define ABSTOL 1e-1// 1e-10
#define RELTOL 1e-1//1e-8        // relative tolerance reltol = 1e-3 in matlab default

// gsl_blas defs
#define DAGGAR CblasConjTrans
#define TRANS CblasTrans
#define NOTRANS CblasNoTrans
#define RIGHT CblasRight
#define LEFT CblasLeft
#define UPPER CblasUpper
#define LOWER CblasLower

using namespace std;

// --- This is the propagatof that steps over the pulses.

KickPropagator::KickPropagator(Rotor &rotorref) :  rotorRef(rotorref),
						   dim( rotorRef.getsizej() - rotorRef.getm() ),
						   UmatPtr(gsl_matrix_complex_calloc(dim,dim)),
						   yPtr(gsl_vector_complex_alloc(dim))
{
  voidPtr = (void * const)&rotorref;
  build();
}

KickPropagator::~KickPropagator()
{
  gsl_matrix_complex_free(UmatPtr);
  gsl_vector_complex_free(yPtr);  
}


bool KickPropagator::build()
{
  //  clog << "\t\t\t ... Building kicker matrix ...\n\t\t\t\t\t.\n\t\t\t\t\t.\n\t\t\t\t\t.\n" << endl;
  kickstepsize = 2 * rotorRef.getCtau();
  
  // allocate ode workspace
  const gsl_odeiv_step_type *T = gsl_odeiv_step_rkf45;  // for a smooth step // gsl_odeiv_step_rk8pd;  // for a fast step
  gsl_odeiv_step *s = gsl_odeiv_step_alloc(T,2*dim);
  gsl_odeiv_control *c = gsl_odeiv_control_y_new(ABSTOL,RELTOL);
  gsl_odeiv_evolve *e = gsl_odeiv_evolve_alloc(2*dim);
  //  gsl_odeiv_system sys = {kickfunc,kickjac,2*dim,voidPtr}; 
  gsl_odeiv_system sys = {kickfunc,kickjac,2*dim,rotorRef.getvoidPtr()}; 

  for (size_t j=0; j<dim; j++){
    gsl_vector_complex_set_basis(yPtr,j);
    rotorRef.setjstart(j+rotorRef.getm());
    clog << "building m = " << rotorRef.getm() << "\t j = " << rotorRef.getjstart() << flush;

    double t = rotorRef.gett0() - rotorRef.getCtau();
    double endtime = rotorRef.gett0() + rotorRef.getCtau();
    double h = 1/fsPau;                                   // initial step size
    int status = GSL_SUCCESS;
    //    clog << "evolveing system with yPtr->size = " << yPtr->size << endl;
      clog << t << "\n" << flush;
    while (t<endtime && status==GSL_SUCCESS){
      //      clog << t << "\n" << flush;
      status = gsl_odeiv_evolve_apply(e,c,s,&sys,&t,endtime,&h,yPtr->data);
      Sqrnorm_Vector_Complex(yPtr);
      
      if (status !=GSL_SUCCESS){
	cerr << "evolve status: " << gsl_strerror(status) << endl;
	return false;
      }
    }
    clog << endl;

    for (int i = 0 ; i< dim ; i++){
      clog << setprecision(2) << yPtr->data[2*i] << "\t";
	}
    clog << endl;    

    // build kickpropagator matrix
    // --- return a vector view of the jth column of proagator matrix U(t,t0) and set to coeffs at t=endtime --- //

      Ucol = gsl_matrix_complex_column(UmatPtr,j);  

      /*
      if (rotorRef.getm() == 0 && rotorRef.getv() == 0 && (j == 5 || j==6))
	sampleyPtr();
      */
      for (size_t i=0; i < Ucol.vector.size; i++){
	gsl_vector_complex_set( &Ucol.vector,i,
				gsl_vector_complex_get(yPtr,i));
      }
      
  }
  
  // free the ode workspace 
  
  gsl_odeiv_evolve_free(e);
  gsl_odeiv_control_free(c);
  gsl_odeiv_step_free(s);
  if (rotorRef.getm() == 0 && rotorRef.getv() == 0)
    {
      printcornerUmat_real();
      printcornerUmat_imag();
    }
  return true;
}

void KickPropagator::printUmat()
{
  Ucol = gsl_matrix_complex_column(UmatPtr,0);  
  for (unsigned l=0;l<Ucol.vector.size;l++){
    for (unsigned k = 0; k<Ucol.vector.size;k++){
      clog <<  setprecision(3) << gsl_complex_abs(gsl_matrix_complex_get(UmatPtr,l,k)) << "\t";
    }
    clog << endl;
  }
  
}
void KickPropagator::printcornerUmat()
{
  Ucol = gsl_matrix_complex_column(UmatPtr,0);  
  unsigned quit =  GSL_MIN(Ucol.vector.size,10);
  for (unsigned l=0; l < quit;l++){
    for (unsigned k = 0; k<quit;k++){
      clog <<  setprecision(3) << gsl_complex_abs(gsl_matrix_complex_get(UmatPtr,l,k)) << "\t";
    }
    clog << endl;
  }
  
}
void KickPropagator::printcornerUmat_real()
{
  clog << "\t\t\t --- Corner of Umat_real ---\n";
  Ucol = gsl_matrix_complex_column(UmatPtr,0);  
  unsigned quit =  GSL_MIN(Ucol.vector.size,10);
  for (unsigned l=0; l < quit;l++){
    for (unsigned k = 0; k<quit;k++){
      clog <<  setprecision(3) << GSL_REAL(gsl_matrix_complex_get(UmatPtr,l,k)) << "\t";
    }
    clog << endl;
  }
  
}
void KickPropagator::printcornerUmat_imag()
{
  clog << "\t\t\t --- Corner of Umat_imag ---\n";
  Ucol = gsl_matrix_complex_column(UmatPtr,0);  
  unsigned quit =  GSL_MIN(Ucol.vector.size,10);
  for (unsigned l=0; l < quit;l++){
    for (unsigned k = 0; k<quit;k++){
      clog <<  setprecision(3) << GSL_IMAG(gsl_matrix_complex_get(UmatPtr,l,k)) << "\t";
    }
    clog << endl;
  }
  
}

void KickPropagator::sampleyPtr(){
  for (unsigned k = 0; k < yPtr->size; k++){
    clog << setprecision(3) << gsl_complex_abs(gsl_vector_complex_get(yPtr,k)) << "\n";
  }
  clog << endl;
}

// --- FreePropagator --- //



FreePropagator::FreePropagator(Rotor &rotorref, const double dtinau) :
  rotorRef(rotorref),
  dim(rotorRef.getsizej() - rotorRef.getm()),
  stepvec(gsl_vector_calloc(dim)),
  variablestep(gsl_vector_calloc(dim))

{
  dt = dtinau;
  build();
}

FreePropagator::~FreePropagator(){ gsl_vector_free(stepvec); gsl_vector_free(variablestep);}

void FreePropagator::apply(double *t, const double delta_tin, gsl_vector_complex *y)
{
  //  clog << "variable freestep.apply(t,dt,yPTr) is in question ..." << endl;
  bool success = false;
  success = rotorRef.ejinau(variablestep);

  if (!success){
    /*
    cerr << "Something went wrong in FreePropagator.apply(t,dt,y)" 
	 << "\nvariablestep->size = " << variablestep->size 
	 << "\ny->size = " << y->size 
	 << endl;
    */
    return;
  }
  gsl_vector_scale(stepvec,-delta_tin);    
  for (size_t i = 0; i < variablestep->size; i++){
    gsl_vector_complex_set(y,
			   i,
			   gsl_complex_mul(gsl_vector_complex_get(y,i),
					   gsl_complex_polar(1.0,fmod(gsl_vector_get(variablestep,i),2*M_PI) ) 
					   )
			   );
    // add phase = -(Ej-Ejstart)dt to coeffs
  }
  *t += delta_tin;
  //  clog << "\t\t\t... variable freestep.apply(t,dt,yPTr) is OK" << endl;
}


// private:
void FreePropagator::build()
{
  //  clog << "\t\t\t ... Building freestep vector ...\n\t\t\t\t\t.\n\t\t\t\t\t.\n\t\t\t\t\t.\n" << endl; 
  rotorRef.ejinau(stepvec);
  gsl_vector_scale(stepvec,-dt);
}
