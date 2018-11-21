

// --- these classes are for the various coupled and free propagators for the rotor system --- //

// standard includes
#include <cmath>

// my headers
#include <Propagators.hpp>


// --- This is the propagator that steps over the pulses.

jKickPropagator::jKickPropagator(jEnsemble & jens) 
: rotorRef(rotorref)
, dim( jens.getsizej())
	UmatPtr(gsl_matrix_complex_calloc(dim,dim)),
	yPtr(gsl_vector_complex_alloc(dim))
{
	Umat = HERE HERE HERE HERE;
	voidPtr = (void * const)&rotorref;
	build();
}

jKickPropagator::~jKickPropagator()
{
}

double & jKickPropagator::apply(double &t, std::vector< std::complex<double> > &y)
{
	// add phase = -(Ej-Ejstart)dt to coeffs
	y = std::exp();
	std::transform(y.begin,y.end,jens.ej.begin + jstart);
	boost::numeric::ublas::prod(y);
	for (size_t i = 0; i < stepvec->size; i++){
		gsl_vector_complex_set(y,
				i,
				gsl_complex_mul(gsl_vector_complex_get(y,i),
					gsl_complex_polar(1.0,fmod(gsl_vector_get(stepvec,i),2*M_PI) ) 
					)
				);
	}
	t += dt;
	return t;
}


bool jKickPropagator::build(jEnsemble & jens,PulseTime & pulse)
{
	//  clog << "\t\t\t ... Building kicker matrix ...\n\t\t\t\t\t.\n\t\t\t\t\t.\n\t\t\t\t\t.\n" << endl;
	kickstepsize = 2 * pulse.getCtau();

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

void jKickPropagator::printUmat()
{
	Ucol = gsl_matrix_complex_column(UmatPtr,0);  
	for (unsigned l=0;l<Ucol.vector.size;l++){
		for (unsigned k = 0; k<Ucol.vector.size;k++){
			clog <<  setprecision(3) << gsl_complex_abs(gsl_matrix_complex_get(UmatPtr,l,k)) << "\t";
		}
		clog << endl;
	}

}
void jKickPropagator::printcornerUmat()
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
void jKickPropagator::printcornerUmat_real()
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
void jKickPropagator::printcornerUmat_imag()
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

void jKickPropagator::sampleyPtr(){
	for (unsigned k = 0; k < yPtr->size; k++){
		clog << setprecision(3) << gsl_complex_abs(gsl_vector_complex_get(yPtr,k)) << "\n";
	}
	clog << endl;
}

// --- FreePropagator --- //



jFreePropagator::jFreePropagator(const double dtinau = 4.)
: m_dt(dtinau)
{
	build();
}

jFreePropagator::~jFreePropagator(){}

inline double jFreePropagator::apply(jEnsemble & jens,std::vector<std::complex> > & yin, double &t)
{
	return apply(jens,yin,t,m_dt);
}
inline double jFreePropagator::apply(jEnsemble & jens,std::vector<std::complex> > & yin, double &t, const double & delta_tin)
{
	std::transform(jens.begin()+jens.m,jens.end(),y.begin(),y.begin(),
		[delta_tin](double & z, std::complex & y){return y * std::polar(1.0,-dt*z)}
	);
	return t;
}

// private:
void jFreePropagator::build(jEnsemble & jens,const double & dt)
{
	Uvec.resize(jens.ej.size()-jens.m);
	std::transform(jens.begin()+jens.m,jens.end(),Uvec.begin()
			[dt](double &z) -> std::complex<double> { return std::polar(1.0,-dt*z) ;}
		      );
	return true;
}
bool jFreePropagator::build(jEnsemble & jens)
{
	return build(jens,m_dt);
}
