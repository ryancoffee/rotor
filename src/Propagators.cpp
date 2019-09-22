

// --- these classes are for the various coupled and free propagators for the rotor system --- //

// standard includes
#include <cmath>

// my headers
#include <Propagators.hpp>





// --- jFreePropagator --- //



jFreePropagator::jFreePropagator(jEnsemble & jens,const double dtinau = 4.)
: m_dt(dtinau)
{
	build(jens);
}

jFreePropagator::~jFreePropagator(){}

inline double jFreePropagator::apply(jEnsemble & jens,cvec_t & yin, double &t)
{
	return apply(jens,yin,t,m_dt);
}
inline double jFreePropagator::apply(jEnsemble & jens,cvec_t & yin, double &t, const double & delta_tin)
{
	std::transform(jens.ej.begin()+jens.m,jens.ej.end(),yin.begin(),yin.begin(),
		[&delta_tin](double & z, std::complex<double> & y){return y * std::polar(1.0,-delta_tin*z);}
	);
	return t;
}

// private:
bool jFreePropagator::build(jEnsemble & jens,const double & dt)
{
	Uvec.resize(jens.ej.size()-jens.m);
	std::transform(jens.ej.begin()+jens.m,jens.ej.end(),Uvec.begin(),
			[&dt](double &z) -> std::complex<double> { return std::polar(1.0,-dt*z) ;}
		      );
	return true;
}
bool jFreePropagator::build(jEnsemble & jens)
{
	return build(jens,m_dt);
}






// --- jKickPropagator --- //



// --- This is the propagator that steps over the pulses.

jKickPropagator::jKickPropagator(jEnsemble & jens,const double & pulsetime)
: dim(jens.ej.size())
, kickstepsize( pulsetime )
{
	Umat.resize(dim,dim);
	m_state.resize(dim,0.);
	m_stateHistory.resize(0);
	m_timesHistory.resize(0);
}
jKickPropagator::jKickPropagator(jEnsemble & jens,PulseTime & pulse) 
: dim( jens.ej.size())
, kickstepsize( pulse.duration() )
{
	Umat.resize(dim,dim);
	m_state.resize(dim,0.);
	m_stateHistory.resize(0);
	m_timesHistory.resize(0);
}

jKickPropagator::~jKickPropagator()
{
}

bool jKickPropagator::build(jEnsemble & jens,PulseTime & pulse)
{
	clog << "\t\t\t ... Building kicker matrix ...\n\t\t\t\t\t.\n\t\t\t\t\t.\n\t\t\t\t\t.\n" << endl;

	HERE HERE HERE HERE

	// allocate ode workspace
	const gsl_odeiv_step_type *T = gsl_odeiv_step_rkf45;  // for a smooth step // gsl_odeiv_step_rk8pd;  // for a fast step
	gsl_odeiv_step *s = gsl_odeiv_step_alloc(T,2*dim);
	gsl_odeiv_control *c = gsl_odeiv_control_y_new(ABSTOL,RELTOL);
	gsl_odeiv_evolve *e = gsl_odeiv_evolve_alloc(2*dim);
	//  gsl_odeiv_system sys = {kickfunc,kickjac,2*dim,voidPtr}; 
	//gsl_odeiv_system sys = {kickfunc,kickjac,2*dim,rotorRef.getvoidPtr()}; 

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
	return true;
}

