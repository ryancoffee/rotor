

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
	t += delta_tin;
	return t;
}

// private:
bool jFreePropagator::build(jEnsemble & jens,const double & dt)
{
	Uvec.resize(jens.ej.size()-jens.getm());
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
// --- Onlu ever for signle m state hard baked into dimensions --- //


// --- This is the propagator that steps over the pulses.

jKickPropagator::jKickPropagator(jEnsemble & jens,const double & pulsetime)
: dim(jens.ej.size() - jens.getm())
, kickstepsize( pulsetime )
{
	Umat.resize(dim,dim);
}
jKickPropagator::jKickPropagator(jEnsemble & jens,PulseTime & pulse) 
: dim( jens.ej.size() - jens.getm())
, kickstepsize( pulse.duration() )
{
	Umat.resize(dim,dim);
}

jKickPropagator::~jKickPropagator()
{
}

bool jKickPropagator::build(jEnsemble & jens,PulseTime & pulse)
{
	std::clog << "\t\t\t ... Building kicker matrix ...\n\t\t\t\t\t.\n\t\t\t\t\t.\n\t\t\t\t\t.\n" << std::endl;

	double evolvetime = pulse.getduration();


	/*

	HERE HERE HERE HERE

	// allocate ode workspace
	const gsl_odeiv_step_type *T = gsl_odeiv_step_rkf45;  // for a smooth step // gsl_odeiv_step_rk8pd;  // for a fast step
	gsl_odeiv_step *s = gsl_odeiv_step_alloc(T,2*dim);
	gsl_odeiv_control *c = gsl_odeiv_control_y_new(ABSTOL,RELTOL);
	gsl_odeiv_evolve *e = gsl_odeiv_evolve_alloc(2*dim);
	//  gsl_odeiv_system sys = {kickfunc,kickjac,2*dim,voidPtr}; 
	//gsl_odeiv_system sys = {kickfunc,kickjac,2*dim,rotorRef.getvoidPtr()}; 
	*/

	cvec_t state(dim,0+0i);
	std::clog << "building m = " << jens.getm() << std::flush;
	for (size_t j=jens.getm(); j<dim; ++j){
		double t = 0;
		state.at(j) = 1. + 0i;
		std::clog << "\t j = " << j << std::flush;

		double t = rotorRef.gett0() - rotorRef.getCtau();
		double endtime = rotorRef.gett0() + rotorRef.getCtau();
		double h = 1/fsPau;                                   // initial step size
		int status = GSL_SUCCESS;
		//    std::clog << "evolveing system with yPtr->size = " << yPtr->size << std::endl;
		std::clog << t << "\n" << std::flush;
		//size_t steps = boost::numeric::odeint::integrate( ho , x , start_time , stop_time , step_time );
		size_t steps = boost::numeric::odeint::integrate( kicker , state , 0 , evolvetime , step_time );
		while (t<endtime && status==GSL_SUCCESS){
			//      std::clog << t << "\n" << std::flush;
			status = gsl_odeiv_evolve_apply(e,c,s,&sys,&t,endtime,&h,yPtr->data);

			Sqrnorm_Vector_Complex(yPtr);

			if (status !=GSL_SUCCESS){
				std::cerr << "evolve status: " << gsl_strerror(status) << std::endl;
				return false;
			}
		}
		std::clog << std::endl;

		for (int i = 0 ; i< dim ; i++){
			std::clog << std::setprecision(2) << yPtr->data[2*i] << "\t";
		}
		std::clog << std::endl;    

		// build kickpropagator matrix
		// --- return a vector view of the jth column of proagator matrix U(t,t0) and set to coeffs at t=endtime --- //

		Ucol = gsl_matrix_complex_column(UmatPtr,j);  

		/*
		   if (rotorRef.getm() == 0 && rotorRef.getv() == 0 && (j == 5 || j==6))
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


Kicker::Kicker(const jEnsemble &jens,const PulseTime &pulse)
{ 
}
/*
	int func (double t,const double y[],double f[],void *paraPtrvoid){
		PARAMS *paraPtr;
		paraPtr = (PARAMS *)paraPtrvoid;
		const int m=paraPtr->m;
		const int dim=paraPtr->dim;
		const int jstart = paraPtr->jstart;

		int realj;
		int *realjPtr=&realj;

		static double FF;
		FF = 0.0;

		static int i;

		for (i=0;i<dim;i++){
			f[i]=0.0;
		}

		if( inpulse(t,paraPtr,&FF) ){  // is coupling

			double aa,bb,cc;

			for (i=0;i<dim;i+=2){
				setrealj(realjPtr,&i,paraPtr);
				//      if (paraPtr->m==1 && paraPtr->jstart==0 && static_cast<int>(t) == 0){cout << *realjPtr << " ";}
				if (*realjPtr != -1){
					aa = (paraPtr->ejPtr[*realjPtr]-paraPtr->ejPtr[jstart+m]) - (1+paraPtr->aajmPtr[*realjPtr])*FF;
					bb = -1.0 * (paraPtr->bbjmPtr[*realjPtr])*FF;
					cc = -1.0 * (paraPtr->ccjmPtr[*realjPtr])*FF;
					f[i] = aa*y[i+1];
					f[i+1] = -1.0*aa*y[i];
					if (i>0) {
						f[i] += bb*y[i-2+1];
						f[i+1] -= bb*y[i-2];
					}
					if (i<dim-2) {
						f[i] += cc*y[i+3];
						f[i+1] -= cc*y[i+2];
					}
				}
			}

		} else {  // no coupling

			for (i=0;i<dim;i+=2){
				setrealj(realjPtr,&i,paraPtr);
				if (*realjPtr != -1){
					f[i]= (paraPtr->ejPtr[*realjPtr]-paraPtr->ejPtr[jstart+m])*y[i+1];
					f[i+1]= -(paraPtr->ejPtr[*realjPtr]-paraPtr->ejPtr[jstart+m])*y[i];
					// y' = E/i y; y = exp(-iEt)
				}
			}
		} 
		return GSL_SUCCESS;
	}

	int jac(double t,const double y[],double *dfdy,double dfdt[],void *paraPtrvoid){
		PARAMS *paraPtr;
		paraPtr = (PARAMS *)paraPtrvoid;
		const int m=paraPtr->m;
		const int dim = paraPtr->dim;
		const int jstart = paraPtr->jstart;

		static int realj;
		int *realjPtr=&realj;

		static int i;
		for (i=0;i<gsl_pow_2(dim);i++){
			dfdy[i]=0.0;
		}
		gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy,dim,dim);
		gsl_matrix *mat = &dfdy_mat.matrix;


		for (int i=0;i<dim;i++){
			dfdt[i]=0.0;
		}

		double FF=0.0, dFFdt = 0.0;

		if( inpulse(t,paraPtr,&FF,&dFFdt) ){  // is coupling


			gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy,dim,dim);
			gsl_matrix *mat = &dfdy_mat.matrix;


			for (i=0;i<dim;i+=2){
				setrealj(realjPtr,&i,paraPtr);
				if (*realjPtr != -1){
					gsl_matrix_set(mat,i,i+1,(paraPtr->ejPtr[*realjPtr]-paraPtr->ejPtr[jstart+m])-(1+paraPtr->aajmPtr[*realjPtr])*FF);
					gsl_matrix_set(mat,i,i-2+1,-(paraPtr->bbjmPtr[*realjPtr])*FF);
					gsl_matrix_set(mat,i,i+2+1,-(paraPtr->ccjmPtr[*realjPtr])*FF);
					gsl_matrix_set(mat,i+1,i,-((paraPtr->ejPtr[*realjPtr]-paraPtr->ejPtr[jstart+m])-(1+paraPtr->aajmPtr[*realjPtr])*FF));
					gsl_matrix_set(mat,i+1,i-2,(paraPtr->bbjmPtr[*realjPtr])*FF);
					gsl_matrix_set(mat,i+1,i+2,(paraPtr->ccjmPtr[*realjPtr])*FF);
				}
			}



			double dadt,dbdt,dcdt;


			for (i=0;i<dim;i+=2){
				setrealj(realjPtr,&i,paraPtr);
				if (*realjPtr != -1){
					dadt=-(1+paraPtr->aajmPtr[*realjPtr])*dFFdt;
					dbdt=-(paraPtr->bbjmPtr[*realjPtr])*dFFdt;
					dcdt=-(paraPtr->ccjmPtr[*realjPtr])*dFFdt;

					dfdt[i] = dadt*y[i+1];
					dfdt[i+1] = -1.0*dadt*y[i];
					if (i>0) {
						dfdt[i] += dbdt*y[i-1];
						dfdt[i+1] -= dbdt*y[i-2];
					}
					if (i<dim-2) {
						dfdt[i] += dcdt*y[i+3];
						dfdt[i+1] -= dcdt*y[i+2];
					}

				}
			}
		} else { // no coupling


			for (i=0;i<dim;i+=2){
				setrealj(realjPtr,&i,paraPtr);
				if (*realjPtr != -1){
					gsl_matrix_set(mat,i,i+1,paraPtr->ejPtr[*realjPtr]-paraPtr->ejPtr[jstart+m]);
					gsl_matrix_set(mat,i+1,i,-(paraPtr->ejPtr[*realjPtr]-paraPtr->ejPtr[jstart+m]));

				}
			}

		} 
		return GSL_SUCCESS;
	*/
