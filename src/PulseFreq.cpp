// Pulse class implimentation
// standard includes
#include <math.h>


// my headers
#include "PulseFreq.hpp"

/*
PulseFreq & PulseFreq::operator=(const PulseFreq &rhs)
	this->omega_center(rhs.omega_center),
        this->omega_width(rhs.omega_width),
	this->omega_high(rhs.omega_high),
	this->omega_onwidth(rhs.omega_onwidth),
        this->domega(rhs.domega),
        this->intime(rhs.intime),
        this->infreq(rhs.infreq),
        this->omega(rhs.omega), // these are static vectors... i'm trying not to make copies just yet
        this->time(rhs.time), // these are static vectors... i'm trying not to make copies just yet
	this->parent(false),
	this->child(true),
	this->i_low( static_cast<unsigned>(static_cast<double>( atof( getenv("nu_low") ) )* 2.0*M_PI*fsPau/domega) ),
	this->i_high( static_cast<unsigned>(static_cast<double>( atof( getenv("nu_high") ) )* 2.0*M_PI*fsPau/domega) )
	{
	this->samples = rhs.sample
	this->startind = rhs.startind;this->stopind=rhs.stopind;this->onwidth=rhs.onwidth;this->offwidth=rhs.offwidth;
	this->tspan = rhs.tspan;
	this->lambda_center=rhs.lambda_center;this->lambda_width=rhs.lambda_width;
	this->phase_GDD=rhs.phase_GDD;this->phase_TOD=rhs.phase_TOD;this->phase_4th=rhs.phase_4th;this->phase_5th=rhs.phase_5th;

        this->dtime = rhs.dtime;this->time_center=rhs.time_center;this->time_wdith=rhs.time_wdith;

	this->cwave=rhs.cwave;this->cspace=rhs.cspace;
	this->cvec = gsl_vector_complex_alloc(samples);
	this->rhovec = gsl_vector_alloc(samples);
	this->phivec = gsl_vector_alloc(samples);
        this->nu0=rhs.nu0;
		
	// point tables to parent //
	// copy vectors explicitly //
	gsl_vector_memcpy(this->rhovec,rhs.rhovec);
	gsl_vector_memcpy(this->phivec,rhs.phivec);
	gsl_vector_complex_memcpy(this->cvec,rhs.cvec);

        return *this;
}
*/


PulseFreq & PulseFreq::operator+=(const PulseFreq &rhs){
	gsl_vector_complex_add(cvec,rhs.cvec);
	for (unsigned i=0;i<samples;i++){
		gsl_vector_set(rhovec,i,gsl_complex_abs(gsl_vector_complex_get(cvec,i)));
		gsl_vector_set(phivec,i,gsl_complex_arg(gsl_vector_complex_get(cvec,i)));
	}
	return *this;
}

PulseFreq & PulseFreq::operator-=(const PulseFreq &rhs){
        gsl_vector_complex_sub(cvec,rhs.cvec);
        for (unsigned i=0;i<samples;i++){
                gsl_vector_set(rhovec,i,gsl_complex_abs(gsl_vector_complex_get(cvec,i)));
                gsl_vector_set(phivec,i,gsl_complex_arg(gsl_vector_complex_get(cvec,i)));
        }
        return *this;
}

PulseFreq & PulseFreq::operator*=(const PulseFreq &rhs){
	gsl_vector_complex_mul(cvec,rhs.cvec);
	for (unsigned i=0;i<samples;i++){
		gsl_vector_set(rhovec,i,gsl_complex_abs(gsl_vector_complex_get(cvec,i)));
		gsl_vector_set(phivec,i,gsl_complex_arg(gsl_vector_complex_get(cvec,i)));
	}
	return *this;
}

PulseFreq & PulseFreq::operator/=(const PulseFreq &rhs){
        gsl_vector_complex_div(cvec,rhs.cvec);
        for (unsigned i=0;i<samples;i++){
                gsl_vector_set(rhovec,i,gsl_complex_abs(gsl_vector_complex_get(cvec,i)));
                gsl_vector_set(phivec,i,gsl_complex_arg(gsl_vector_complex_get(cvec,i)));
        }
        return *this;
}

PulseFreq & PulseFreq::diffamps(const PulseFreq &rhs){
	gsl_complex z;
        gsl_vector_sub(rhovec,rhs.rhovec);
        for (unsigned i=0;i<samples;i++){
		z = gsl_complex_rect(gsl_vector_get(rhovec,i) , 0.0);
                gsl_vector_complex_set(cvec,i,z);
                gsl_vector_set_zero(phivec);
        }
        return *this;
}

PulseFreq & PulseFreq::normamps(const PulseFreq &rhs){
	gsl_complex z;
        gsl_vector_div(rhovec,rhs.rhovec);
        for (unsigned i=0;i<samples;i++){
		GSL_SET_COMPLEX(&z,gsl_vector_get(rhovec,i),0.0);
                gsl_vector_complex_set(cvec,i,z);
                gsl_vector_set_zero(phivec);
        }
        return *this;
}


void PulseFreq::attenuate(double attenfactor){
	gsl_complex z;
	for( unsigned i=0;i<samples;i++){
		rhovec->data[i] *= attenfactor;
		z = gsl_complex_polar(rhovec->data[i],phivec->data[i]);
		gsl_vector_complex_set(cvec,i,z);
	}
}
void PulseFreq::delay(double delayin){
	gsl_complex z;
	if(intime){
		fft_tofreq();
	}
	for (unsigned i=0;i<samples;i++){
		phivec->data[i] += omega[i]*delayin/fsPau;
		z = gsl_complex_polar(rhovec->data[i],phivec->data[i]);
		gsl_vector_complex_set(cvec,i,z);
	}
	if(intime){
		fft_totime();
	}
}

void PulseFreq::printfrequency(ofstream * outfile){
	double nu,lambda;
	(*outfile) << ("#omega[fs^-1]\tlambda[nm]\trho\tphi\n");
	for (unsigned i = i_low;i<i_high;i+=(unsigned)(atoi(getenv("sampleinterval")))){
		nu = (omega[i]/(2.0*M_PI)/fsPau);
		lambda = C_nmPfs/nu;
		(*outfile) << nu << "\t" << lambda << "\t" << rhovec->data[i] << "\t" << phivec->data[i] << "\n";
	}
}
void PulseFreq::appendfrequency(ofstream * outfile){
	double nu;
        for (unsigned i = i_low;i<i_high;i+=(unsigned)(atoi(getenv("sampleinterval")))){
       		(*outfile) << rhovec->data[i] << "\t" ;
        }
	(*outfile) << "\n";
}

void PulseFreq::appendnoisy(ofstream * outfile){
	time_t timeseed = std::time(NULL);
	double nu,outval;
	const gsl_rng_type * rngType = gsl_rng_taus;
	gsl_rng * rngPtr = gsl_rng_alloc(rngType);
	gsl_rng_set(rngPtr,timeseed);
	double noisescale = static_cast<double>( atof( getenv("noisescale") ) );
        for (unsigned i = i_low;i<i_high;i+=(unsigned)(atoi(getenv("sampleinterval")))){
		outval = rhovec->data[i];
		outval += gsl_ran_gaussian(rngPtr,noisescale);
       		(*outfile) << outval << "\t" ;
        }
	(*outfile) << "\n";
	gsl_rng_free(rngPtr);
}

void PulseFreq::appendfrequencybins(ofstream * outfile){
	double nu;
        for (unsigned i = i_low;i<i_high;i+=(unsigned)(atoi(getenv("sampleinterval")))){
		nu = (omega[i]/(2.0*M_PI)/fsPau);
       		(*outfile) << nu << "\t" ;
        }
	(*outfile) << "\n";
}
void PulseFreq::printfrequencydelay(ofstream * outfile, const double *delay){
	double nu,lambda,thisdelay;
	thisdelay = (*delay)*fsPau;
	(*outfile) << ("#omega[fs^-1]\tlambda[nm]\trho\tphi\tdelay[fs]\n");
	for (unsigned i = i_low;i<i_high;i+=(unsigned)(atoi(getenv("sampleinterval")))){
		nu = (omega[i]/(2.0*M_PI)/fsPau);
		lambda = C_nmPfs/nu;
		(*outfile) << nu << "\t" << lambda << "\t" << rhovec->data[i] << "\t" << phivec->data[i] << "\t" << thisdelay << "\n";
	}
	(*outfile) << "\n";
}
void PulseFreq::printfrequencydelaychirp(ofstream * outfile, const double *delay,const double *chirp){
	double nu,lambda,thisdelay,thischirp,reltime;
	thisdelay = (*delay)*fsPau;
	thischirp = (*chirp)*gsl_pow_2(fsPau);
	(*outfile) << ("#omega[fs^-1]\tlambda[nm]\trho\tphi\tdelay[fs]\tchirp[fs^2]\treldelays[fs]\n");
	for (unsigned i = i_low;i<i_high;i+=(unsigned)(atoi(getenv("sampleinterval")))){
		nu = (omega[i]/(2.0*M_PI)/fsPau);
		lambda = C_nmPfs/nu;
		reltime = (nu-nu0)*thischirp*10.0;
		(*outfile) << nu << "\t" << lambda << "\t" << rhovec->data[i] << "\t" << phivec->data[i] << "\t" << thisdelay << "\t" << thischirp << "\t" << reltime << "\n";
	}
	(*outfile) << "\n";
}

void PulseFreq::printwavelength(ofstream * outfile,const double *delay){
	double nu,lambda,lamlast;
	lamlast=0;
        for (unsigned i = i_high;i>i_low;i-=(unsigned)(atoi(getenv("sampleinterval")))){
                nu = (omega[i]/(2.0*M_PI)/fsPau);
                lambda = static_cast<double>(static_cast<int>((C_nmPfs/nu)*10.0))/10.0;
		if (lambda>lamlast & static_cast<int>(lambda*10)%5 == 0){
                	(*outfile) << lambda << "\t" << (*delay) << "\t" << rhovec->data[i] << "\n";
			lamlast = lambda;
		}
        }
	(*outfile) << "\n";
}

void PulseFreq::printtime(ofstream * outfile){
	(*outfile) << ("#time\treal\timag\n");
	/*
	for (unsigned i = 0;i<samples;i++){
		(*outfile) << (time[i]*fsPau) << "\t" << GSL_REAL(gsl_vector_complex_get(cvec,i)) << "\t" << GSL_IMAG(gsl_vector_complex_get(cvec,i)) << "\n";
	}
	*/
	for (unsigned i = samples/2;i<samples;i++){
		(*outfile) << (time[i]*fsPau) << "\t" << GSL_REAL(gsl_vector_complex_get(cvec,i)) << "\t" << GSL_IMAG(gsl_vector_complex_get(cvec,i)) << "\n";
	}
	for (unsigned i = 0;i<samples/2;i++){
		(*outfile) << (time[i]*fsPau) << "\t" << GSL_REAL(gsl_vector_complex_get(cvec,i)) << "\t" << GSL_IMAG(gsl_vector_complex_get(cvec,i)) << "\n";
	}

} 


void PulseFreq::allocatetables(void){
	cwave = gsl_fft_complex_wavetable_alloc(samples);
	while ( (cwave->factor[ (cwave->nf) -1 ] > cwave->factor[ (cwave->nf)-2 ]) ) {
		gsl_fft_complex_wavetable_free(cwave);
		samples += SAMPLEROUND;
		cwave = gsl_fft_complex_wavetable_alloc(samples);
	}
	cspace = gsl_fft_complex_workspace_alloc(samples);
}
void PulseFreq::buildvectors(void){
	static gsl_complex z;
	//factorization();

	cvec = gsl_vector_complex_alloc(samples);
	rhovec = gsl_vector_alloc(samples);
	phivec = gsl_vector_calloc(samples); // calloc since we will start with flat zero phase for a pulse centered about vector index 0 and samples
	modamp = gsl_vector_alloc(samples);
	modphase = gsl_vector_calloc(samples);
	omega = new double[samples];
	time = new double[samples];
	omega[0] = 0.0;
	time[0] = 0.0;
	rhovec->data[0] = 0.0;
	z = gsl_complex_polar(rhovec->data[0],phivec->data[0]);
	gsl_vector_complex_set(cvec,0,z);
	omega[samples/2] = -static_cast<double>(samples/2)*domega;
	time[samples/2] = -static_cast<double>(samples/2)*dtime;
	rhovec->data[samples/2] = 0.0;
	z = gsl_complex_polar(rhovec->data[samples/2],phivec->data[samples/2]);
	gsl_vector_complex_set(cvec,samples/2,z);
	startind = static_cast<unsigned>((omega_center-(omega_width/2.0))/domega);
	stopind = static_cast<unsigned>((omega_center+(omega_width/2.0))/domega);
	onwidth = static_cast<unsigned>(omega_onwidth/domega); // 2.0)/domega);//  /10.0)/domega);// sin^2 goes from0..1 in 0..pi/2
	offwidth = static_cast<unsigned>(omega_offwidth/domega); // 2.0)/domega);//  /10.0)/domega);// sin^2 goes from0..1 in 0..pi/2
	for (unsigned i = 1; i<startind;i++){
		omega[i] = domega*static_cast<double>(i);
		time[i] = dtime*static_cast<double>(i);
		rhovec->data[i] = 0.0;
		z = gsl_complex_polar(rhovec->data[i],phivec->data[i]);
		gsl_vector_complex_set(cvec,i,z);
		omega[samples-i] = -domega*static_cast<double>(i);
		time[samples-i] = -dtime*static_cast<double>(i);
		rhovec->data[samples-i] = 0.0;
		z = gsl_complex_polar(rhovec->data[samples-i],phivec->data[samples-i]);
		gsl_vector_complex_set(cvec,samples-i,z);
	}
	for (unsigned i = startind;i<startind+onwidth; i++){
		omega[i] = domega*static_cast<double>(i);
		time[i] = dtime*static_cast<double>(i);
		rhovec->data[i] = rising(i);
		z = gsl_complex_polar(rhovec->data[i],phivec->data[i]);
		gsl_vector_complex_set(cvec,i,z);
		omega[samples-i] = -domega*static_cast<double>(i);
		time[samples-i] = -dtime*static_cast<double>(i);
		rhovec->data[samples-i] = rising(i);
		z = gsl_complex_polar(rhovec->data[samples-i],phivec->data[samples-i]);
		gsl_vector_complex_set(cvec,samples-i,z);
	}
	for (unsigned i = startind+onwidth;i<stopind-offwidth; i++){
		omega[i] = domega*static_cast<double>(i);
		time[i] = dtime*static_cast<double>(i);
		rhovec->data[i] = 1.0;
		z = gsl_complex_polar(rhovec->data[i],phivec->data[i]);
		gsl_vector_complex_set(cvec,i,z);
		omega[samples-i] = -domega*static_cast<double>(i);
		time[samples-i] = -dtime*static_cast<double>(i);
		rhovec->data[samples-i] = 1.0;
		z = gsl_complex_polar(rhovec->data[samples-i],phivec->data[samples-i]);
		gsl_vector_complex_set(cvec,samples-i,z);
	}
	for (unsigned i = stopind-offwidth;i<stopind; i++){
		omega[i] = domega*static_cast<double>(i);
		time[i] = dtime*static_cast<double>(i);
		rhovec->data[i] = falling(i);
		z = gsl_complex_polar(rhovec->data[i],phivec->data[i]);
		gsl_vector_complex_set(cvec,i,z);
		omega[samples-i] = -domega*static_cast<double>(i);
		time[samples-i] = -dtime*static_cast<double>(i);
		rhovec->data[samples-i] = falling(i);
		z = gsl_complex_polar(rhovec->data[samples-i],phivec->data[samples-i]);
		gsl_vector_complex_set(cvec,samples-i,z);
	}
	for (unsigned i = stopind;i<samples/2; i++){
		omega[i] = domega*static_cast<double>(i);
		time[i] = dtime*static_cast<double>(i);
		rhovec->data[i] = 0.0;
		z = gsl_complex_polar(rhovec->data[i],phivec->data[i]);
		gsl_vector_complex_set(cvec,i,z);
		omega[samples-i] = -domega*static_cast<double>(i);
		time[samples-i] = -dtime*static_cast<double>(i);
		rhovec->data[samples-i] = 0.0;
		z = gsl_complex_polar(rhovec->data[samples-i],phivec->data[samples-i]);
		gsl_vector_complex_set(cvec,samples-i,z);
	}
	
}
void PulseFreq::killvectors(void){
	gsl_vector_free(modamp);
	gsl_vector_free(modphase);
	gsl_vector_free(rhovec);
	gsl_vector_free(phivec);
	gsl_vector_complex_free(cvec);
	gsl_fft_complex_wavetable_free(cwave);
	gsl_fft_complex_workspace_free(cspace);
	delete omega;
	delete time;
	omega = time = NULL;
}
void PulseFreq::killtheseonly(void){
	gsl_vector_free(modamp);
	gsl_vector_free(modphase);
	gsl_vector_free(rhovec);
	gsl_vector_free(phivec);
	gsl_vector_complex_free(cvec);
	gsl_fft_complex_workspace_free(cspace);
}
void PulseFreq::factorization(void){
        /* inspecting factorization of wavetable */
        clog << samples << " samples in " << cwave->nf << " factors:\t";
        for (unsigned i = 0;i<(cwave->nf);i++){
                clog << (cwave->factor)[i] << "\t";
        }
        clog << "\n";

}
