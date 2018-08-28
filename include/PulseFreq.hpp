// Pulse class definition

#ifndef PULSEFREQ_H
#define PULSEFREQ_H


// standard includes
#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>

using namespace std;

// gsl includes
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// my headers
#include "Constants.hpp"

// my definitions
#define SAMPLEROUND 1000

#define NU_LOW 0.3
#define NU_HIGH 0.7
#define LAMBDA_LOW 400
#define LAMBDA_HIGH 900

class PulseFreq {

public:
PulseFreq & operator=(const PulseFreq &rhs); // function definitions must follow the PulseFreq definition
PulseFreq & operator+=(const PulseFreq &rhs); // function definitions must follow the PulseFreq definition
PulseFreq & operator-=(const PulseFreq &rhs); // function definitions must follow the PulseFreq definition
PulseFreq & operator*=(const PulseFreq &rhs); // function definitions must follow the PulseFreq definition
PulseFreq & operator/=(const PulseFreq &rhs); // function definitions must follow the PulseFreq definition
PulseFreq & diffamps(const PulseFreq &rhs);
PulseFreq & normamps(const PulseFreq &rhs);

public:
	PulseFreq(const double omcenter_in=(0.55*fsPau),const double omwidth_in=(0.15*fsPau),const double omonoff_in=(0.1*fsPau), double tspan_in=(10000.0/fsPau)): // default constructor
		omega_center(omcenter_in),
		omega_width(omwidth_in ),
		omega_high( GSL_MAX_DBL(4.0*(omcenter_in + omwidth_in),10.0*omcenter_in) ),
		domega( 2.0*M_PI/tspan_in),
		omega(NULL),
		intime(false),
		infreq(true),
		time(NULL),
		parent(true),
		child(false),
		i_low( static_cast<unsigned>(static_cast<double>( atof( getenv("nu_low") ) )* 2.0*M_PI*fsPau/domega) ),
		i_high( static_cast<unsigned>(static_cast<double>( atof( getenv("nu_high") ) )* 2.0*M_PI*fsPau/domega) )
	{
		samples = (( static_cast<unsigned>(2.0 * omega_high / domega))/SAMPLEROUND + 1 ) *SAMPLEROUND;// dt ~ .1fs, Dt ~ 10000fs, dom = 2*pi/1e4, omhigh = 2*pi/.1fs, samples = omhigh/dom*2
		allocatetables();
		dtime = tspan_in/static_cast<double>(samples);
		omega_onwidth = omega_offwidth = omega_width/2.0; // forcing sin2 gaussian spectrum
		buildvectors();
		nu0=omcenter_in/(2.0*M_PI)*fsPau;
		phase_GDD=phase_TOD=phase_4th=phase_5th=0.0;

	}
	PulseFreq(PulseFreq &rhs): // copy constructor
		omega_center(rhs.omega_center),
                omega_width(rhs.omega_width),
                omega_high(rhs.omega_high),
		omega_onwidth(rhs.omega_onwidth),
                domega(rhs.domega),
                intime(rhs.intime),
                infreq(rhs.infreq),
                omega(rhs.omega), // these are static vectors... i'm trying not to make copies just yet
                time(rhs.time), // these are static vectors... i'm trying not to make copies just yet
		parent(false),
		child(true),
		i_low( static_cast<unsigned>(static_cast<double>( atof( getenv("nu_low") ) )* 2.0*M_PI*fsPau/domega) ),
		i_high( static_cast<unsigned>(static_cast<double>( atof( getenv("nu_high") ) )* 2.0*M_PI*fsPau/domega) )
	{
		samples = rhs.samples;
		startind = rhs.startind;stopind=rhs.stopind;onwidth=rhs.onwidth;offwidth=rhs.offwidth;
		tspan = rhs.tspan;
		lambda_center=rhs.lambda_center;lambda_width=rhs.lambda_width;
		phase_GDD=rhs.phase_GDD;phase_TOD=rhs.phase_TOD;phase_4th=rhs.phase_4th;phase_5th=rhs.phase_5th;

                dtime = rhs.dtime;time_center=rhs.time_center;time_wdith=rhs.time_wdith;

		cwave=rhs.cwave;
		cspace = gsl_fft_complex_workspace_alloc(samples);
		cvec = gsl_vector_complex_alloc(samples);
		rhovec = gsl_vector_alloc(samples);
		phivec = gsl_vector_alloc(samples);
		modamp = gsl_vector_alloc(samples);
                modphase = gsl_vector_alloc(samples);

                nu0=rhs.nu0;
		
		// point tables to parent //
		// copy vectors explicitly //
		gsl_vector_memcpy(rhovec,rhs.rhovec);
		gsl_vector_memcpy(phivec,rhs.phivec);
		gsl_vector_complex_memcpy(cvec,rhs.cvec);
		gsl_vector_memcpy(modamp,rhs.modamp);
		gsl_vector_memcpy(modphase,rhs.modphase);

	}

	~PulseFreq(void){
		if(child){
			killtheseonly();
		} else {
			killvectors();
		}
	}
	inline unsigned getsamples(void) {
		return samples;
	}
	inline unsigned getdt(void) {
		return dtime;
	}
	inline void fft_totime(void) {
		gsl_fft_complex_inverse(cvec->data,cvec->stride,cvec->size,cwave,cspace);
		for (unsigned i = 0;i<(cvec->size);i++){
			gsl_vector_set(rhovec,i,gsl_complex_abs(gsl_vector_complex_get(cvec,i)));
			gsl_vector_set(phivec,i,gsl_complex_arg(gsl_vector_complex_get(cvec,i)));
		}
		infreq = false;
		intime = true;
	}
	inline void fft_tofreq(void) {
		if (infreq)
			cerr << "died here at fft_tofreq()" << endl;
		gsl_fft_complex_forward(cvec->data,cvec->stride,cvec->size,cwave,cspace);
		for (unsigned i = 0;i<(cvec->size);i++){
			gsl_vector_set(rhovec,i,gsl_complex_abs(gsl_vector_complex_get(cvec,i)));
			gsl_vector_set(phivec,i,gsl_complex_arg(gsl_vector_complex_get(cvec,i)));
		}
		infreq=true;
		intime=false;
	}
	
	inline int addchirp(double chirp_in) {
		if (intime){
			cerr << "whoops, trying to add phase in the time domain" << endl;
			return 1;
		}
		phase_GDD = chirp_in;
		addGDDtoindex(0,1);
		addGDDtoindex(samples/2,-1);
		for (unsigned i = 1; i<samples/2;i++){
			addGDDtoindex(i,1);
			addGDDtoindex(samples-i,-1);
		}
		return 0;
	}

	inline int addchirp(double* chirp_in) {
		if (intime){
			cerr << "whoops, trying to add phase in the time domain" << endl;
			return 1;
		}
		phase_GDD = chirp_in[0];
		phase_TOD = chirp_in[1];
		phase_4th = chirp_in[2];
		phase_5th = chirp_in[3];
		addGDDtoindex(0,1);
		addGDDtoindex(samples/2,-1);
		addTODtoindex(0,1);
		addTODtoindex(samples/2,-1);
		add4thtoindex(0,1);
		add4thtoindex(samples/2,-1);
		add5thtoindex(0,1);
		add5thtoindex(samples/2,-1);
		for (unsigned i = 1; i<samples/2;i++){
			addGDDtoindex(i,1);
			addGDDtoindex(samples-i,-1);
			addTODtoindex(i,1);
			addTODtoindex(samples-i,-1);
			add4thtoindex(i,1);
			add4thtoindex(samples-i,-1);
			add5thtoindex(i,1);
			add5thtoindex(samples-i,-1);
		}
		return 0;
	}


	void attenuate(double attenfactor);
	void delay(double delayin);

	inline int modulateamp_time(const gsl_vector *modulation) {
		if (infreq){
			cerr << "whoops, trying time modulation but in frequency domain" << endl;
			return 1; }
		if (modulation->size < cvec->size){
			cerr << "size mismatch, out of range in modulateamp_time()" << endl;
			return 2;
		}
		modampatindx(0,modulation);
		modampatindx(samples/2,modulation);
		for (unsigned i = 1;i<samples/2;i++){
			modampatindx(i,modulation);
			modampatindx(samples-i,modulation);
		}
		return 0;
	}
        inline int modulateamp_time(void) {
                if (infreq){
                        cerr << "whoops, trying time modulation but in frequency domain" << endl;
                        return 1; }
                modampatindx(0);
                modampatindx(samples/2);
                for (unsigned i = 1;i<samples/2;i++){
                        modampatindx(i);
                        modampatindx(samples-i);
                }
                return 0;
        }

	inline int modulatephase_time(const gsl_vector *modulation) {
                if (infreq){
                        cerr << "whoops, trying time modulation but in frequency domain" << endl;
                        return 1;
                }
                if (modulation->size < cvec->size){
                        cerr << "size mismatch, out of range in modulateamp_time()" << endl;
                        return 2;
                }
                modphaseatindx(0,modulation);
                modphaseatindx(samples/2,modulation);
                for (unsigned i = 1;i<samples/2;i++){
                        modphaseatindx(i,modulation);
                        modphaseatindx(samples-i,modulation);
                }
                return 0;
	}
	inline int modulatephase_time(void) {
                if (infreq){
                        cerr << "whoops, trying time modulation but in frequency domain" << endl;
                        return 1;
                }
                modphaseatindx(0);
                modphaseatindx(samples/2);
                for (unsigned i = 1;i<samples/2;i++){
                        modphaseatindx(i);
                        modphaseatindx(samples-i);
                }
                return 0;
        }


	inline int modulateamp_freq(const gsl_vector *modulation){
		if (intime){
			cerr << "whoops, trying to freq modulate but in time domain" << endl;
			return 1;
		}
		if (modulation->size < cvec->size){
			cerr << "size mismatch, out of range in modulateamp_freq()" << endl;
			return 2;
		}

		return 0;
	}
	void appendfrequency(ofstream * outfile);
	void appendnoisy(ofstream * outfile);
	void appendfrequencybins(ofstream * outfile);
	void printfrequency(ofstream * outfile);
	void printfrequencydelay(ofstream * outfile, const double *delay);
	void printfrequencydelaychirp(ofstream * outfile, const double *delay,const double *chirp);
	void printtime(ofstream * outfile);
	void printwavelength(ofstream * outfile,const double *delay);
	inline double gettime(unsigned ind){return (time[ind]*fsPau);}

	gsl_vector * modamp;
	gsl_vector * modphase;

private:
	const bool parent,child;
	bool intime,infreq;
	unsigned samples;
	unsigned startind,stopind,onwidth,offwidth;
	double tspan;
	double domega,lambda_center,lambda_width,omega_center,omega_width,omega_high;
	double omega_onwidth,omega_offwidth;
	double phase_GDD,phase_TOD,phase_4th,phase_5th;
	double dtime,time_center,time_wdith;
	gsl_fft_complex_wavetable * cwave;
	gsl_fft_complex_workspace * cspace;
	gsl_vector_complex * cvec;
	gsl_vector * rhovec;
	gsl_vector * phivec;
	double * omega, *time;
	double nu0;

	const unsigned i_low, i_high;


	void allocatetables(void);
	void buildvectors(void);
	void factorization(void);
	void killvectors(void);
	void killtheseonly(void);

	inline void addGDDtoindex(const unsigned indx,const int omega_sign) {
		static gsl_complex z;
		phivec->data[indx] += omega_sign*phase_GDD*gsl_pow_2(omega[indx]-(static_cast<double>(omega_sign)*omega_center));
		z = gsl_complex_polar(rhovec->data[indx],phivec->data[indx]);
		gsl_vector_complex_set(cvec,indx,z);
	}	
	inline void addTODtoindex(const unsigned indx,const int omega_sign) {
		static gsl_complex z;
		phivec->data[indx] += omega_sign*phase_TOD*gsl_pow_3(omega[indx]-(static_cast<double>(omega_sign)*omega_center));
		z = gsl_complex_polar(rhovec->data[indx],phivec->data[indx]);
		gsl_vector_complex_set(cvec,indx,z);
	}	
	inline void add4thtoindex(const unsigned indx,const int omega_sign) {
		static gsl_complex z;
		phivec->data[indx] += omega_sign*phase_4th*gsl_pow_4(omega[indx]-(static_cast<double>(omega_sign)*omega_center));
		z = gsl_complex_polar(rhovec->data[indx],phivec->data[indx]);
		gsl_vector_complex_set(cvec,indx,z);
	}	
	inline void add5thtoindex(const unsigned indx,const int omega_sign) {
		static gsl_complex z;
		phivec->data[indx] += omega_sign*phase_5th*gsl_pow_5(omega[indx]-(static_cast<double>(omega_sign)*omega_center));
		z = gsl_complex_polar(rhovec->data[indx],phivec->data[indx]);
		gsl_vector_complex_set(cvec,indx,z);
	}	

	inline void modampatindx(const unsigned indx,const gsl_vector * modvec) {
		static gsl_complex z;
		rhovec->data[indx] *= gsl_vector_get(modvec,indx);
		z = gsl_complex_polar(rhovec->data[indx],phivec->data[indx]);
		gsl_vector_complex_set(cvec,indx,z);
	}
	inline void modampatindx(const unsigned indx) {
		static gsl_complex z;
		rhovec->data[indx] *= gsl_vector_get(modamp,indx);
		z = gsl_complex_polar(rhovec->data[indx],phivec->data[indx]);
		gsl_vector_complex_set(cvec,indx,z);
	}
/*
E(t) = E0exp(iwt)=E0 exp(i w(t) (t-t0))
		= E0 exp(i (w0+GDD*t) (t-t0) )
		= E0 exp(i ((w0 * t) - (GDD*t*t0) + (GDD*t**2) - (w0*t0)) so adding a delay means adding an arbitrary phase of w0*t0 and subtracting GDD*t*t0 ( a linear term in phase which we know only dials the pulse position in z, or likewise phase in time... unless t0 is actually t0(t), then you have to include this one as well as w0*t0(t).
		= E0 exp(i ((w0 * t) + (GDD*t**2) - (w0*t0) - (GDD*t*t0) ) )
		= E0 exp(i w(t)*t) exp(i -w(t)*t0(t))
		= E(t) exp(i -(w0 + 1/GDD*t + 1/TOD*t**2 + 1/FOD*t**3)*t0(t))
*/
	inline void modphaseatindx(const unsigned indx,const gsl_vector * modvec) {
		gsl_complex z;
		if (gsl_vector_get(modvec,indx)!=0){
			phivec->data[indx] += gsl_vector_get(modvec,indx);
			z = gsl_complex_polar(rhovec->data[indx],phivec->data[indx]);
			gsl_vector_complex_set(cvec,indx,z);
		}
	}
	inline void modphaseatindx(const unsigned indx) {
		gsl_complex z;
		double val;
		val = gsl_vector_get(modphase,indx);
		if (val != 0){
			phivec->data[indx] += val;
			z = gsl_complex_polar(rhovec->data[indx],phivec->data[indx]);
			gsl_vector_complex_set(cvec,indx,z);
		}
	}
	inline double rising(const unsigned indx) {
		return gsl_pow_2(sin(static_cast<double>(M_PI_2*(indx-startind)/onwidth)));
	}
	inline double falling(const unsigned indx) {
		return gsl_pow_2(sin(static_cast<double>(M_PI_2*(stopind-indx)/offwidth)));
	}
};

/*
PulseFreq & PulseFreq::operator=(const PulseFreq &rhs){
	if(this != &rhs){
		*this = rhs;
	}
	return *this;
}
*/

/*
PulseFreq & PulseFreq::operator*=(const double rhs){
	gsl_complex z;
	for (unsigned i=0;i<samples;i++){
		rhovec->data[i] *= rhs;
		z=gsl_complex_polar(rhovec->data[i],phivec->data[i]);
		gsl_vector_complex_set(cvec,i,z);
	}
	return *this;
}
*/

#endif
