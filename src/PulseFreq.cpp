// Pulse class implimentation
// standard includes
#include <cmath>
#include <iterator>

// my headers
#include <PulseFreq.hpp>
#include <boost/math/interpolators/barycentric_rational.hpp>
#include <fftw3.h>
#include <algorithm>
#include <Constants.hpp>
#include <boost/lexical_cast.hpp>
#include <DataOps.hpp>
#include <random>
#include <cassert>

using namespace Constants;
using namespace DataOps;

PulseFreq::PulseFreq(const double omcenter_in=(0.55*fsPau<double>()),const double omwidth_in=(0.15*fsPau<double>()),const double omonoff_in=(0.1*fsPau<double>()), double tspan_in=(10000.0/fsPau<double>())):
	omega_center(omcenter_in),
	omega_width(omwidth_in ),
	omega_high( std::max(4.0*(omcenter_in + omwidth_in),10.0*omcenter_in) ),
	domega( 2.0*pi<double>()/tspan_in),
	intime(false),
	infreq(true),
	m_noisescale(1e-3),
	m_sampleinterval(2),
	m_saturate(4096),
	m_gain(1000000),
	m_lamsamples(1024),
	sampleround(1000),
	cvec(NULL),
	r_vec(NULL),
	hc_vecFT(NULL),
	r_vec_2x(NULL),
	hc_vec_2xFT(NULL)
{
	std::cerr << "In constructor PulseFreq()" << std::endl;
	i_low =  (unsigned)(double( atof( getenv("nu_low") ) )* twopi<double>()*fsPau<double>()/domega);
	i_high =  (unsigned)(double( atof( getenv("nu_high") ) )* twopi<double>()*fsPau<double>()/domega);
	samples = (( (unsigned)(2.0 * omega_high / domega))/sampleround + 1 ) *sampleround;// dt ~ .1fs, Dt ~ 10000fs, dom = 2*pi/1e4, omhigh = 2*pi/.1fs, samples = omhigh/dom*2
	dtime = tspan_in/double(samples);
	omega_onwidth = omega_offwidth = omega_width/2.0; // forcing sin2 gaussian spectrum
	buildvectors();
	nu0=omcenter_in/(2.0*pi<double>())*fsPau<double>();
	phase_GDD=phase_TOD=phase_4th=phase_5th=0.0;
	m_lamsamples = (size_t)atoi(getenv("lamsamples"));
	m_gain = boost::lexical_cast<float>( getenv("gain"));
	m_noisescale = boost::lexical_cast<double>( getenv("noisescale") ) ;
	m_sampleinterval = boost::lexical_cast<size_t>(getenv("sampleinterval"));
	m_saturate = uint16_t( boost::lexical_cast<int>( getenv("saturate")));
}


PulseFreq::PulseFreq(PulseFreq &rhs): // deep-ish copy constructor
	omega_center(rhs.omega_center),
	omega_width(rhs.omega_width),
	omega_high(rhs.omega_high),
	omega_onwidth(rhs.omega_onwidth),
	domega(rhs.domega),
	intime(rhs.intime),
	infreq(rhs.infreq),
	i_low(rhs.i_low), 
	i_high(rhs.i_high), 
	m_noisescale(rhs.m_noisescale),
	m_sampleinterval(rhs.m_sampleinterval),
	m_saturate(rhs.m_saturate),
	m_gain(rhs.m_gain),
	m_lamsamples(rhs.m_lamsamples),
	sampleround(1000)
{
	//std::cerr << "\t\t\t+++++  Copy constructor of PulseFreq::PulseFreq(PulseFreq &rhs)\n" << std::flush;
	DataOps::clone(omega,rhs.omega);
	DataOps::clone(time,rhs.time);

	samples = rhs.samples;
	startind = rhs.startind;stopind=rhs.stopind;onwidth=rhs.onwidth;offwidth=rhs.offwidth;
	tspan = rhs.tspan;
	lambda_center=rhs.lambda_center;lambda_width=rhs.lambda_width;
	phase_GDD=rhs.phase_GDD;phase_TOD=rhs.phase_TOD;phase_4th=rhs.phase_4th;phase_5th=rhs.phase_5th;

	dtime = rhs.dtime;time_center=rhs.time_center;time_wdith=rhs.time_wdith;


	nu0=rhs.nu0;
	FTplan_forwardPtr = rhs.FTplan_forwardPtr; 
	FTplan_backwardPtr = rhs.FTplan_backwardPtr; 
	buildvectors();

	DataOps::clone(rhovec,rhs.rhovec);
	DataOps::clone(phivec,rhs.phivec);
	DataOps::clone(cvec,rhs.cvec,samples);
	DataOps::clone(r_vec,rhs.r_vec,samples);
	DataOps::clone(hc_vecFT,rhs.hc_vecFT,samples);
	DataOps::clone(r_vec_2x,rhs.r_vec_2x,2*samples);
	DataOps::clone(hc_vec_2xFT,rhs.hc_vec_2xFT,2*samples);

	DataOps::clone(modamp,rhs.modamp);
	DataOps::clone(modphase,rhs.modphase);
}

PulseFreq & PulseFreq::operator=(const PulseFreq & rhs) // shallow-ish assignment
{
	//std::cerr << "\t\t\t+++++  Shallow copy of PulseFreq::operator=\n" << std::flush;
	omega_center=rhs.omega_center;
	omega_width=rhs.omega_width;
	omega_high=rhs.omega_high;
	omega_onwidth=rhs.omega_onwidth;
	domega=rhs.domega;
	intime=rhs.intime;
	infreq=rhs.infreq;
	omega=rhs.omega; // these are static vectors... i'm trying not to make copies just yet
	time=rhs.time; // these are static vectors... i'm trying not to make copies just yet
	i_low=rhs.i_low; 
	i_high=rhs.i_high; 
	m_noisescale=rhs.m_noisescale;
	m_sampleinterval=rhs.m_sampleinterval;
	m_saturate=rhs.m_saturate;
	m_gain=rhs.m_gain;
	m_lamsamples=rhs.m_lamsamples;

	samples = rhs.samples;
	startind = rhs.startind;stopind=rhs.stopind;onwidth=rhs.onwidth;offwidth=rhs.offwidth;
	tspan = rhs.tspan;
	lambda_center=rhs.lambda_center;lambda_width=rhs.lambda_width;
	phase_GDD=rhs.phase_GDD;phase_TOD=rhs.phase_TOD;phase_4th=rhs.phase_4th;phase_5th=rhs.phase_5th;

	dtime = rhs.dtime;time_center=rhs.time_center;time_wdith=rhs.time_wdith;

	nu0=rhs.nu0;
	FTplan_forwardPtr = rhs.FTplan_forwardPtr; 
	FTplan_backwardPtr = rhs.FTplan_backwardPtr; 

	DataOps::clone(rhovec,rhs.rhovec);
	DataOps::clone(phivec,rhs.phivec);
	DataOps::clone(cvec,rhs.cvec,samples);
	DataOps::clone(r_vec,rhs.r_vec,samples);
	DataOps::clone(hc_vecFT,rhs.hc_vecFT,samples);
	DataOps::clone(r_vec_2x,rhs.r_vec_2x,2*samples);
	DataOps::clone(hc_vec_2xFT,rhs.hc_vec_2xFT,2*samples);

	DataOps::clone(modamp,rhs.modamp);
	DataOps::clone(modphase,rhs.modphase);
	return *this;
}

PulseFreq::~PulseFreq(void){
	killvectors();
}


void PulseFreq::rhophi2cvec(void)
{
	for (size_t i=0;i<samples;i++){
		cvec[i] = std::polar(rhovec[i],phivec[i]);
	}
}
void PulseFreq::cvec2rhophi(void)
{
	for (size_t i=0;i<samples;i++){
		rhovec[i] = std::abs(cvec[i]);
		phivec[i] = std::arg(cvec[i]);
	}
}

PulseFreq & PulseFreq::operator+=(const PulseFreq &rhs){
	DataOps::sum(cvec,rhs.cvec,samples);
	cvec2rhophi();
	return *this;
}

PulseFreq & PulseFreq::operator-=(const PulseFreq &rhs){
	DataOps::diff(cvec,rhs.cvec,samples);
	cvec2rhophi();
        return *this;
}

PulseFreq & PulseFreq::operator*=(const PulseFreq &rhs){
	DataOps::mul(cvec,rhs.cvec,samples);
	cvec2rhophi();
	return *this;
}

PulseFreq & PulseFreq::operator*=(const double s){
	DataOps::mul(cvec,s,samples);
	cvec2rhophi();
	return *this;
}

PulseFreq & PulseFreq::operator/=(const PulseFreq &rhs){
	DataOps::div(cvec,rhs.cvec,samples);
	cvec2rhophi();
        return *this;
}

PulseFreq & PulseFreq::interfere(const PulseFreq &rhs){
	*this += rhs;
        return *this;
}

PulseFreq & PulseFreq::diffamps(const PulseFreq &rhs){
	rhovec -= rhs.rhovec;
	std::fill(phivec.begin(),phivec.end(),0);
	rhophi2cvec();
        return *this;
}

PulseFreq & PulseFreq::normamps(const PulseFreq &rhs){
	rhovec /= rhs.rhovec;
	std::fill(phivec.begin(),phivec.end(),0);
	rhophi2cvec();
        return *this;
}

void PulseFreq::print_amp(std::ofstream & outfile)
{
	outfile << "# amp\n";
	outfile << rhovec << std::endl;
}
void PulseFreq::print_phase(std::ofstream & outfile)
{
	outfile << "# phase\n";
	outfile << phivec << std::endl;
}
void PulseFreq::print_phase_powerspectrum(std::ofstream & outfile)
{
/*
	double * phase = (double *) fftw_malloc(sizeof(double) * samples);
	double * phaseFT = (double *) fftw_malloc(sizeof(double) * samples);
	fftw_plan plan_r2hc = fftw_plan_r2r_1d(samples,
			phase,
			phaseFT,
			FFTW_R2HC,
			FFTW_MEASURE
			);
*/
	std::copy(phivec.begin(),phivec.end(),r_vec);
	fftw_execute_r2r(*FTplan_r2hcPtr.get(),r_vec,hc_vecFT);
	

	outfile << "# power spectrum of the Fourier phase\n";
	outfile << std::pow(hc_vecFT[0],int(2)) << "\n";
	for (size_t i = 1; i<samples/2;++i){
		outfile << std::pow(hc_vecFT[i],int(2)) + std::pow(hc_vecFT[samples-i],int(2)) << "\n";
	}
	outfile << std::pow(hc_vecFT[samples/2],int(2)) << std::endl;
	outfile << std::endl;
}

bool PulseFreq::addrandomphase(void)
{
	if (!infreq){
		std::cerr << "died here at addrandomphase()" << std::endl;
		return false;
	}
	size_t sz = samples*2; // doubling the vector to mirror it so that DFT hansles the phase well

	double * randphase = (double *) fftw_malloc(sizeof(double) * sz);
	double * randphaseFT = (double *) fftw_malloc(sizeof(double) * sz);
/*

	fftw_plan plan_r2hc = fftw_plan_r2r_1d(sz,
			randphase,
			randphaseFT,
			FFTW_R2HC,
			FFTW_MEASURE
			);
	fftw_plan plan_hc2r = fftw_plan_r2r_1d(sz,
			randphaseFT,
			randphase,
			FFTW_HC2R,
			FFTW_MEASURE
			);

*/


	std::uniform_real_distribution<double> distribution(
		(double(atof(getenv("randphase_mean")))-double(atof(getenv("randphase_std"))))*Constants::pi<double>(),
		(double(atof(getenv("randphase_mean")))+double(atof(getenv("randphase_std"))))*Constants::pi<double>()
		);

	double phase = distribution(rng);
	randphase[0] = phase;
	randphase[samples/2] = -phase;
	for (size_t i = 1; i<samples/2;i++){
		phase = distribution(rng);
		randphase[i] = phase;
		randphase[samples-i] = -phase;
	}
	for (size_t i=sz-1;i>sz/2-1;--i){
		randphase[i] = randphase[sz-i];
	}

	size_t lowpass = boost::lexical_cast<size_t>(atoi(getenv("phaseNoiseLowpass")));
	std::cerr << "\n======== lowpass is " << lowpass << " =======\n" << std::flush;

	fftw_execute_r2r(*FTplan_r2hc_2xPtr.get(),randphase,randphaseFT);
	std::fill(randphaseFT+lowpass,randphaseFT+sz-lowpass,0.);
	for (size_t i=1;i<lowpass;++i){
		double filter = std::pow(std::cos(double(i)/(double(lowpass)) * Constants::half_pi<double>() ),int(2));
		randphaseFT[i] *= filter;
		randphaseFT[sz-i] *= filter;
	}
	randphaseFT[sz/2] = 0.;
	fftw_execute_r2r(*FTplan_hc2r_2xPtr.get(),randphaseFT,randphase);

	for (size_t i=0;i<samples;++i){
		phivec[i] += randphase[i]/samples;
	}
	rhophi2cvec();
	return true;
}

void PulseFreq::attenuate(double attenfactor){
	rhovec *= attenfactor;
	rhophi2cvec();
}
void PulseFreq::phase(double phasein){ // expects delay in units of pi , i.e. 1.0 = pi phase flip 
	if(intime){
		fft_tofreq();
	}
	phivec += phasein*Constants::pi<double>();
	rhophi2cvec();
	if(intime){
		fft_totime();
	}
}
void PulseFreq::delay(double delayin){ // expects delay in fs
	if(intime){
		fft_tofreq();
	}
	for (unsigned i=0;i<samples;i++){
		phivec[i] += omega[i]*delayin/fsPau<double>();
	}
	rhophi2cvec();
	if(intime){
		fft_totime();
	}
}


void PulseFreq::printfrequency(std::ofstream * outfile){
	double nu,lambda;
	(*outfile) << ("#omega[fs^-1]\tlambda[nm]\trho\tphi\n");
	for (unsigned i = i_low;i<i_high;i+=(unsigned)(atoi(getenv("sampleinterval")))){
		nu = (omega[i]/(2.0*pi<double>())/fsPau<double>());
		lambda = C_nmPfs<double>()/nu;
		(*outfile) << nu << "\t" << lambda << "\t" << std::pow(rhovec[i],int(2)) << "\t" << phivec[i] << "\n";
	}
}
void PulseFreq::printwavelengthbins(std::ofstream * outfile)
{
	std::vector<double> x(2);
	x.front() = C_nmPfs<double>()*2.0*pi<double>()*fsPau<double>()/omega[i_low];
	x.back() = C_nmPfs<double>()*2.0*pi<double>()*fsPau<double>()/omega[i_high-1];
	double dlam = (x.front()-x.back())/double(m_lamsamples);
        for (size_t i = 0;i<m_lamsamples;++i){
		(*outfile) << x.back() + i*dlam << "\t";
        }
	(*outfile) << "\n";
	return;
}
void PulseFreq::appendwavelength(std::ofstream * outfile)
{
	std::vector<double> x(i_high-i_low);
	std::vector<double> y(i_high-i_low);	
	for (size_t i=0;i<y.size();++i){
		x[i] = C_nmPfs<double>()*2.0*pi<double>()*fsPau<double>()/omega[i_low+i];
		//y[i] = std::pow(rhovec[i_low+i],int(2)) * 200000000000;
		y[i] = std::min(std::pow(rhovec[i_low+i],int(2)) * m_gain,double(m_saturate));
	}
	double dlam = (x.front()-x.back())/double(m_lamsamples);
	boost::math::barycentric_rational<double> interpolant(x.data(), y.data(), y.size());
	for (size_t i=0;i<m_lamsamples;++i){
		(*outfile) << uint16_t(interpolant(x.back()+i*dlam)) << "\t";
	}
	(*outfile) << std::endl;
	return;
}
void PulseFreq::appendfrequency(std::ofstream * outfile){
        for (unsigned i = i_low;i<i_high;i+=m_sampleinterval){ 
		uint16_t val = std::min(uint16_t(rhovec[i] * m_gain),uint16_t(m_saturate));
       		(*outfile) << std::pow(val,int(2)) << "\t";
        }
	(*outfile) << std::endl;
}

void PulseFreq::appendnoisy(std::ofstream * outfile){
	std::normal_distribution<double> norm_dist( 0.0, m_noisescale);
        for (unsigned i = i_low;i<i_high;i+=m_sampleinterval){ 
		double outval = std::pow(rhovec[i],int(2)) + norm_dist(rng);
       		(*outfile) << outval << "\t" ;
        }
	(*outfile) << std::endl;
}

void PulseFreq::printfrequencybins(std::ofstream * outfile){
	double nu,lam;
        for (unsigned i = i_low;i<i_high;i+=(unsigned)(atoi(getenv("sampleinterval")))){
		nu = (omega[i]/(2.0*pi<double>())/fsPau<double>());
		lam = C_nmPfs<double>()/nu;
       		(*outfile) << nu << "\t" << lam << "\n";
        }
	(*outfile) << "\n";
}
void PulseFreq::appendfrequencybins(std::ofstream * outfile){
	double nu,lam;
        for (unsigned i = i_low;i<i_high;i+=(unsigned)(atoi(getenv("sampleinterval")))){
		nu = (omega[i]/(2.0*pi<double>())/fsPau<double>());
       		(*outfile) << nu << "\t";
        }
	(*outfile) << "\n";
}
void PulseFreq::printfrequencydelay(std::ofstream * outfile, const double *delay){
	double nu,lambda,thisdelay;
	thisdelay = (*delay)*fsPau<double>();
	(*outfile) << ("#omega[fs^-1]\tlambda[nm]\trho\tphi\tdelay[fs]\n");
	for (unsigned i = i_low;i<i_high;i+=(unsigned)(atoi(getenv("sampleinterval")))){
		nu = (omega[i]/(2.0*pi<double>())/fsPau<double>());
		lambda = C_nmPfs<double>()/nu;
		(*outfile) << nu << "\t" << lambda << "\t" << rhovec[i] << "\t" << phivec[i] << "\t" << thisdelay << "\n";
	}
	(*outfile) << "\n";
}
void PulseFreq::printfrequencydelaychirp(std::ofstream * outfile, const double *delay,const double *chirp){
	double nu,lambda,thisdelay,thischirp,reltime;
	thisdelay = (*delay)*fsPau<double>();
	thischirp = (*chirp)*std::pow(fsPau<double>(),int(2));
	(*outfile) << ("#omega[fs^-1]\tlambda[nm]\trho\tphi\tdelay[fs]\tchirp[fs^2]\treldelays[fs]\n");
	for (unsigned i = i_low;i<i_high;i+=(unsigned)(atoi(getenv("sampleinterval")))){
		nu = (omega[i]/(2.0*pi<double>())/fsPau<double>());
		lambda = C_nmPfs<double>()/nu;
		reltime = (nu-nu0)*thischirp*10.0;
		(*outfile) << nu << "\t" << lambda << "\t" << rhovec[i] << "\t" << phivec[i] << "\t" << thisdelay << "\t" << thischirp << "\t" << reltime << "\n";
	}
	(*outfile) << "\n";
}

void PulseFreq::printwavelength(std::ofstream * outfile,const double *delay){
	double nu,lambda,lamlast;
	lamlast=0;
        for (unsigned i = i_high;i>i_low;i-=(unsigned)(atoi(getenv("sampleinterval")))){
                nu = (omega[i]/(2.0*pi<double>())/fsPau<double>());
                lambda = (double)((int)((C_nmPfs<double>()/nu)*10.0))/10.0;
		if (lambda>lamlast & (int)(lambda*10)%5 == 0){
                	(*outfile) << lambda << "\t" << (*delay) << "\t" << std::pow(rhovec[i],int(2)) << "\n";
			lamlast = lambda;
		}
        }
	(*outfile) << "\n";
}

void PulseFreq::printtime(std::ofstream * outfile){
	(*outfile) << ("#time\treal\timag\n");
	for (unsigned i = samples/2;i<samples;i++){
		(*outfile) << (time[i]*fsPau<double>()) << "\t" << cvec[i].real() << "\t" << cvec[i].imag() << "\n";
	}
	for (unsigned i = 0;i<samples/2;i++){
		(*outfile) << (time[i]*fsPau<double>()) << "\t" << cvec[i].real() << "\t" << cvec[i].imag() << "\n";
	}

} 


void PulseFreq::buildvectors(void){
	//std::cerr << "allocating with fftw_malloc with samples = " << samples << std::endl;
	cvec = (std::complex<double> *) fftw_malloc(sizeof(std::complex<double>) * samples);
        std::fill(cvec,cvec + samples,std::complex<double>(0));
	r_vec = (double *) fftw_malloc(sizeof(double) * samples);
        std::fill(r_vec,r_vec + samples,double(0));
	//std::cerr << "\t\t...allocated with fftw_malloc with samples = " << samples << std::endl;
	hc_vecFT = (double *) fftw_malloc(sizeof(double) * samples);
        std::fill(hc_vecFT,hc_vecFT + samples,double(0));
	//std::cerr << "allocating with fftw_malloc with samples = " << (2*samples) << std::endl;
	r_vec_2x = (double *) fftw_malloc(sizeof(double) * samples * 2);
        //std::fill(r_vec_2x,r_vec_2x + 2*samples,double(0));
	hc_vec_2xFT = (double *) fftw_malloc(sizeof(double) * samples * 2);
        std::fill(hc_vec_2xFT,hc_vec_2xFT + 2*samples,double(0));
	//std::cerr << "\t\t...allocated with fftw_malloc with 2*samples = " << (2*samples) << std::endl;

	rhovec.resize(samples,0.0);
	phivec.resize(samples,0.0);
	modamp.resize(samples,1.0);
	modphase.resize(samples,0.0);
	omega.resize(samples);
	time.resize(samples);

	omega[0] = 0.0;
	time[0] = 0.0;
	omega[samples/2] = -(double)(samples/2)*domega;
	time[samples/2] = -(double)(samples/2)*dtime;

	startind = (unsigned)((omega_center-(omega_width/2.0))/domega);
	stopind = (unsigned)((omega_center+(omega_width/2.0))/domega);
	onwidth = (unsigned)(omega_onwidth/domega); // 2.0)/domega);//  /10.0)/domega);// sin^2 goes from0..1 in 0..pi/2
	offwidth = (unsigned)(omega_offwidth/domega); // 2.0)/domega);//  /10.0)/domega);// sin^2 goes from0..1 in 0..pi/2
	for (unsigned i = 1; i<startind;i++){
		omega[i] = domega*i;
		time[i] = dtime*i;
		omega[samples-i] = -domega*i;
		time[samples-i] = -dtime*i;
	}
	for (unsigned i = startind;i<startind+onwidth; i++){
		omega[i] = domega*i;
		time[i] = dtime*i;
		rhovec[i] = rising(i);
		cvec[i] = std::polar(rhovec[i],phivec[i]);
		omega[samples-i] = -domega*(double)i;
		time[samples-i] = -dtime*(double)i;
		rhovec[samples-i] = rising(i);
		cvec[samples-i] = std::polar(rhovec[samples-i],phivec[samples-i]);
	}
	for (unsigned i = startind+onwidth;i<stopind-offwidth; i++){
		omega[i] = domega*i;
		time[i] = dtime*i;
		rhovec[i] = 1.0;
		cvec[i] = std::polar(rhovec[i],phivec[i]);
		omega[samples-i] = -domega*(double)i;
		time[samples-i] = -dtime*(double)i;
		rhovec[samples-i] = 1.0;
		cvec[samples-i] = std::polar(rhovec[samples-i],phivec[samples-i]);
	}
	for (unsigned i = stopind-offwidth;i<stopind; i++){
		omega[i] = domega*i;
		time[i] = dtime*i;
		rhovec[i] = falling(i);
		cvec[i] = std::polar(rhovec[i],phivec[i]);
		omega[samples-i] = -domega*i;
		time[samples-i] = -dtime*i;
		rhovec[samples-i] = falling(i);
		cvec[samples - i] = std::polar(rhovec[samples - i],phivec[samples - i]);
	}
	for (unsigned i = stopind;i<samples/2; i++){
		omega[i] = domega*i;
		time[i] = dtime*i;
		omega[samples-i] = -domega*i;
		time[samples-i] = -dtime*i;
	}
	
}
void PulseFreq::killvectors(void){
	fftw_free(cvec);
	fftw_free(r_vec);
	fftw_free(hc_vecFT);
	fftw_free(r_vec_2x);
	fftw_free(hc_vec_2xFT);
	cvec = NULL;
	r_vec = hc_vecFT = r_vec_2x = hc_vec_2xFT = NULL;
}

void PulseFreq::setplans(const PulseFreq & rhs)
{
	FTplan_forwardPtr = rhs.FTplan_forwardPtr;
	FTplan_backwardPtr = rhs.FTplan_backwardPtr;
}
void PulseFreq::setmasterplans(fftw_plan * const forward,fftw_plan * const backward)
{
	assert(FTplan_forwardPtr.use_count()==0 && FTplan_backwardPtr.use_count()==0);
	*forward = fftw_plan_dft_1d(samples, 
			reinterpret_cast<fftw_complex*>(cvec),
			reinterpret_cast<fftw_complex*>(cvec), 
			FFTW_FORWARD, FFTW_ESTIMATE);
	*backward = fftw_plan_dft_1d(samples, 
			reinterpret_cast<fftw_complex*>(cvec), 
			reinterpret_cast<fftw_complex*>(cvec), 
			FFTW_BACKWARD, FFTW_ESTIMATE);
	FTplan_forwardPtr = std::make_shared<fftw_plan> (*forward);
	FTplan_backwardPtr = std::make_shared<fftw_plan> (*backward);
}
void PulseFreq::setancillaryplans(fftw_plan * const r2hc,fftw_plan * const hc2r,fftw_plan * const r2hc_2x,fftw_plan * const hc2r_2x)
{

	assert(FTplan_r2hcPtr.use_count()==0
			&& FTplan_hc2rPtr.use_count()==0
			&& FTplan_r2hc_2xPtr.use_count()==0
			&& FTplan_hc2r_2xPtr.use_count()==0);
	*r2hc = fftw_plan_r2r_1d(samples,
			r_vec,
			hc_vecFT,
			FFTW_R2HC,
			FFTW_MEASURE
			);
	*hc2r = fftw_plan_r2r_1d(samples,
			hc_vecFT,
			r_vec,
			FFTW_HC2R,
			FFTW_MEASURE
			);
	*r2hc_2x = fftw_plan_r2r_1d(2*samples,
			r_vec_2x,
			hc_vec_2xFT,
			FFTW_R2HC,
			FFTW_MEASURE
			);
	*hc2r_2x = fftw_plan_r2r_1d(2*samples,
			hc_vec_2xFT,
			r_vec_2x,
			FFTW_HC2R,
			FFTW_MEASURE
			);
	FTplan_r2hcPtr = std::make_shared<fftw_plan> (*r2hc);
	FTplan_hc2rPtr = std::make_shared<fftw_plan> (*hc2r);
	FTplan_r2hc_2xPtr = std::make_shared<fftw_plan> (*r2hc_2x);
	FTplan_hc2r_2xPtr = std::make_shared<fftw_plan> (*hc2r_2x);

}

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
