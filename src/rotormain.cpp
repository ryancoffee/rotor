#include <rotormain.hpp>

// standard includes
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <string>
#include <vector>

/*
// gsl includes
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_odeiv.h>
*/

// my headers
#include <Constants.hpp>
#include <Molecules.hpp>
#include <jEnsemble.hpp>
#include <vEnsemble.hpp>
#include <Propagator.hpp>

// my preprocessor defs

int main(int argc, char *argv[]) {

	// Checking that environment variables have been set
	if (! (getenv("njs") && getenv("nvibs"))) {
		std::cerr << "first source set_vars rather than original syntax of \n" 
		<< "$ main njs T(K) twinstart(ps) tend(ps) tstep (ps) numpulses strength(10^16W/cm^2) decay(>0) width(fs) delay(fs)" << std::endl;
		std::cerr << "Writes files pj_nn.dat, cossq_nn.dat, field_nn.dat" << std::endl;
		return 1;
	}

	using namespace Constants;

	// generate filenames //
	std::string datadir = getenv("datadir");
	std::string filename;
	std::string pjstarttail="_";
	std::string filetail = "_";
	filetail += getenv("TinKelvin");//argv[2];
	filetail += "K_";
	filetail += getenv("npulses"); //argv[6];
	filetail += "pls_";
	filetail += getenv("intensity");//argv[7];
	filetail += "inten_";
	filetail += getenv("fswidth");//argv[9];
	filetail += "fswidth_";
	filetail += getenv("fsdelay");//argv[10];
	filetail += "fsdelay.dat";

	pjstarttail += getenv("TinKelvin");//argv[2];
	pjstarttail += "K.dat";

	// Start the clock //
	time_t start,end,rawtime;
	double diff;
	time(&start);
	filename = "kicklog";
	filename += filetail;
	std::ofstream logout(filename.c_str(),std::ios::app);
	if (!logout){
		std::cerr << "Couldn't append to " << filename << std::endl;
		return 1;
	}

	time(&rawtime);
	logout << ctime(&rawtime) << " $ here we still need to write the parameters that are now being passed via environment variables" ;
	logout << std::endl;

	const size_t nthreads((size_t)atoi(getenv("nthreads")));
	const double kTinau((float)(atof(getenv("TinKelvin")))*kb<float>()/Eh<float>());
	Molecules molecule(double(kTinau),Molecules::nno); 
	std::cout << "The molecules ID is:\t" << molecule.printMolID() << std::endl;

	size_t nvibs((size_t)atoi(getenv("nvibs"))); 
	std::cerr << "\n\n\t\t\t --- HERE HERE HERE HERE --- \n\n" << std::flush;
	return 0;
}

/*
	vEnsemble vsystem(nvibs,molecule);
	nvibs = vsystem.initdist();

	const size_t njs(atoi(getenv("njs")));
	std::vector< jEnsemble<double> > jsystem;
	jsystem.reserve(nvibs);
	for (unsigned v=0; v<jsystem.size();++v){
		jsystem.push_back(njs,v,molecule);
		jsystem.fill6js();
	}

	const unsigned npulses = (unsigned)std::stoi(argv[6]);
	std::vector<PulseTime> pulses<double>(npulses);
	for (unsigned i=0;i<pulses.size();++i){
		pulses[i].setstrength(molecule.delta_alpha()*auenergy<double>()/Eh<double>() * std::pow(aufor10PW<double>(),int(2)) * float(atof(getenv("intensity"))) );
		pulses[i].scalestrength(std::pow(float(atof(getenv(scalestrength))),//argv[8]),
					i));
		pulses[i].setwidth(root_pi<double>() * boost::lexical_cast<double>(atof(getenv("pulsewidth")))/fsPau/2.0;);
		pulses[i].sett0(double(i+1) * double(-atof(getenv("delay")))/fsPau<double>()); // run pulses back to negative time so 1st full is near 0
	}




	// setting up local variables //
	int j[sizej],vibs[nvibs],newvibsize;
	double vibens[nvibs],pv[nvibs];
	double ej[sizej],pj[sizej],pjnew[sizej];//,pjmnew[sizej][2*sizej-1];
	//  double sqrt2jp1[sizej+2];
	double aajm[sizej],bbjm[sizej],ccjm[sizej];
	
	//	HERE HERE HERE HERE 	//
	//	Still here... working on 6j symbols 	//

	// I GOT HERE //
	// need to make jEnsemble class and vEnsemble class //
	// Then make those two friend classes of Params and Propagator //
	// UGH... need to move all this stuff into the classes //
	std::vector<PulseTime> pulsesVec;
	for (unsigned i =0;i< (unsigned) atoi(argv[6]) ; i++){
		pulsesVec.emplace_back();
	}
	vEnsemble vensemble(nvibs);


	std::vector<Ensemble> jensembleVec;
	for (unsigned i=0;i<nvibs;i++){
		jensembleVec.emplace_back( (unsigned)atoi(argv[1]) ); // this is guaranteed to fail... trying to construct at list addition
	}
	std::vector<Ensemble>::const_iter eItr;
	std::vector<Ensemble>::const_iter eEnd;
	for (eItr = jensembleVec.begin(); eItr != jensembleVec.end(); eItr++) {
		eItr.setkT((float)atof(argv[2])); // this too must be wrong, somehow trying no tto use pointers is likely a bad mistake.
	}
	Params params;
	Propagator propagator(params, jensemble, vensemble, pulsesPtr);

	// this would likely be better if we did virtual functions and had j and vensembles inherit from a general ensemble class //
	// OK here is how we could handle this.  Each ensemble j and v have common interfaces jus tlike in leap stuff.
	// then we make a list of j ensembles, one for each m state that we want to worry about.
	// somehow, kT tells us which m's we need to worry about... count m's up to pj < 1e-3 or something... set the min jpop in its header
	jensemble.setkT((float)atof(argv[2]));
	vensemble.setkT((float)atof(argv[2]));


	params.strengths = strengths;
	params.Ctaus = Ctaus;
	params.t0s = t0s;
	params.phases = phases;
	params.omegas = omegas;
	params.npulses = npulses;


	params.kTinau=kTinau;
	params.jPtr = jPtr;
	params.ejPtr = ejPtr;
	params.pjPtr = pjPtr;
	params.pjnewPtr = pjnewPtr;
	params.aajmPtr=aajmPtr;
	params.bbjmPtr=bbjmPtr;
	params.ccjmPtr=ccjmPtr;

	params.vibsPtr=vibsPtr;
	params.vibensPtr=vibensPtr;
	params.pvPtr=pvPtr;
	params.nvibs=nvibs;
	params.newvibsize=newvibsize;

	// Init vibs //
	setvibsvibens(&params);
	setpvvec(&params);

	params.maxj=maxj;
	params.sizej=sizej;

	const double twinstart = static_cast<double>(atof(argv[3]))*1e3/fsPau;
	const double tend = static_cast<double>(atof(argv[4]))*1e3/fsPau;
	if (tend<twinstart){
		cerr << "Whoops, twinstart is greater than tend... avoid a segmentation fault, fix this!" << endl;
		return 1;
	};

	const double tstepsize = static_cast<double>(atof(argv[5]))*1e3/fsPau;
	double tstart = twinstart;
	for (int i=0;i<npulses;i++){
		tstart = GSL_MIN_DBL(params.t0s[i] - params.Ctaus[i],tstart);
	}
	double tmp  =  (( tend - twinstart )/ tstepsize) +1;
	int ntsteps = static_cast<int>(tmp);
	params.ntsteps = ntsteps;
	params.tstepsize = tstepsize;
	clog << twinstart*fsPau << " fs: " << tstepsize*fsPau << " fs: " << tend*fsPau << " fs in " << ntsteps << " steps." << endl;

	// ------------- setting up time and signal vectors ---------------- //

	int tstep=0;
	double times[2*ntsteps],signal[2*ntsteps],imsignal[2*ntsteps];//,cossq[2]={0.0,0.0};
	for(tstep=0;tstep<ntsteps;tstep++){
		signal[tstep]=0.0;
		imsignal[tstep]=0.0;
		times[tstep]=twinstart+tstep*tstepsize;
	}
	tstep=0;


	// ----- testing the setrealj function vvv
	   //params.m=2;
	   //params.jstart=3;
	   //int i=2;
	   //int testj;
	   //setrealj(&testj,&i,&params);
	   //cout << "m=" << params.m << "\tjstart=" << params.jstart << "\ti=" << i << "\t gives:\trealj =" << testj << "\n";
	// ------ testing setrealj function ^^^

	// ---------- setting up the ODE system ---------- //

	// initializing evolution system
	double passt=tstart;                     // - delay -  (1/e half-widths) of pulse1  
	//  clog << "t starts at " << passt*fsPau << "fs" << endl;


	// --------- Init pjnew vec --------//
	for (int i=0;i<sizej;i++){
		pjnewPtr[i] = 0.0;
		//      for (int m = 0 ;m<(2*sizej-1);m++){
		//        pjmnew[i][m] = 0.0;
		//      }
	}

	//  params.newvibsize=1; // use to evaluate only v=0
	for (int v=0;v<params.newvibsize;v++) {    // ultimately loop on v
		params.currentv=v;
		setjejvecs(&params);
		setpjvec(&params);


		double pop=1.0; //for testing the preservation of square magnitude
		const double poptol=0.5; // 	used in if (pop > 1+poptol || pop < 1-poptol){ ... 

		for(int m=0;m<sizej-2;m++) {                                               // ultimately loop on m
			//    clog << "working on m = " << m << "  "; 

			params.m=m;

			setelemsjm(m,sizej,aajmPtr,bbjmPtr,ccjmPtr);   // ,coeffaamPtr,coeffbbmPtr,coeffccmPtr,aaj0Ptr,bbj0Ptr,ccj0Ptr,sqrt2jp1Ptr);

			// ------- now I want to make a function to set dim ----------- //
			int dim=getmaxdim(&params); //(sizej-m)*2;// (jwin*2+1)*2;  // here's where we remove jwin, it seems dim is twice too big since only even(odd) Js couple

			params.dim=dim;
			params.ejPtr=ejPtr;
			const gsl_odeiv_step_type *T = gsl_odeiv_step_rk8pd;  // fast step
			//    const gsl_odeiv_step_type *T = gsl_odeiv_step_rkf45;  // smooth step
			gsl_odeiv_step *s = gsl_odeiv_step_alloc(T,dim);
			gsl_odeiv_control *c = gsl_odeiv_control_y_new(abstol,reltol);
			gsl_odeiv_evolve *e = gsl_odeiv_evolve_alloc(dim);
			gsl_odeiv_system sys = {func,jac,dim,paraPtrvoid};

			double h = .1/fsPau;                                   // initial step size
			double y[dim];        

			for(int jstart = 0;jstart<sizej-m;jstart++) {           // ultimately loop on j
				params.jstart = jstart;
				if (jstart+m<10 || pjPtr[jstart+m]>5e-4) {

					double t=tstart;
					inity(y,&params);



					int status=GSL_SUCCESS,tind=0;
					for(tind=0;tind<ntsteps;tind++){
						if (!nearpulse(t,&params)){             // somehow do a freepropstep in here

							freepropstep(&t,&tstepsize,y,&params);

						} else {                              // solve ODEs for coeffs

							while (t<times[tind] && status==GSL_SUCCESS){
								status = gsl_odeiv_evolve_apply(e,c,s,&sys,&t,times[tind],&h,y);
								sqrnormalizey(y,&params);
								//	    fieldout << t*fsPau/1e3 << "\t" << params.field << endl;
								if (status !=GSL_SUCCESS){
									cerr << "evolve status: " << gsl_strerror(status) << endl;
									break;
								}
							}
						}

						*/

						size_t oldsize = jens.size();
						if (jens.checkgrowth()){molecule.fill(jens,oldsize)}

						/*

						addtosignal(y,signal,imsignal,tind,&params);

						passt = t;
					}
					//	clog << "exited while loop, status: " << status << " = " << gsl_strerror(status) << endl;

					addtopjnew(y,&params);
					if (pop> 1+poptol || pop < 1-poptol){
						cerr << "not preserving square magnitude for j" << params.jstart+params.m << ", m" << params.m << endl;
						cerr << "square magnitude = " << pop << endl;
						return 1;
					}

					if (m==8 && jstart==m && v==0 ){ // Test coupling for j,m=0,0 and j=26,m=13
						samplecoupling(y,&params);
					}
				} else {
					passtosignal(signal,imsignal,&params);
					passtopjnew(&params);
				}
			}

			// ---------- free the ODE memory ---------- //

			// needs to be in the m-loop since m determines the dimension... maybe later see if always keeping track of sizej equations is better.
			gsl_odeiv_evolve_free(e);
			gsl_odeiv_control_free(c);
			gsl_odeiv_step_free(s);

			//    clog << "done." << endl;    

		}

		time(&end);
		diff=difftime(end,start);
		clog << "just finished v = " << params.currentv << " at " << diff << "s." << endl;
		//  clog << "t ends at " << passt*fsPau << "fs" << endl;
	}




	// stop the clock
	time(&end);
	diff=difftime(end,start);
	clog << ntsteps << " timesteps, " << diff << " s runtime.\a" << endl;
	logout << ntsteps << " timesteps, " << diff << " s runtime." << endl;

	string pjstartname,pjnewname;
	pjstartname = "pjstart";
	pjstartname += pjstarttail;
	pjnewname = "pj";
	pjnewname += filetail;

	int wrote = writepops(&params,pjstartname,pjnewname);
	if (wrote !=0)
		return 1;


	//open output file
	filename = "cossq";
	filename += filetail;
	ofstream cossqout(filename.c_str(),ios::out);
	if (!cossqout){
		cerr << filename << " could not be opened" << endl;
		return 1;
	}

	for(int i=0;i<ntsteps;i++){
		cossqout << times[i]*fsPau/1e3 << " " << signal[i] << " " << imsignal[i] << endl;
	}


	// ---------- exit main ---------- //


	return 0;
}
*/







/*


	// ---------- member functions ---------- //



	void print2col(const int size, const int *intvec,  const double *doublevec){
		for(int i=0;i<size;i++){
			cout << *intvec << " " << *doublevec << endl;
			intvec++;
			doublevec++;
		}
		intvec-=size;
		doublevec-=size;
	}

	void print2col(const int size, const double *vec1, const double *vec2){
		for(int i=0;i<size;i++){
			cout << *vec1 << " " << *vec2 << endl;
			vec1++;
			vec2++;
		}
		vec1-=size;
		vec2-=size;
	}

	double sumvec(const int size,const double *vec){
		double sum=0;
		for(int i=0;i<size;i++){
			sum+=*vec;
			vec++;
		}
		vec -=size;
		return sum;
	}

	void scalevec(const int size,double *vec,double scale){
		for(int i=0;i<size;i++){
			*vec*=scale;
			vec++;
		}
		vec-=size;
	}

	void sumnormvec(const int size,double *vec){
		double sum=0;
		for(int i=0;i<size;i++){
			sum+=*vec;
			vec++;
		}
		vec -=size;
		scalevec(size,vec,1/sum);
	}

	// Initialize vib and rot vectors

	void setvibsvibens(PARAMS *paraPtr){
		// first set vibs via vibsPtr
		for(int vv=0;vv<paraPtr->nvibs;vv++){
			paraPtr->vibsPtr[vv]=vv;
		}
		// now set vibens via vibensPtr
		vibrationalenergies_nn(paraPtr);
	}

	void setjejvecs(PARAMS *paraPtr){
		// first set j via jPtr
		for(int i=0;i<paraPtr->sizej;i++){
			paraPtr->jPtr[i]=i;
		}
		// now set ej via ejPtr
		rotationalenergies_nn(paraPtr);
	}

	void setpvvec(PARAMS *paraPtr){//mark
		double *pvPtr = paraPtr->pvPtr;
		double frac;
		for(int i=0;i<paraPtr->nvibs;i++){
			frac=(paraPtr->vibensPtr[i]) / paraPtr->kTinau;
			pvPtr[i]=exp(-frac);
		}
		sumnormvec(paraPtr->nvibs,pvPtr);
		paraPtr->newvibsize=getnewsize(paraPtr->nvibs,pvPtr);
		sumnormvec(paraPtr->newvibsize,pvPtr);
		clog << "relevent vib pops = ";
		for(int i=0;i<paraPtr->newvibsize;i++){
			clog << pvPtr[i] << " ";
		}
		clog << endl;
	}


	}

	// Calculate vibrational energies:

	void vibrationalenergies_ii(PARAMS *paraPtr){
		const int nterms=10;
		int v=0,n=0;
		double vv=0;
		const double cc[]={214.5481, -0.616259, 7.507e-5, -1.263643e-4, 6.198129e-6, -2.0255975e-7, 3.9662824e-9, -4.6346554e-11, 2.9330755e-13, -7.61000e-16};
		const double *ccPtr=cc;
		for(v=0;v<paraPtr->nvibs;v++){
			vv=static_cast<double>(v);
			paraPtr->vibensPtr[v]=0;
			for(n=0;n<nterms;n++){
				paraPtr->vibensPtr[v]+=*ccPtr*gsl_pow_int(vv+0.5,n+1);
				ccPtr++;
			}
			ccPtr-=nterms;
			paraPtr->vibensPtr[v] /= icmPau;
		}
	}

	void vibrationalenergies_nn(PARAMS *paraPtr){
		const int n=min(11,paraPtr->nvibs);
		const double ens[]={1175.5, 3505.2, 5806.5, 8079.2, 10323.3, 12538.8, 14725.4, 16883.1, 19011.8, 21111.5, 23182.0};
		for(int v=0;v<n;v++){
			paraPtr->vibensPtr[v]=ens[v];
			paraPtr->vibensPtr[v] /=icmPau;
		}
	}
	void vibrationalenergies_oo(PARAMS *paraPtr){
		// from JOURNAL OF MOLECULAR SPECTROSCOPY 154,372-382 ( 1992) G. ROUILLE "High-Resolution Stimulated Raman Spectroscopy of O2"
		const int n=min(4,paraPtr->nvibs);
		const double env_0 = (1556.38991/2.0);
		const double dens[]={1556.38991, 1532.86724, 1509.5275};
		double ens[4];
		ens[0]=env_0;
		ens[1] = ens[0] + dens[0];
		ens[2] = ens[1] + dens[1];
		ens[3] = ens[2] + dens[2];
		for(int v=0;v<n;v++){
			paraPtr->vibensPtr[v]=ens[v];
			paraPtr->vibensPtr[v] /= icmPau;
		}
	}


	// Calculate rotational energies:
	void rotationalenergies_oo(PARAMS *paraPtr){
		// from JOURNAL OF MOLECULAR SPECTROSCOPY 154,372-382 ( 1992) G. ROUILLE "High-Resolution Stimulated Raman Spectroscopy of O2"
		const int vv = min(3,paraPtr->currentv);
		const double Bv[] = {1.437676476, 1.42186454, 1.4061199, 1.39042};
		const double Dv[] = {4.84256e-6, 4.8418e-6, 4.8410e-6, 4.8402e-6};
		const double Hv = 2.8e-12;
		double jj;
		for (int i=0;i<paraPtr->sizej;i++){
			jj= static_cast<double>(paraPtr->jPtr[i]);
			paraPtr->ejPtr[i] = Bv[vv]*jj*(jj+1)-Dv[vv]*gsl_pow_2(jj)*gsl_pow_2(jj+1) + Hv*gsl_pow_3(jj)*gsl_pow_3(jj+1); // in cm^-1
			paraPtr->ejPtr[i] /= icmPau; // in atomic units
		}
	}
	void rotationalenergies_nn(PARAMS *paraPtr){
		// numbers come from Loftus and Kuprienie, J. Phys. Chem. ref. Data, Vol. 6, No. 1, 1977. p242
		const int vv = min(16,paraPtr->currentv);
		const double Bv[]={1.98957, 1.972, 1.9548, 1.9374, 1.9200, 1.9022, 1.8845, 1.8666, 1.8488, 1.8310, 1.8131, 1.7956, 1.7771, 1.7590, 1.7406, 1.7223};
		double Dv = 1e-6*5.75, jj;
		//  Dv*=10; // artificially increase the distortion rather than temperature.
		//  Dv *= 0.0; // artificially kill distortion
		for (int i=0;i<paraPtr->sizej;i++){
			jj= static_cast<double>(paraPtr->jPtr[i]);
			paraPtr->ejPtr[i] = Bv[vv]*jj*(jj+1)-Dv*gsl_pow_2(jj)*gsl_pow_2(jj+1); // in cm^-1
			paraPtr->ejPtr[i] /= icmPau; // in atomic units
		}
	}

	void rotationalenergies_nn_nodistortion(PARAMS *paraPtr){
		// numbers come from Loftus and Kuprienie, J. Phys. Chem. ref. Data, Vol. 6, No. 1, 1977. p242
		const int vv = min(16,paraPtr->currentv);
		const double Bv[]={1.98957, 1.972, 1.9548, 1.9374, 1.9200, 1.9022, 1.8845, 1.8666, 1.8488, 1.8310, 1.8131, 1.7956, 1.7771, 1.7590, 1.7406, 1.7223};
		double jj;
		for (int i=0;i<paraPtr->sizej;i++){
			jj= static_cast<double>(paraPtr->jPtr[i]);
			paraPtr->ejPtr[i] = Bv[vv]*jj*(jj+1);// in cm^-1
			paraPtr->ejPtr[i] /= icmPau; // in atomic units
		}
	}

	void rotationalenergies_ii(PARAMS *paraPtr){
		double Bv, Dv, vv, jj;
		vv = static_cast<double>(paraPtr->currentv);
		for (int i=0;i<paraPtr->sizej;i++){
			jj= static_cast<double>(paraPtr->jPtr[i]);
			Bv=3.7395e-2 - 1.2435e-4*(vv+0.5) + 4.498e-7*gsl_pow_2(vv+0.5) - 1.482e-8*gsl_pow_3(vv+0.5) - 3.64e-11*gsl_pow_4(vv+0.5);
			Dv=4.54e-9 + 1.7e-11*(vv+0.5) + 7e-12*gsl_pow_2(vv+0.5);
			paraPtr->ejPtr[i]=Bv*jj*(jj+1)-Dv*gsl_pow_2(jj)*gsl_pow_2(jj+1); // in cm^-1
			paraPtr->ejPtr[i] /= icmPau; // in atomic units
		}
	}


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
	}


	void addtosignal(double *y,double *signal,double *imsignal,int tind,PARAMS *paraPtr){

		const int m=paraPtr->m;
		const int dim = paraPtr->dim;
		const int jstart = paraPtr->jstart;
		const int v = paraPtr->currentv;

		double realcossq=0.0, imagcossq=0.0;
		int realj;
		int *realjPtr=&realj;

		//  sqrnormalizey(y,paraPtr);

		for (int i=0;i<dim;i+=2){
			setrealj(realjPtr,&i,paraPtr);
			if (*realjPtr != -1){
				realcossq += ( gsl_pow_2(y[i]) + gsl_pow_2(y[i+1]) ) * paraPtr->aajmPtr[*realjPtr];
				imagcossq += 0.0;
				if(i<(dim-2)){
					realcossq += (y[i] * y[i+2] + y[i+1] * y[i+3])  * paraPtr->ccjmPtr[*realjPtr];
					imagcossq += (-y[i+3] * y[i] + y[i+2] * y[i+1]) * paraPtr->ccjmPtr[*realjPtr];
				}
				if(i>2){
					realcossq += (y[i] * y[i-2] + y[i+1] * y[i-1]) * paraPtr->bbjmPtr[*realjPtr];
					imagcossq += (-y[i-1] * y[i] + y[i-2] * y[i+1]) * paraPtr->bbjmPtr[*realjPtr];
				}
			}
		}

		if(m==0){
			signal[tind] += realcossq * paraPtr->pvPtr[v]  * paraPtr->pjPtr[jstart+m]/(2*(jstart+m)+1);
			imsignal[tind] += imagcossq * paraPtr->pvPtr[v]  * paraPtr->pjPtr[jstart+m]/(2*(jstart+m)+1);
		} else {
			signal[tind] += 2*realcossq * paraPtr->pvPtr[v]  * paraPtr->pjPtr[jstart+m]/(2*(jstart+m)+1);
			imsignal[tind] += 2*imagcossq * paraPtr->pvPtr[v]  * paraPtr->pjPtr[jstart+m]/(2*(jstart+m)+1);
		}
	}

	// here we can add addtocos() rather than addtocossq().  That would include terms that look forward and backward only by one element //


	void passtosignal(double *signal,double *imsignal,PARAMS *paraPtr){
		const int m = paraPtr->m;
		const int jstart = paraPtr->jstart;
		const int v = paraPtr->currentv;

		for (int tind=0;tind<paraPtr->ntsteps;tind++){
			if(m==0){
				signal[tind] += paraPtr->aajmPtr[jstart+m] * paraPtr->pvPtr[v] * paraPtr->pjPtr[jstart+m]/(2*(jstart+m)+1); // same as adding 1/3 * pj
			} else {
				signal[tind] += 2 * paraPtr->aajmPtr[jstart+m] * paraPtr->pvPtr[v] * paraPtr->pjPtr[jstart+m]/(2*(jstart+m)+1); // same as adding 2/3 * pj
			}
		}
	}

	void addtopjnew(const double *y,PARAMS *paraPtr){
		const int m=paraPtr->m;
		const int v=paraPtr->currentv;
		const int dim = paraPtr->dim;
		const int jstart = paraPtr->jstart;
		double *pjnewPtr = paraPtr->pjnewPtr;
		const double *pjPtr = paraPtr->pjPtr;
		const double *pvPtr = paraPtr->pvPtr;
		double sqrmag;

		int realj;
		int *realjPtr=&realj;

		for (int i=0;i<dim;i+=2){
			setrealj(realjPtr,&i,paraPtr);
			sqrmag = gsl_pow_2(y[i]) + gsl_pow_2(y[i+1]);
			if (m==0){
				pjnewPtr[*realjPtr] += sqrmag  * pvPtr[v] * pjPtr[jstart+m]/(2*(jstart+m)+1);
			} else {
				pjnewPtr[*realjPtr] += 2 * sqrmag * pvPtr[v] * pjPtr[jstart+m]/(2*(jstart+m)+1);
			}
		}
	}

	void passtopjnew(PARAMS *paraPtr){

		const int m=paraPtr->m;
		const int v=paraPtr->currentv;
		const int jstart = paraPtr->jstart;
		double *pjnewPtr = paraPtr->pjnewPtr;
		const double *pjPtr = paraPtr->pjPtr;
		const double *pvPtr = paraPtr->pvPtr;

		if (m==0){
			pjnewPtr[jstart+m] += pvPtr[v]* pjPtr[jstart+m]/(2*(jstart+m)+1);
		} else {
			pjnewPtr[jstart+m] += pvPtr[v]*2*pjPtr[jstart+m]/(2*(jstart+m)+1);
		}
	}

	void setrealj(int *realjPtr,const int *i,PARAMS *paraPtr){

		const int *m=&(paraPtr->m);
		const int *jstart = &(paraPtr->jstart);
		const int *maxj = &(paraPtr->maxj);

		*realjPtr = *i + *m; // + *jstart; is not right.
		*realjPtr += *jstart%2;

		if (*realjPtr<0 || *realjPtr>*maxj){
			*realjPtr=-1;
		}
	}

	void inity(double *y,PARAMS *paraPtr){
		const int *dim = &(paraPtr->dim);
		const int *jstart = &(paraPtr->jstart);
		//  const int *m = &(paraPtr->m);
		for(int i=0;i< *dim;i++){
			y[i]=0.0;                                // set ys to 0;
		}
		//  Now, we need to set the proper y[i] = 1.0;
		int indx;
		indx = 2 * ( *jstart / 2); // Grabs the integer part of j/2 so that jstart = 1 lands on index 0, and 3 on 2.
		y[indx] = 1.0;
	}


	void samplecoupling(double *y,PARAMS *paraPtr){
		int realj;
		clog << "\nlog10(population times 1e10):\n";
		for (int i=0;i<paraPtr->dim;i+=2){
			setrealj(&realj,&i,paraPtr); 
			double value = gsl_pow_2(y[i]) + gsl_pow_2(y[i+1]);
			value *= 1e10;
			int dots = static_cast<int>(log10(value));
			clog << realj << ":\t";
			for (int i = 0; i<dots; i++)
				clog << ".";
			clog << "\n";
		}
		clog << "\n";
	}

	void samplefield(double t,double FF){
		cout << t*fsPau/1e3 << "\t" << FF << "\n";   
	}

	int writepops(PARAMS *paraPtr,string pjstartname,string pjnewname){
		//open output file
		ofstream out(pjstartname.c_str(),ios::out);
		if (!out){
			cerr << pjstartname << " could not be opened" << endl;
			return 1;
		}
		for (int i=0;i<paraPtr->sizej;i++){
			out << i << " " << paraPtr->pjPtr[i] << endl;
		}

		//open output file
		ofstream pjout(pjnewname.c_str(),ios::out);
		if (!pjout){
			cerr << pjnewname << " could not be opened" << endl;
			return 1;
		}
		double pjnewsum = 0;
		for (int i=0;i<paraPtr->sizej;i++){
			pjout << i << " " << paraPtr->pjnewPtr[i] << endl;
			pjnewsum +=  paraPtr->pjnewPtr[i];
		}

		//  clog << "total pop is " << pjnewsum << endl;

		return 0;
	}

	void sqrnormalizey(double *y,PARAMS *paraPtr){
		double pop=0.0;
		for (int i=0;i<paraPtr->dim;i+=2){
			pop +=  gsl_pow_2(y[i]) + gsl_pow_2(y[i+1]);
		}
		scalevec(paraPtr->dim,y,1.0/sqrt(pop));
	}

	int getmaxdim(PARAMS *paraPtr){
		const int sizej = paraPtr->sizej;
		const int m = paraPtr->m;
		int out = sizej - m;
		out += out%2;
		return out;
	}


	void freepropstep(double *t,const double * dt, double y[],void *paraPtrvoid){
		PARAMS * paraPtr = (PARAMS *)paraPtrvoid;
		int jstart = paraPtr->jstart;
		int m = paraPtr->m;
		gsl_vector_view yview = gsl_vector_view_array(y,paraPtr->dim);
		size_t cdim = (paraPtr->dim)/2;
		gsl_vector_view yrealview = gsl_vector_subvector_with_stride(&yview.vector,0,2,cdim);
		gsl_vector_view yimagview = gsl_vector_subvector_with_stride(&yview.vector,1,2,cdim);

		gsl_complex z;

		gsl_vector_view ejvec = gsl_vector_view_array(paraPtr->ejPtr,paraPtr->sizej);

		//    gsl_vector_view energiesPart = gsl_vector_subvector_with_stride( &ejvec.vector,
		//  								   m+(jstart%2),
		//  								   2,
		//  								   cdim );
		double phase;
		for (size_t i = 0; i< cdim-1; i++){
			z = gsl_complex_rect(gsl_vector_get(&yrealview.vector,i),
					gsl_vector_get(&yimagview.vector,i) );
			phase = gsl_vector_get( &ejvec.vector , i*2 +jstart%2+m );
			phase -= gsl_vector_get( &ejvec.vector, (jstart+m) );
			phase *= -*dt;
			z = gsl_complex_mul(z,gsl_complex_polar(1.0,phase));
			gsl_vector_set(&yrealview.vector,i,GSL_REAL(z));
			gsl_vector_set(&yimagview.vector,i,GSL_IMAG(z));
			// add phase = -(Ej-Ejstart)dt to coeffs
			// I got here, now I need to wrap the phase by the energy.
		}
		*t += *dt;
		//  clog << "made it here: j/m = " << jstart << " / " << m;
		//  clog << " dim/cdim/sizeyrealview = " << paraPtr->dim << " / " << cdim << " / " << yrealview.vector.size << endl;
	}

	bool nearpulse(const double t,PARAMS * paraPtr){
		int p=0;
		for (p=0;p<paraPtr->npulses;p++){
			if (t >= paraPtr->t0s[p]-paraPtr->Ctaus[p] - 3*paraPtr->tstepsize && t <= paraPtr->t0s[p]+paraPtr->Ctaus[p] + 3*paraPtr->tstepsize){
				return true;
			}
		}
		return false;
	}

	bool inpulse(const double t,PARAMS * paraPtr){
		int p=0;
		for (p=0;p<paraPtr->npulses;p++){
			if (t >= paraPtr->t0s[p]-paraPtr->Ctaus[p] && t <= paraPtr->t0s[p]+paraPtr->Ctaus[p]){
				return true;
			}
		}
		return false;
	}
	bool inpulse(const double t,PARAMS * paraPtr,double *FF){
		bool isinpulse = false;
		int p=0;
		for (p=0;p<paraPtr->npulses;p++){
			if (t >= paraPtr->t0s[p]-paraPtr->Ctaus[p] && t <= paraPtr->t0s[p]+paraPtr->Ctaus[p]){
				isinpulse = true;
				*FF += paraPtr->strengths[p]* ( gsl_pow_2( cos(M_PI_2*(t-paraPtr->t0s[p])/paraPtr->Ctaus[p]) ) );
			}
		}
		return isinpulse;
	}
	bool inpulse(const double t,PARAMS * paraPtr,double *FF,double *dFFdt){
		bool isinpulse = false;
		int p=0;
		for (p=0;p<paraPtr->npulses;p++){
			if (t >= paraPtr->t0s[p]-paraPtr->Ctaus[p] && t <= paraPtr->t0s[p]+paraPtr->Ctaus[p]){
				isinpulse = true;
				*FF += paraPtr->strengths[p]* ( gsl_pow_2( cos(M_PI_2*(t-paraPtr->t0s[p])/paraPtr->Ctaus[p]) ) );
				*dFFdt += -paraPtr->strengths[p]/2 * ( M_PI/paraPtr->Ctaus[p] * sin(M_PI*(t-paraPtr->t0s[p])/paraPtr->Ctaus[p]));
			}
		}
		return isinpulse;
	}


*/
