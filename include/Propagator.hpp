#include "Propagator.hpp"

Params::Params(void)
{

}

Params::~Params(void)
{
	delete jPtr;
	delete vibsptr;
	delete vibensPtr;
	delete pvPtr;
	delete ejPtr;
	delete pjPtr;
	delete pjnewPtr;
	delete aajmPtr;
	delete bbjmPtr;
	delete ccjmPtr;
}

Propagator::Propagator(Params & params,jEnsemble & jensemble, vEnsemble & vensemble, PulseTime * pulsesPtr)
: abstol(1e-10)
, reltol(1e-8)
{

	paramsPtr = new PARAMS;

  	// setting up variables //
	// Ugh... now we should move all this over to the classes //
	paramsptr->jPtr = new unsigned[sizej];
 	paramsPtr->vibsPtr = new unsigned[nvibs];
	paramsPtr->vibensPtr = new float[nvibs];
	paramsPtr->pvPtr = new float[nvibs];
	paramsPtr->ejPtr = new double[sizej];
	paramsPtr->pjPtr = new float[sizej];
	paramsPtr->pjnewPtr = new float[sizej];
	
	paramsPtr->aajmPtr = new double[sizej];
	paramsPtr->bbjmPtr = new double[sizej];
	paramsPtr->ccjmPtr = new double[sizej];
	paramsPtr->pulsesPtr = new Pulses[npulses];

}

Propagator::~Propagator(void)
{

	delete paramsPtr;

}
