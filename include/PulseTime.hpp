// Pulse class definition

#ifndef PULSETIME_H
#define PULSETIME_H


// standard includes
#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>

// my headers
#include <Constants.hpp>

// my definitions

class PulseTime {

	public:
		PulseTime(double strength_in = 1e-3 * 0.696, double width_in = 50, double t0_in = 0.0) : 
			strength(strength_in * Constants::auenergy<double>()/Constants::Eh<double>() * std::pow(Constants::aufor10PW<double>(),unsigned(2))), 
			Ctau(width_in * Constants::root_pi<double>() / Constants::fsPau<double>() / 2.0),
			t0(t0_in / Constants::fsPau<double>())
	{
		//    clog << "Creating Pulse " << this << endl;
	}
		~PulseTime()
		{
			//    clog << "Destroying Pulse " << this << endl;
		}

		void setstrength(const double in);
		void scalestrength(const double in);
		void setwidth(const double in);
		void sett0(const double in);
		void adddelay(const double in);

		inline double getstrength() { return strength; }
		inline double getCtau() { return Ctau; }
		inline double gett0() { return t0; }

		inline bool getenvelope(const double t,double *FF,double *dFFdt) 
		{
			if ( inpulse(t) ){
				*FF = strength * ( gsl_pow_2( cos(Constants::half_pi<double>()*(t-t0)/Ctau) ) );
				*dFFdt = -strength/2 * ( Constants::pi<double>()/Ctau * sin(Constants::pi<double>()*(t-t0)/Ctau));
				return true;
			} else {
				*FF = 0.0;
				*dFFdt = 0.0;
				return false;
			}
		}
		inline bool getenvelope(const double t,double *FF) 
		{
			if ( inpulse(t) ){
				*FF = strength * ( gsl_pow_2( cos(Constants::half_pi<double>()*(t-t0)/Ctau) ) );
				return true;
			} else {
				*FF = 0.0;
				return false;
			}
		}

	private:
		double strength, Ctau, t0;

		inline bool inpulse(const double t) 
		{
			if (t >= -Ctau && t <= Ctau){
				return true;
			} else {
				return false;
			}
		}
};

#endif
