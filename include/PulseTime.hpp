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

using namespace Constants;

class PulseTime {

	public:
		PulseTime(double strength_in = 1e-3 * 0.696, double width_in = 50, double t0_in = 0.0) : 
			strength(strength_in * auenergy<double>()/Eh<double>() * std::pow(aufor10PW<double>(),int(2))), 
			Ctau(width_in * root_pi<double>() / fsPau<double>() / 2.0),
			t0(t0_in / fsPau<double>())
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

		inline double getstrength(void) { return strength; }
		inline double getCtau(void) { return Ctau; }
		inline double gett0(void) { return t0; }
		inline double getduration(void) { return 2.0*Ctau; }

		inline bool getenvelope(const double t,double &FF,double &dFFdt) 
		{
			if ( !inpulse(t) ){
				FF = 0.0;
				dFFdt = 0.0;
				return false;
			} else {
				FF = strength * ( std::pow( std::cos(half_pi<double>()*(t-t0)/Ctau) , int(2)) );
				dFFdt = -strength/2 * ( pi<double>()/Ctau * std::sin(pi<double>()*(t-t0)/Ctau));
				return true;
			}
		}
		inline bool getenvelope(const double t,double &FF) 
		{
			if ( !inpulse(t) ){
				FF = 0.0;
				return false;
			} else {
				FF = strength * ( std::pow( std::cos(half_pi<double>()*(t-t0)/Ctau),int(2) ) );
				return true;
			}
		}

	private:
		double strength, Ctau, t0;

		inline bool inpulse(const double t) { return ((t-t0) >= -Ctau && (t-t0) <= Ctau); }
};

#endif
