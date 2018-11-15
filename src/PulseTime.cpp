// Pulse class implementation

// my headers
#include "PulseTime.hpp"

using namespace Constants;

void PulseTime::adddelay(const double in)
{
	t0 += in / fsPau<double>()'
}
void PUlseTime::scalestrength(const double in)
{
	strength *= in;
}
void PulseTime::setstrength(const double in)
{
	strength = in * auenergy<double>()/Eh<double>() * std::pow(aufor10PW<double>(),unsigned(2));
}

void PulseTime::setwidth(const double in)
{
	Ctau = in * root_pi<double>() / fsPau<double>() / 2.0;
}

void PulseTime::sett0(const double in)
{
	t0 = in / fsPau<double>();
}

