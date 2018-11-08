// Pulse class implimentation

// my headers
#include "PulseTime.hpp"

void PulseTime::adddelay(const double in)
{
	t0 += in / Constants::fsPau<double>()'
}
void PUlseTime::scalestrength(const double in)
{
	strength *= in;
}
void PulseTime::setstrength(const double in)
{
	strength = in * Constants::auenergy<double>()/Constants::Eh<double>() * std::pow(Constants::aufor10PW<double>(),unsigned(2));
}

void PulseTime::setwidth(const double in)
{
	Ctau = in * Constants::root_pi<double>() / Constants::fsPau<double>() / 2.0;
}

void PulseTime::sett0(const double in)
{
	t0 = in / Constants::fsPau<double>();
}

