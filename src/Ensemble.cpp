#include "Ensemble.hpp"

Ensemble::Ensemble(unsigned nin = 50)
: nstates(nin)
{

}

void Ensemble::setkT(float kTinkelvin)
{
	setkTinau(kTinkelvin*Constants::kb<float>()/Constants::Eh<float>());
}

