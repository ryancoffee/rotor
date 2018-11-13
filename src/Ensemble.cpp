#include "Ensemble.hpp"

#include <Constants.hpp> // for kb and Eh

Ensemble::Ensemble(const unsigned nstatesin = 50,const float kTinkelvinin = 300)
: nstates(nstatesin)
, kTinkelvin(kTin)
{
	setkTinau(kTinkelvin*Constants::kb<float>()/Constants::Eh<float>()):
}

Ensemble::Ensemble(unsigned nin = 50)
: nstates(nin)
{

}

void Ensemble::setkT(float kTinkelvin)
{
	setkTinau(kTinkelvin*Constants::kb<float>()/Constants::Eh<float>());
}

