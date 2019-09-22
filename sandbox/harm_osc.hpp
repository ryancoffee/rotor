#ifndef HARM_OSC_H
#define HARM_OSC_H

#include <vector>
#include <iostream>
#include <fstream>
#include <complex>

typedef std::vector< double > state_type;
typedef std::vector< std::complex<double> > cstate_type;
const double gam = 0.15;

/* The rhs of x' = f(x) */
//[
void harmonic_oscillator( const state_type &x , state_type &dxdt , const double t )
{
    dxdt[0] = x[1];
    dxdt[1] = -x[0] - gam*x[1];
}
//]

//[ rhs_class
/* The rhs of x' = f(x) defined as a class */
class harm_osc {

	double m_gam;

	public:
	harm_osc( double gam ) : m_gam(gam) { }

	void operator() ( const state_type &x , state_type &dxdt , const double t )
	{
		dxdt[0] = x[1];
		dxdt[1] = -x[0] - m_gam*x[1];
	}
};

class charm_osc {

	double m_gam;

	public:
	charm_osc( double gam ) : m_gam(gam) { }

	void operator() ( const cstate_type &x , cstate_type &dxdt , const double t )
	{
		dxdt[0] = x[1];
		dxdt[1] = -x[0] - m_gam*x[1];
	}
};
//]

struct push_back_state_and_time
{
    std::vector< state_type >& m_states;
    std::vector< double >& m_times;

    push_back_state_and_time( std::vector< state_type > &states , std::vector< double > &times )
    : m_states( states ) , m_times( times ) { }

    void operator()( const state_type &x , double t )
    {
        m_states.push_back( x );
        m_times.push_back( t );
    }
};
struct push_back_cstate_and_time
{
    std::vector< cstate_type >& m_cstates;
    std::vector< double >& m_times;

    push_back_cstate_and_time( std::vector< cstate_type > &cstates , std::vector< double > &times )
    : m_cstates( cstates ) , m_times( times ) { }

    void operator()( const cstate_type &z , double t )
    {
        m_cstates.push_back( z );
        m_times.push_back( t );
    }
};

struct write_state
{
    void operator()( const state_type &x ) const
    {
        std::cout << x[0] << "\t" << x[1] << "\n";
    }
};

#endif
