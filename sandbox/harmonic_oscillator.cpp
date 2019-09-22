#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <boost/numeric/odeint.hpp>
#include <harm_osc.hpp>
#include <string>


int main(int argc, char* argv[]){
	std::cout << "OK, I'm here and alive\n\t ;-)\n" << std::flush;

	double start_time = 0.0;
	double stop_time = 60.0;
	double step_time = 1.1;

	harm_osc ho(0.15);
	
	if (argc > 1){
		stop_time = double(atof(argv[1]));
	}

	state_type x(2);
	x[0] = 1.;
	x[1] = 0.;
	size_t steps = boost::numeric::odeint::integrate( ho , x , start_time , stop_time , step_time );
	std::cout << "steps = " << steps << "\n" << std::flush;

	std::vector<state_type> x_vec;
	std::vector<double> times;
	size_t guessnsteps = int((stop_time-start_time)/step_time+1.);
	x_vec.reserve(guessnsteps);
	times.reserve(guessnsteps);

	// define const stepper
	boost::numeric::odeint::runge_kutta4< state_type > stepper;

	typedef boost::numeric::odeint::runge_kutta_cash_karp54< state_type > error_stepper_type;
	typedef boost::numeric::odeint::controlled_runge_kutta< error_stepper_type > controlled_stepper_type;
	controlled_stepper_type controlled_stepper;

	x[0] = 1.;
	x[1] = 0.;
	steps = boost::numeric::odeint::integrate_adaptive( controlled_stepper , ho , x , start_time , stop_time , step_time , 
			push_back_state_and_time( x_vec , times ) );

	std::string filename("./ho_output.dat");
	std::cout << "filename = " << filename << std::endl;
	std::ofstream outfile(filename.c_str(),std::ios::out);
	outfile << "#here is something?" << std::endl;
	for(size_t i; i<x_vec.size(); ++i){
	//	std::cout << times[i] << "\t" << x_vec[i][0] << "\t" << x_vec[i][1] << "\n";
		outfile << times[i] << "\t" << x_vec[i][0] << "\t" << x_vec[i][1] << "\n";
	}
	std::cout << std::flush;
	outfile.close();

        std::vector<cstate_type> z_vec;
	std::vector<double> newtimes;
        z_vec.reserve(guessnsteps);
        newtimes.reserve(guessnsteps);

        // define const stepper
        boost::numeric::odeint::runge_kutta4< cstate_type > newstepper;

        typedef boost::numeric::odeint::runge_kutta_cash_karp54< cstate_type > newerror_stepper_type;
        typedef boost::numeric::odeint::controlled_runge_kutta< newerror_stepper_type > newcontrolled_stepper_type;
        newcontrolled_stepper_type newcontrolled_stepper;


	charm_osc zho(0.15);
	cstate_type z(2);
	z[0] = 1.+ 0.j;
	z[1] = 0. + 1.j;
	steps = boost::numeric::odeint::integrate_adaptive( newcontrolled_stepper , zho , z , start_time , stop_time , step_time , 
			push_back_cstate_and_time( z_vec , newtimes ) );

	filename = "./zho_output.dat";
	std::cout << "filename = " << filename << std::endl;
	outfile.open(filename.c_str(),std::ios::out);
	outfile << "#here is something?" << std::endl;
	for(size_t i; i<z_vec.size(); ++i){
	//	std::cout << newtimes[i] << "\t" << z_vec[i][0].real() << "\t" << z_vec[i][1].real() 
	//		<< "\t" << z_vec[i][0].imag()<< "\t" << z_vec[i][1].imag() << "\n";
		outfile << newtimes[i] << "\t" << z_vec[i][0].real() << "\t" << z_vec[i][1].real() 
			<< "\t" << z_vec[i][0].imag()<< "\t" << z_vec[i][1].imag() << "\n";
	}
	std::cout << std::flush;
	outfile.close();


	return 0;
}
