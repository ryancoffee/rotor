#include <boost/numeric/ublas/hermitian.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <algorithm>
#include <vector>

using namespace std::complex_literals;

typedef boost::numeric::ublas::matrix<double> boost_mat;
typedef boost::numeric::ublas::vector<double> boost_vec;
typedef boost::numeric::ublas::matrix<std::complex<double> > boost_cmat;
typedef boost::numeric::ublas::vector<std::complex<double> > boost_cvec;


double f(void){
	static double e = 0.;
	return e += 0.125;
}

int main () {
	using namespace boost::numeric::ublas;
	boost_cmat m (3, 3);
	hermitian_adaptor<matrix<std::complex<double> >, lower> hal (m);
	for (unsigned i = 0; i < hal.size1 (); ++ i) {
		for (unsigned j = 0; j < i; ++ j)
			hal (i, j) = std::complex<double> (3 * i + j, 3 * i + j);
		hal (i, i) = std::complex<double> (4 * i, 0);
	}
	std::cout << hal << std::endl;
	hermitian_adaptor<matrix<std::complex<double> >, upper> hau (m);
	for (unsigned i = 0; i < hau.size1 (); ++ i) {
		hau (i, i) = std::complex<double> (4 * i, 0);
		for (unsigned j = i + 1; j < hau.size2 (); ++ j)
			hau (i, j) = std::complex<double> (3 * i + j, 3 * i + j);
	}
	std::cout << hau << std::endl;


	boost_cvec Uvec(10);
	boost_vec ens(10);
	std::generate(ens.begin(),ens.end(),f);
	for(size_t i=0;i<ens.size();i++){
		std::cout << ens[i] << " ";
	}
	std::cout << "\n" << std:: flush;

	double dt=0.1;	
	std::transform(ens.begin()+3,ens.end(),Uvec.begin(),
			[dt](double &z) -> std::complex<double> { return std::polar(1.0,-dt*z) ;}
		      );

	std::cout << Uvec << std::endl;

	boost_cvec Pvec(10,std::complex<double>(3.0,0.));
	std::cout << Pvec << std::endl;
	std::transform(ens.begin()+2,ens.end(),Pvec.begin(),Pvec.begin(),
			[dt](double &z,std::complex<double> &p) -> std::complex<double> { return p*std::polar(1.0,-dt*z) ;}
		      );
	std::cout << Pvec << std::endl;


}

