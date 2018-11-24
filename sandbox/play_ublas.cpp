#include <iostream>
#include <vector>
#include <complex>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace std::complex_literals;

int main () {
	namespace boost_ublas = boost::numeric::ublas;
	boost_ublas::matrix<std::complex<double> > m (5, 3);
	std::vector<std::complex< double > > v (3,10+20i);

	std::complex<double> z;
	for (unsigned i = 0; i < m.size2 (); ++ i) {
		for (unsigned j = 0; j < m.size1 (); ++ j){
			z = (3 * i + j) +  1i*j;
			m (j, i) = z; 
		}
		v[i] = i + 1i;
	}
	boost_ublas::vector<std::complex< double> > vin(v.size(),*v.data());
	std::cout << vin << std::endl;
	std::cout << m << std::endl;

	std::cout << "\n\n" << std::endl;


	auto temp = boost_ublas::prod(m,vin);


	std::cout << temp << std::endl;
	std::cout << "\n\n\tv is intermediately:\n";
	v.resize(5,-1-5i);
	for (size_t i = 0; i<v.size();++i)
		std::cout << v[i] << "\t";
	std::cout << std::endl;

	// How to copy a boost blas vec to a standard vector //
	std::copy(temp.begin(),temp.end(),v.begin());
	std::cout << "\n\n\tNow v is:\n";
	for (size_t i = 0; i<v.size();++i)
		std::cout << v[i] << "\t";
	std::cout << std::endl;

}
