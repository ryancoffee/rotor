#include <iostream>
#include <vector>
#include <complex>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace std::complex_literals;

int main () {
    using namespace boost::numeric::ublas;
    boost::numeric::ublas::matrix<std::complex<double> > m (5, 3);
    vector<std::complex< double > > v (3);
    std::complex<double> z;
    for (unsigned i = 0; i < m.size2 (); ++ i) {
        for (unsigned j = 0; j < m.size1 (); ++ j){
	    z = (3 * i + j) +  1i*j;
            m (j, i) = z; 
	}
        v(i) = i + 1i;
    }

    std::cout << v << std::endl;
    std::cout << m << std::endl;
    std::cout << prod (m, v) << std::endl;
    //std::cout << prod (v, m) << std::endl; this fails if sizes mismatch
}
