#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/hermitian.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <complex>
#include <cmath>


int main () {
	using namespace boost::numeric::ublas;
	matrix<double> m (9, 9);
	banded_adaptor<matrix<double> > ba (m, 1, 1);
	for (signed i = 0; i < signed (ba.size1 ()); ++ i)
		for (signed j = std::max (i - 1, 0); j < std::min (i + 2, signed (ba.size2 ())); ++ j)
			ba (i, j) = 3 * i + j;
	std::cout << ba << std::endl;

	banded_matrix<double> bm (5,5, 1, 1);

	hermitian_adaptor<matrix<std::complex<double> >, lower> halm (m);
	for (size_t i=0;i< halm.size1();++i){
		for(size_t j= std::max(i-1,0);j<i;++j){
			halm(i,j) = std::polar(3*i + j,double(int(i-j));
		}
		halm(i,i)=3.*std::pow(double(i),int(2));
	}
	std::cout << halm << std::endl;

	/*
	Uvec *= -dt;
	std::transform(jens.begin+jens.m,jens.end,Uvec.begin,
		[dt](std::complex<double> &z) -> std::complex<double> { return std::polar(1.0,-dt*z) ;}
		);

	*/
}
