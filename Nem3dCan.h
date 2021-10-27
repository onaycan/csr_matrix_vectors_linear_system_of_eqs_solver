//can:boost include files for small sized linear algebraic operations 
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/storage.hpp>




//can:intel mkl lapack header files for large sized & complicated linear operations 
#include "mkl.h" 
#include "mkl_lapack.h" 
#include "mkl_blas.h" 
#include "mkl_service.h"
#include "mkl_spblas.h"
#include "mkl_pardiso.h"
#include "mkl_dss.h"
//#include "mkl_lapacke.h"
//can:standard headers 
#include <complex> 
#include <cstdlib>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <iostream>
#include <iomanip>
#include <iterator>
#include <algorithm>
#include <vector>
#include <fstream> 
#include <list>
#include <windows.h> 


#define BLACK 0
#define BLUE 1
#define GREEN 2
#define CYAN 3
#define RED 4
#define MAGENTA 5
#define BROWN 6
#define LIGHTGREY 7
#define DARKGREY 8
#define LIGHTBLUE 9
#define LIGHTGREEN 10
#define LIGHTCYAN 11
#define LIGHTRED 12
#define LIGHTMAGENTA 13
#define YELLOW 14
#define WHITE 15
#define BLINK 128



//can: definitions of canstants and namespaces 
#define PI 3.14159265
using namespace std;
namespace ublas = boost::numeric::ublas;
//#pragma warning(disable:xxxx) 
