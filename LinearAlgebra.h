

//can: cronecker
double Identity(int i,int j);
//can::the determinant of a matrix of size 3 in static array form 
void Matrixdet(double M[][3], double& detM);
//can:determinant of a 3*3 matrix in ublas form 
void Matrixdet(ublas::c_matrix<double,3,3> M, double& detM);
//can:inverse of a 3*3 matrix in array form 
int MatrixInv(double M[][3], double Minv[][3]);
//can: matrix time vector in array format 
void MatrixDotVector(double _M[], double _V[], double  _Result[], int _size);
//can:inversion of a 3*3 matrix in ublas form 
ublas::matrix<double> MatrixInv(ublas::c_matrix<double,3,3> M, double &_detM);
//can:length of a ublas vector 
double len(ublas::c_vector<double,3> v);
double len4(ublas::c_vector<double,4> v);
//can:returns the cross product of two ublas vector in 3d 
ublas::c_vector<double,3> cross(ublas::c_vector<double,3> v1,ublas::c_vector<double,3> v2);
//can: returns the sign of an integer 
int sign(int value);
//can: mkl explicit inversion of a matrix 
void mkl_lapack_InvMatrix(double _M[], int _r, int _c);
//can: mkl matrix dot vector 
void mkl_blas_Mdotv(double _M[], double _left[], double _right[], int _dim);
//can: matrix times matrix ftom mkl blas 
void mkl_blas_MdotM(double res[], double M1[], double M2[], int _r1, int _c1r2, int _c2);
//can: solvers 
void mkl_lapack_Solve(double _K[], double _rhs[], double _lhs[], int _size);
void mkl_lapack_SparseSolve(csrM _K, double _lhs[], double _rhs[], int _size);
//can:setting an array zero 
void Set_VectorZero(double _v[], int _size);
//can: matrix plus matrix in naive array form 
void MplusM(double res[], double M1[], double M2[], int _m);
//can: Iterative explicit inversion of a square sparse matrix in csr form. 
void Neumann_Inverse(double*& _M,int*& _jM, int*& _iM, int _m);
//void Neumann_IminM_Inverse(csrM &_M, int _m);
//can:addition of two csr matrices 
void Csr_add(double*& _M,int* &_jM, int*& _iM,
             double*& _N,int* &_jN, int*& _iN, 
             double*& _MandN,int*& _jMandN, int*& _iMandN, 
             int _m); 
//can:addition of two csr matrices, the result is assigned on to the first matrix 
void Csr_add_and_assign(double*& _M,int*& _jM, int*& _iM,
                        double*& _N,int*& _jN, int*& _iN, 
                        int _m);
//can: substraction of two csr matrices of m by m 
void Csr_subs(double*& _M,int*& _jM, int*& _iM,
              double*& _N,int*& _jN, int*& _iN, 
              double*& _MorN,int*& _jMorN, int*& _iMorN, 
              int _m);
//can: multiplication of two square m by m sparse matrices of form csr  
void Csr_mult(double*& _M,int*& _jM, int*& _iM,
              double*& _N,int*& _jN, int*& _iN, 
              double*& _MthenN,int*& _jMthenN, int*& _iMthenN, 
              int _m); 
//can: sets an array zero 
void SetArrayZero(double _v[], int _m);
extern int omp_get_max_threads();
/* PARDISO prototype. */
/*
extern void PARDISO
(void *, int *, int *, int *, int *, int *,
 double a[], int ia[], int ja[], int*, int *, int *,
int *, double b[], double x[], int*);
*/ 
//can:sparse matrix solver for unsymmetric case 
int Pardiso_UnsymSolver(csrM _K, double _lhs[], double _rhs[], int _size, char comment);
double ArrayNorm(double _V[], int _size);
void Decimal_to_thertary(int decimal, int thertary[]);
void SpaceAngleBetween(ublas::c_vector<double,3> v1,
                     ublas::c_vector<double,3> v2,
                     double &absangle);
void AbsAngleBetween(ublas::c_vector<double,3> v1,
                     ublas::c_vector<double,3> v2,
                     double &absangle);

void RotateVectorAroundAxis(ublas::c_vector<double,3> &rotated, 
                            ublas::c_vector<double,3> axis,
                            double theta); 
