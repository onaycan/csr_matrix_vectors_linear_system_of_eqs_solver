#include "stdafx.h"
#include "Nem3dCan.h"
#include "ClassHeaders.h"


double Identity(int i,int j)
{
double result; 
if(i==j)
result=1.0; 
else 
result=0.0; 
return result; 
}

void Matrixdet(double M[][3], double& detM)
{
detM=M[0][0]*(M[1][1]*M[2][2]-M[2][1]*M[1][2])+
     M[0][1]*(M[2][0]*M[1][2]-M[1][0]*M[2][2])+
     M[0][2]*(M[1][0]*M[2][1]-M[2][0]*M[1][1]);
}

void Matrixdet(ublas::c_matrix<double,3,3> M, double& detM)
{
detM=M(0,0)*(M(1,1)*M(2,2)-M(2,1)*M(1,2))+
     M(0,1)*(M(2,0)*M(1,2)-M(1,0)*M(2,2))+
     M(0,2)*(M(1,0)*M(2,1)-M(2,0)*M(1,1));
}


int MatrixInv(double M[][3], double Minv[][3])
{
int success=0; 
double detM; 
Minv[0][0]=M[1][1]*M[2][2]-M[2][1]*M[1][2];
Minv[0][1]=M[0][2]*M[2][1]-M[2][2]*M[0][1];
Minv[0][2]=M[0][1]*M[1][2]-M[1][1]*M[0][2];


Minv[1][0]=M[1][2]*M[2][0]-M[2][2]*M[1][0];
Minv[1][1]=M[0][0]*M[2][2]-M[2][0]*M[0][2];
Minv[1][2]=M[0][2]*M[1][0]-M[1][2]*M[0][0];


Minv[2][0]=M[1][0]*M[2][1]-M[2][0]*M[1][1];
Minv[2][1]=M[0][1]*M[2][0]-M[2][1]*M[0][0];
Minv[2][2]=M[0][0]*M[1][1]-M[1][0]*M[0][1]; 
Matrixdet(M,detM);
if(detM!=0)
success=1; 
for(int r=0;r<3;r++)
for(int c=0;c<3;c++)
Minv[r][c]=Minv[r][c]/detM;  
return success; 
}


ublas::matrix<double> MatrixInv(ublas::c_matrix<double,3,3> M, double &_detM)
{
ublas::c_matrix<double,3,3> Minv; 
Minv(0,0)=M(1,1)*M(2,2)-M(2,1)*M(1,2);
Minv(0,1)=M(0,2)*M(2,1)-M(2,2)*M(0,1);
Minv(0,2)=M(0,1)*M(1,2)-M(1,1)*M(0,2);


Minv(1,0)=M(1,2)*M(2,0)-M(2,2)*M(1,0);
Minv(1,1)=M(0,0)*M(2,2)-M(2,0)*M(0,2);
Minv(1,2)=M(0,2)*M(1,0)-M(1,2)*M(0,0);

Minv(2,0)=M(1,0)*M(2,1)-M(2,0)*M(1,1);
Minv(2,1)=M(0,1)*M(2,0)-M(2,1)*M(0,0);
Minv(2,2)=M(0,0)*M(1,1)-M(1,0)*M(0,1); 

Matrixdet(M,_detM);
Minv=Minv/_detM;
return Minv; 
}

//can: returns the length of a ublas vector
double len(ublas::c_vector<double,3> v)
{
return sqrt(inner_prod(v,v));
}

double len4(ublas::c_vector<double,4> v)
{
return sqrt(prec_inner_prod(v,v));
}

void MplusM(double res[], double M1[], double M2[], int _m)
{
for(int i=0;i<_m*_m;i++)
res[i]=0; 
for(int i=0;i<_m*_m;i++)
res[i]=M1[i]+M2[i];

}


//can:returns the cross product of two ublas vector in 3d 
ublas::c_vector<double,3> cross(ublas::c_vector<double,3> v1,ublas::c_vector<double,3> v2)
{
ublas::c_vector<double,3> result; 

result[0]=v1[1]*v2[2]-v2[1]*v1[2];
result[1]=-v1[0]*v2[2]+v2[0]*v1[2];
result[2]=v1[0]*v2[1]-v2[0]*v1[1];

return result; 
}


int sign(int value)
{
int result=0;
if(value>0)
result=1;
if(value<0)
result=-1; 
return result;
}

void mkl_lapack_InvMatrix(double _M[], int _r, int _c)
{

double *work=new double[_r*_c]; 
int m=_r;
int n=_c;
int lda=_r;
int* ipiv=new int[_r]; 
int lwork=_r*_c; 
int info=1; 



dgetrf(&m,&n,_M,&lda,ipiv,&info); 
cout<<"info of LU factorization in inverting: "<<info<<endl; 


dgetri(&n,_M,&lda,ipiv,work,&lwork,&info); 
cout<<"info of inversion: "<<info<<endl; 


delete[] work;
delete[] ipiv; 


}

void mkl_blas_Mdotv(double _M[], double _left[], double _right[], int _dim)
{
char trans='N';
int lda=_dim;
int m=_dim;
int n=_dim;
int incx=1;
int incy=1;
double alpha=1.0;
double beta=0.0; 
dgemv(&trans,&m,&n,&alpha,_M,&lda,_left,&incx,&beta,_right,&incy);

}

//Matrix dot matrix by using mkl_lapack routines 
void mkl_blas_MdotM(double res[], double M1[], double M2[], int _r1, int _c1r2, int _c2)
{

//inputs for blas 
    //change it 
char transa='N'; 
char transb='N'; 
int m=_r1; 
int n=_c2; 
int k=_c1r2; 
double alpha=1.0; 
int lda=m; 
int ldb=k;
double beta=0.0; 
int ldc=m; 



dgemm(&transa,&transb,&m,&n,&k,&alpha,M1,&lda,M2,&ldb,&beta,res,&ldc);


}


void mkl_lapack_Solve(double _K[], double _lhs[], double _rhs[], int _size)
{

double *work=new double[_size*_size]; 
double* _rhstemp=new double[_size];
int m=_size; 
int n=_size; 
int nrhs=1;  
int lda=_size; 
int ldb=_size; 
int* ipiv=new int[_size]; 
int lwork=_size*_size; 
int info=1; 

char trans='N'; 


dgetrf(&m,&n,_K,&lda,ipiv,&info); 
cout<<"info of LU factorization in Solver: "<<info<<endl; 

for(int i=0;i<_size;i++)
_rhstemp[i]=_rhs[i];

dgetrs(&trans,&m,&nrhs,_K,&lda,ipiv,_rhstemp,&ldb,&info); 
cout<<"info of Solver: "<<info<<endl; 

for(int i=0;i<_size;i++)
_lhs[i]=_rhstemp[i];

delete[] work; 
delete[] _rhstemp; 
delete[] ipiv; 
}


//PARDISOPARDISOPARDISO
void mkl_lapack_SparseSolve(csrM _K, double _lhs[], double _rhs[], int _size)
{
//can: this might be a pointer as well 
int pt[64];
for(int i=0;i<64;i++)
pt[i]=0;

int maxfct=1;
int mnum=0; 
int type=-2; 
int phase; 
int n=_K.iM.size()-1; 
cout<<"n:"<<n<<endl;
double* a=new double[_K.vM.size()];
for(int i=0;i<_K.vM.size();i++)
a[i]=_K.vM[i];

int*ia=new int[_K.iM.size()];
for(int i=0;i<_K.iM.size();i++)
ia[i]=_K.iM[i]+1;

int* ja=new int[_K.jM.size()];
for(int i=0;i<_K.jM.size();i++)
ja[i]=_K.jM[i]+1;

int*perm=new int[n];

int nrhs=1; 
int iparm[64];
iparm[0]=0;
int msglvl=1; 

//lhs=x
//rhs=b
_MKL_DSS_HANDLE_t handle;
int error; 
pardiso(&handle, &maxfct, &mnum, &type, &phase, &n, a, ia, ja, perm, &nrhs, iparm, &msglvl, _rhs, _lhs, &error);
cout<<"pardiso ended actually:"<<endl;
/*
char transa='N';
int _m=size; 
double alpha=1.0;
char* matdescra=new char[6]; 

double* val=new double[_K.vM.size()];
for(int i=0;i<_K.vM.size();i++)
val[i]=_K.vM[i]; 
int* indx=new int[_K.jM.size()];
for(int i=0;i<_K.jM.size();i++)
indx[i]=_K.jM[i]; 
int* pntrb=new int[_K.iM.size()-1]; 
for(int i=0;i<_K.iM.size()-1;i++)
pntrb[i]=_K.iM[i];
int* pntre=new int[_K.iM.size()-1]; 
for(int i=1;i<_K.iM.size();i++)
pntre[i-1]=_K.iM[i];

double *work=new double[_size*_size]; 
double* _rhstemp=new double[_size];
int m=_size; 
int n=_size; 
int nrhs=1;  
int lda=_size; 
int ldb=_size; 
int* ipiv=new int[_size]; 
int lwork=_size*_size; 
int info=1; 

char trans='N'; 


dgetrf(&m,&n,_K,&lda,ipiv,&info); 
cout<<"info of LU factorization in Solver: "<<info<<endl; 

for(int i=0;i<_size;i++)
_rhstemp[i]=_rhs[i];

dgetrs(&trans,&m,&nrhs,_K,&lda,ipiv,_rhstemp,&ldb,&info); 
cout<<"info of Solver: "<<info<<endl; 

for(int i=0;i<_size;i++)
_lhs[i]=_rhstemp[i];

delete[] work; 
delete[] _rhstemp; 
delete[] ipiv; 
*/ 
}


double ArrayNorm(double _V[], int _size)
{
double result=0.0; 
for(int i=0;i<_size;i++)
result=result+_V[i]*_V[i]; 
return sqrt(result); 
}


void Decimal_to_thertary(int decimal, int thertary[])
{
int remainder=decimal%3; 
int division=decimal/3; 
thertary[3]=remainder; 
remainder=division%3; 
decimal=division/3; 
thertary[2]=remainder; 
remainder=division%3; 
division=decimal/3; 
thertary[1]=remainder; 
remainder=division%3; 
decimal=division/3; 
thertary[0]=remainder;
}

void Set_VectorZero(double _v[], int _size)
{
for(int i=0;i<_size;i++)
_v[i]=0.0;
}


void Neumann_Inverse(double*& _M,int*& _jM, int*& _iM, int _m)
{
char trans='N'; 
int request;
int sort; 
double beta;
int nzmax; 
int info; 

//can: Identity matrix m_by_m in csr format one indexing 
double* I=new double[_m];
int* jI=new int[_m]; 
int* iI=new int[_m+1];
for(int i=0;i<_m;i++)
{
I[i]=1.0; 
jI[i]=i+1;
iI[i]=i+1;
}
iI[_m]=_m+1;

//can: identity matrix minus the matrix to be inverted...
//can: this is the matrix whose powers to be calculated...
double* IminM; 
int* jIminM; 
int* iIminM=new int[_m+1];
request=1; 
beta=-1.0;
mkl_dcsradd(&trans, &request, &sort, &_m, &_m,
            I, jI, iI,
            &beta,
            _M, _jM, _iM,
            IminM, jIminM, iIminM,
            &nzmax, &info);

IminM=new double[iIminM[_m]-1];
jIminM=new int[iIminM[_m]-1];
request=2;
info=5; 
beta=-1.0;
mkl_dcsradd(&trans, &request, &sort, &_m, &_m,
            I, jI, iI,
            &beta,
            _M, _jM, _iM,
            IminM, jIminM, iIminM,
            &nzmax, &info);

cout<<"IminM"<<endl;
for(int i=0;i<iIminM[_m]-1;i++)
cout<<IminM[i]<<'\t';
cout<<endl; 
for(int i=0;i<iIminM[_m]-1;i++)
cout<<jIminM[i]<<'\t';
cout<<endl; 
for(int i=0;i<_m+1;i++)
cout<<iIminM[i]<<'\t';
cout<<endl; 
 
//can: the ver first previous matrix is equal to the identity matrix
double* Mprev=new double[_m];
int* jMprev=new int[_m]; 
int* iMprev=new int[_m+1];
for(int i=0;i<_m;i++)
{
Mprev[i]=1.0; 
jMprev[i]=i+1;
iMprev[i]=i+1;
}
iMprev[_m]=_m+1;
//can: we do not have any information on the current matrix 
double* Mcur; 
int* jMcur; 
int* iMcur=new int[_m+1]; 
for(int j=0;j<50;j++)
{
request=1; 
sort=7;
mkl_dcsrmultcsr(&trans, &request, &sort,
                &_m, &_m, &_m,
                IminM, jIminM, iIminM,
                Mprev, jMprev, iMprev,
                Mcur, jMcur, iMcur,
                &nzmax, &info);
Mcur=new double[iMcur[_m]-1];
jMcur=new int[iMcur[_m]-1];
request=2;
info=5; 
mkl_dcsrmultcsr(&trans, &request, &sort,
                &_m, &_m, &_m,
                IminM, jIminM, iIminM,
                Mprev, jMprev, iMprev,
                Mcur, jMcur, iMcur,
                &nzmax, &info);

delete[] Mprev; 
delete[] jMprev; 
delete[] iMprev; 
iMprev=new int[_m+1];
request=1; 
beta=1.0;
mkl_dcsradd(&trans, &request, &sort, &_m, &_m,
            I, jI, iI,
            &beta,
            Mcur, jMcur, iMcur,
            Mprev, jMprev, iMprev,
            &nzmax, &info);
Mprev=new double[iMprev[_m]-1];
jMprev=new int[iMprev[_m]-1];
request=2;
info=5; 
mkl_dcsradd(&trans, &request, &sort, &_m, &_m,
            I, jI, iI,
            &beta,
            Mcur, jMcur, iMcur,
            Mprev, jMprev, iMprev,
            &nzmax, &info);

delete[] Mcur; 
delete[] jMcur; 
delete[] iMcur; 
iMcur=new int[_m+1]; 

}

delete[] I;
delete[] jI; 
delete[] iI; 

delete[] IminM; 
delete[] jIminM; 
delete[] iIminM; 


delete[] _M; 
delete[] _jM; 
delete[] _iM; 

_M=new double[iMprev[_m]-1];
_jM=new int[iMprev[_m]-1]; 
_iM=new int[_m+1];

for(int i=0;i<iMprev[_m]-1;i++)
{
_M[i]=Mprev[i];
_jM[i]=jMprev[i]; 
}

for(int i=0;i<_m+1;i++)
_iM[i]=iMprev[i]; 

delete[] Mprev; 
delete[] jMprev; 
delete[] iMprev; 

}
/*

void Neumann_IminM_Inverse(csrM &_M, int _m)
{
char trans='N'; 
int request;
int sort; 
double beta;
int nzmax;
int info; 




double* I=new double[_m];
int* jI=new int[_m]; 
int* iI=new int[_m+1];
for(int i=0;i<_m;i++)
{
I[i]=1.0; 
jI[i]=i+1;
iI[i]=i+1;
}
iI[_m]=_m+1;

//can: the ver first previous matrix is equal to the identity matrix
double* Mprev=new double[_m];
int* jMprev=new int[_m]; 
int* iMprev=new int[_m+1];
for(int i=0;i<_m;i++)
{
Mprev[i]=1.0; 
jMprev[i]=i+1;
iMprev[i]=i+1;
}
iMprev[_m]=_m+1;
//can: we do not have any information on the current matrix 
double* Mcur; 
int* jMcur; 
int* iMcur=new int[_m+1]; 
for(int j=0;j<50;j++)
{
request=1; 
sort=7;
mkl_dcsrmultcsr(&trans, &request, &sort,
                &_m, &_m, &_m,
                _M.Mval, _M.jMval, _M.iMval,
                Mprev, jMprev, iMprev,
                Mcur, jMcur, iMcur,
                &nzmax, &info);

Mcur=new double[iMcur[_m]-1];
jMcur=new int[iMcur[_m]-1];
request=2;
info=5; 
mkl_dcsrmultcsr(&trans, &request, &sort,
                &_m, &_m, &_m,
                _M.Mval, _M.jMval, _M.iMval,
                Mprev, jMprev, iMprev,
                Mcur, jMcur, iMcur,
                &nzmax, &info);

delete[] Mprev; 
delete[] jMprev; 
delete[] iMprev; 
iMprev=new int[_m+1];
request=1; 
beta=1.0;
mkl_dcsradd(&trans, &request, &sort, &_m, &_m,
            I, jI, iI,
            &beta,
            Mcur, jMcur, iMcur,
            Mprev, jMprev, iMprev,
            &nzmax, &info);
Mprev=new double[iMprev[_m]-1];
jMprev=new int[iMprev[_m]-1];
request=2;
info=5; 
mkl_dcsradd(&trans, &request, &sort, &_m, &_m,
            I, jI, iI,
            &beta,
            Mcur, jMcur, iMcur,
            Mprev, jMprev, iMprev,
            &nzmax, &info);

delete[] Mcur; 
delete[] jMcur; 
delete[] iMcur; 
iMcur=new int[_m+1]; 

}

delete[] I;
delete[] jI; 
delete[] iI; 


delete[] _M.Mval; 
delete[] _M.jMval; 
delete[] _M.iMval; 

_M.Mval=new double[iMprev[_m]-1];
_M.jMval=new int[iMprev[_m]-1]; 
_M.iMval=new int[_m+1];

for(int i=0;i<iMprev[_m]-1;i++)
{
_M.Mval[i]=Mprev[i];
_M.jMval[i]=jMprev[i]; 
}

for(int i=0;i<_m+1;i++)
_M.iMval[i]=iMprev[i]; 

delete[] Mprev; 
delete[] jMprev; 
delete[] iMprev; 

}


*/ 

void Csr_add(double*& _M,int*& _jM, int*& _iM,
             double*& _N,int*& _jN, int*& _iN, 
             double*& _MandN,int*& _jMandN, int*& _iMandN, 
             int _m)
{
char trans='N'; 
int request=1;
int sort; 
double beta=1.0;
int nzmax;
delete[] _iMandN;
_iMandN=new int[_m+1]; 
int info; 


mkl_dcsradd(&trans, &request, &sort, &_m, &_m,
            _M, _jM, _iM,
            &beta,
            _N, _jN, _iN,
            _MandN, _jMandN, _iMandN,
            &nzmax, &info);

_MandN=new double[_iMandN[_m]-1];
_jMandN=new int[_iMandN[_m]-1];
request=2;
info=5; 

mkl_dcsradd(&trans, &request, &sort, &_m, &_m,
            _M, _jM, _iM,
            &beta,
            _N, _jN, _iN,
            _MandN, _jMandN, _iMandN,
            &nzmax, &info);
}


void Csr_add_and_assign(double*& _M,int*& _jM, int*& _iM,
                        double*& _N,int*& _jN, int*& _iN, 
                        int _m)
{
char trans='N'; 
int request=1;
int sort; 
double beta=1.0;
int nzmax;
double* MandN; 
int* jMandN;
int* iMandN;
int info; 


Csr_add(_M,_jM,_iM,_N,_jN,_iN,MandN,jMandN,iMandN,_m);


delete []_M; 
delete []_jM; 
delete []_iM; 

_M=new double[iMandN[_m]-1];
_jM=new int[iMandN[_m]-1];
_iM=new int[_m+1];


for(int i=0;i<iMandN[_m]-1;i++)
{
_M[i]=MandN[i];
_jM[i]=jMandN[i];
}


for(int i=0;i<_m+1;i++)
{
_iM[i]=iMandN[i];
}

delete []MandN; 
delete []jMandN; 
delete []iMandN; 

}


void Csr_subs(double*& _M,int*& _jM, int*& _iM,
              double*& _N,int*& _jN, int*& _iN, 
              double*& _MorN,int*& _jMorN, int*& _iMorN, 
              int _m)
{
char trans='N'; 
int request=1;
int sort; 
double beta=-1.0;
int nzmax;
delete[] _iMorN;
_iMorN=new int[_m+1]; 
int info; 


mkl_dcsradd(&trans, &request, &sort, &_m, &_m,
            _M, _jM, _iM,
            &beta,
            _N, _jN, _iN,
            _MorN, _jMorN, _iMorN,
            &nzmax, &info);

_MorN=new double[_iMorN[_m]-1];
_jMorN=new int[_iMorN[_m]-1];
request=2;
info=5; 

mkl_dcsradd(&trans, &request, &sort, &_m, &_m,
            _M, _jM, _iM,
            &beta,
            _N, _jN, _iN,
            _MorN, _jMorN, _iMorN,
            &nzmax, &info);
}


void Csr_mult(double*& _M,int*& _jM, int*& _iM,
              double*& _N,int*& _jN, int*& _iN, 
              double*& _MthenN,int*& _jMthenN, int*& _iMthenN, 
              int _m)
{

char trans='N'; 
int request=1;
int sort=1; 
int nzmax;
_iMthenN=new int[_m+1]; 
int info; 


mkl_dcsrmultcsr(&trans, &request, &sort,
                &_m, &_m, &_m,
                _M, _jM, _iM,
                _N, _jN, _iN,
                _MthenN, _jMthenN, _iMthenN,
                &nzmax, &info);

_MthenN=new double[_iMthenN[_m]-1];
_jMthenN=new int[_iMthenN[_m]-1];
request=2;
info=5; 

mkl_dcsrmultcsr(&trans, &request, &sort,
                &_m, &_m, &_m,
                _M, _jM, _iM,
                _N, _jN, _iN,
                _MthenN, _jMthenN, _iMthenN,
                &nzmax, &info);


}






int Pardiso_UnsymSolver(csrM _K, double _lhs[], double _rhs[], int _size, char comment)
{

int n=_K.iM.size()-1; 
//cout<<"n:"<<n<<endl;
double* a=new double[_K.vM.size()];
//cout<<"a:"<<endl;
for(int i=0;i<_K.vM.size();i++)
{
a[i]=_K.vM[i];
//cout<<a[i]<<'\t';
}
//cout<<endl;
//cout<<"ia:"<<endl;
int*ia=new int[_K.iM.size()];
for(int i=0;i<_K.iM.size();i++)
{
ia[i]=_K.iM[i]+1;
//cout<<ia[i]<<'\t';
}
//cout<<endl;
//cout<<"ja:"<<endl;
int* ja=new int[_K.jM.size()];
for(int i=0;i<_K.jM.size();i++)
{
ja[i]=_K.jM[i]+1;
//cout<<ja[i]<<'\t';
}
//cout<<endl;

/*
int n = 8;
int ia[ 9] = { 1, 5, 8, 10, 12, 15, 17, 18, 19 };
int ja[18] = { 1, 3, 6, 7,
2, 3, 5,
3, 8,
4, 7,
5, 6, 7,
6, 8,
7,
8 };
double a[18] = { 7.0, 1.0, 2.0, 7.0,
-4.0, 8.0, 2.0,
1.0, 5.0,
7.0, 9.0,
5.0, 1.0, 5.0,
-1.0, 5.0,
11.0,
5.0 };
double  b[8], x[8];
*/
int mtype = 11; /* Real unsymmetric matrix */
/* RHS and solution vectors.*/
int nrhs = 1; /* Number of right hand sides. */
/* Internal solver memory pointer pt, */
/* 32-bit: int pt[64]; 64-bit: long int pt[64] */
/* or void *pt[64] should be OK on both architectures */
void *pt[64];
/* Pardiso control parameters.*/
int iparm[64];
int maxfct, mnum, phase, error, msglvl;
/* Auxiliary variables. */
int i;
double ddum; /* Double dummy*/
int idum; /* Integer dummy.*/
/* --------------------------------------------------------------------*/
/* .. Setup Pardiso control parameters.*/
/* --------------------------------------------------------------------*/
for (i = 0; i < 64; i++) {
iparm[i] = 0;
}
iparm[0] = 1; /* No solver default*/
iparm[1] = 2; /* Fill-in reordering from METIS */
/* Numbers of processors, value of MKL_NUM_THREADS */
iparm[2] = mkl_get_max_threads();
iparm[3] = 0; /* No iterative-direct algorithm */
iparm[4] = 0; /* No user fill-in reducing permutation */
iparm[5] = 0; /* Write solution into x */
iparm[6] = 16; /* Default logical fortran unit number for output */
iparm[7] = 2; /* Max numbers of iterative refinement steps */
iparm[8] = 0; /* Not in use*/
iparm[9] = 13; /* Perturb the pivot elements with 1E-13 */
iparm[10] = 1; /* Use nonsymmetric permutation and scaling MPS */
iparm[11] = 0; /* Not in use*/
iparm[12] = 0; /* Not in use*/
iparm[13] = 0; /* Output: Number of perturbed pivots */
iparm[14] = 0; /* Not in use*/
iparm[15] = 0; /* Not in use*/
iparm[16] = 0; /* Not in use*/
iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
iparm[18] = -1; /* Output: Mflops for LU factorization */
iparm[19] = 0; /* Output: Numbers of CG Iterations */
maxfct = 1; /* Maximum number of numerical factorizations. */
mnum = 1; /* Which factorization to use. */
msglvl = 0; /* Don't print statistical information in file */
error = 0; /* Initialize error flag */
/* --------------------------------------------------------------------*/
/* .. Initialize the internal solver memory pointer. This is only */
/* necessary for the FIRST call of the PARDISO solver. */
/* --------------------------------------------------------------------*/
for (i = 0; i < 64; i++) {
pt[i] = 0;
}
/* --------------------------------------------------------------------*/
/* .. Reordering and Symbolic Factorization. This step also allocates */
/* all memory that is necessary for the factorization. */
/* --------------------------------------------------------------------*/
phase = 11;
PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
&n, a, ia, ja, &idum, &nrhs,
iparm, &msglvl, &ddum, &ddum, &error);
if (error != 0) {
printf("\nERROR during symbolic factorization: %d", error);
exit(1);
}
if(comment=='Y')
{
printf("\nReordering completed ... ");
printf("\nNumber of nonzeros in factors = %d", iparm[17]);
printf("\nNumber of factorization MFLOPS = %d", iparm[18]);
}
/* --------------------------------------------------------------------*/
/* .. Numerical factorization.*/
/* --------------------------------------------------------------------*/
phase = 22;
PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
&n, a, ia, ja, &idum, &nrhs,
iparm, &msglvl, &ddum, &ddum, &error);
if (error != 0) {
printf("\nERROR during numerical factorization: %d", error);
exit(2);
}
if(comment=='Y')
printf("\nFactorization completed ... ");
/* --------------------------------------------------------------------*/
/* .. Back substitution and iterative refinement. */
/* --------------------------------------------------------------------*/
phase = 33;
iparm[7] = 2; /* Max numbers of iterative refinement steps. */
/* Set right hand side to one.*/

/*
for (i = 0; i < n; i++) {
b[i] = 1;
}
*/ 
PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
&n, a, ia, ja, &idum, &nrhs,
iparm, &msglvl, _rhs, _lhs, &error);
if (error != 0) {
printf("\nERROR during solution: %d", error);
exit(3);
}
if(comment=='Y')
printf("\nSolve completed ... ");


/*
printf("\nThe solution of the system is: ");
for (i = 0; i < n; i++) {
printf("\n x [%d] = % f", i, _lhs[i] );
}
printf ("\n");
*/ 
/* --------------------------------------------------------------------*/
/* .. Termination and release of memory. */
/* --------------------------------------------------------------------*/
phase = -1; /* Release internal memory. */
PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
&n, &ddum, ia, ja, &idum, &nrhs,
iparm, &msglvl, &ddum, &ddum, &error);
delete []a;
delete []ia;
delete []ja;
return 0;
cout<<endl;
}


void SetArrayZero(double _v[], int _m)
{
for(int i=0;i<_m;i++)
_v[i]=0.0; 
}


void MatrixDotVector(double _M[], double _V[], double  _Result[], int _size)
{
for(int r=0;r<_size;r++)
_Result[r]=0.0; 

for(int i=0;i<_size;i++)
for(int j=0;j<_size;j++)
_Result[i]=_Result[i]+_M[i+j*_size]*_V[j]; 

}


void SpaceAngleBetween(ublas::c_vector<double,3> v1,
                       ublas::c_vector<double,3> v2,
                       double &absangle)
{
double cosangle=inner_prod(v1/norm_2(v1),v2/norm_2(v2));
if(cosangle>1.0)
cosangle=1.0;
if(cosangle<-1.0)
cosangle=-1.0;
absangle=acos(cosangle)*180.0/PI;
//absangle=acos(inner_prod(v1/norm_2(v1),v2/norm_2(v2)))*180.0/PI;
//if(absangle>90)
//absangle=absangle-90;
//absangle=len(v2);
}




void AbsAngleBetween(ublas::c_vector<double,3> v1,
                     ublas::c_vector<double,3> v2,
                     double &absangle)
{

double cosangle=inner_prod(v1/norm_2(v1),v2/norm_2(v2));
if(cosangle>1.0)//can: sounds weird but precision gives some value greater then 1 
cosangle=1.0;
if(cosangle<-1.0)
cosangle=-1.0;
absangle=acos(cosangle)*180.0/PI;

if(absangle>90)
absangle=180-absangle;


/*
//added later on 
if(absangle<-90)
absangle=180-fabs(absangle);

absangle=fabs(absangle);
*/ 
//absangle=len(v2);
}


void RotateVectorAroundAxis(ublas::c_vector<double,3> &rotated, 
                            ublas::c_vector<double,3> axis,
                            double theta)
{
ublas::c_vector<double,3> unrotated; 
for(int i=0;i<3;i++)
unrotated[i]=rotated[i]; 

//R[i][j]=cronecker[i][j]cos(theta)+K[i][j]sin(theta)+(1-cos(theta))kk[i][j];
ublas::c_matrix<double,3,3> K; 
K(0,0)=0.0;
K(1,1)=0.0;
K(2,2)=0.0; 
K(0,1)=-axis[2];
K(1,0)=axis[2];
K(0,2)=axis[1];
K(2,0)=-axis[1];
K(1,2)=-axis[0];
K(2,1)=axis[0];

ublas::c_matrix<double,3,3> kk; 
for(int i=0;i<3;i++)
for(int j=0;j<3;j++)
kk(i,j)=axis[i]*axis[j];

theta=theta*PI/180.0;
ublas::c_matrix<double,3,3> R;
for(int i=0;i<3;i++)
for(int j=0;j<3;j++)
R(i,j)=Identity(i,j)*cos(theta)+K(i,j)*sin(theta)+kk(i,j)*(1.0-cos(theta));

rotated=prod(R,unrotated);

}