// sepp_solver.cpp : Definiert den Einstiegspunkt für die Konsolenanwendung.
//
#include "stdafx.h"
#include "Nem3dCan.h"
#include "ClassHeaders.h"


int _tmain(int argc, _TCHAR* argv[])
{
//cout<<"Yes It works"<<endl;


//Stiffness matrix decleration (2 elements with 2 degrees of freedom)
csrM K;
K.CreateOfHollow(4);//#dof 


//START OF ELEMENT ASSEMBLY
//csr Form
K.writeCsrM(); 
//writing it out in the usual form 
cout<<endl; 
K.writeCsrM_arrayform(); 

int row; 
int column; 
 
//Matrix1 (bzw. element eins)
row=0; 
column=0;
K.Push_Value_incsrM(1.0,row,column); 
row=0; 
column=2;
K.Push_Value_incsrM(3.0,row,column); 
row=1; 
column=0;
K.Push_Value_incsrM(2.0,row,column); 
row=1; 
column=2;
K.Push_Value_incsrM(4.0,row,column); 

//csr Form
K.writeCsrM(); 
//writing it out in the usual form 
cout<<endl; 
K.writeCsrM_arrayform(); 

//Matrix2 (bzw. element 2)
row=0; 
column=0;
K.Push_Value_incsrM(2.0,row,column); 
row=0; 
column=3;
K.Push_Value_incsrM(4.0,row,column); 
row=3; 
column=0;
K.Push_Value_incsrM(3.0,row,column); 
row=3; 
column=3;
K.Push_Value_incsrM(6.0,row,column); 

cout<<endl;
//csr Form
K.writeCsrM(); 
//writing it out in the usual form 
cout<<endl; 
K.writeCsrM_arrayform(); 

//Matrix3 (bzw. element 3)
row=1; 
column=0;
K.Push_Value_incsrM(1.0,row,column); 
row=1; 
column=1;
K.Push_Value_incsrM(1.0,row,column); 
row=2; 
column=0;
K.Push_Value_incsrM(2.0,row,column); 
row=2; 
column=1;
K.Push_Value_incsrM(1.0,row,column); 

cout<<endl;
//csr Form
K.writeCsrM(); 
//writing it out in the usual form 
cout<<endl; 
K.writeCsrM_arrayform(); 

//Matrix4 (bzw. element 4)
row=0; 
column=2;
K.Push_Value_incsrM(3.0,row,column); 
row=0; 
column=3;
K.Push_Value_incsrM(2.0,row,column); 
row=3; 
column=2;
K.Push_Value_incsrM(3.0,row,column); 
row=3; 
column=3;
K.Push_Value_incsrM(1.0,row,column); 

cout<<endl;
//csr Form
K.writeCsrM(); 
//writing it out in the usual form 
cout<<endl; 
K.writeCsrM_arrayform(); 

//Matrix5 (bzw. element 5)
row=2; 
column=1;
K.Push_Value_incsrM(2.0,row,column); 
row=2; 
column=2;
K.Push_Value_incsrM(1.0,row,column); 
row=3; 
column=1;
K.Push_Value_incsrM(1.0,row,column); 
row=3; 
column=2;
K.Push_Value_incsrM(7.0,row,column); 

cout<<endl;
//csr Form
K.writeCsrM(); 
//writing it out in the usual form 
cout<<endl; 
K.writeCsrM_arrayform(); 

//END OF ELEMENT ASSEMBLY 

//definition and decleration of the force vector 
double* f;
double* u;
f=new double[4];
u=new double[4];

for(int i=0;i<4;i++)
u[i]=0.0; 

f[0]=-6.3;
f[1]=6.1;
f[2]=12.4;
f[3]=-4.65; 


cout<<"Start: Solving Sparse system by Pardiso unsymmetric solver"<<endl;
Pardiso_UnsymSolver(K,u,f,4,'Y');
cout<<"End: Solving Sparse system by Pardiso unsymmetric solver"<<endl;
cout<<u[0]<<endl;
cout<<u[1]<<endl;
cout<<u[2]<<endl;
cout<<u[3]<<endl;









return 0;
}

