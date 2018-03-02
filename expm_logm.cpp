#include <iostream>
using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////////////////
//Notice:
//  You need to include some packages to complete the two functions to calculate matrix exponent and matrix logarithm.
//	There are several softwares or freely-available packges (e.g. GNU Octave, MATLAB) calculates the two matrix functions, 
//and they can be included in c/c++ projects or can generate libs that can be included in c/c++ projects.
//  Because these packages and libs are somehow related to the operating system and runtime environment, (e.g. MATLAB lib requires MATLAB platform)
//as well as the problem of copyright, we didn't include them in our source code. 
//  Please import any packages that can calculate matrix exponent and matrix logarithm, 
//and rewrite the following two functions using these packages.
///////////////////////////////////////////////////////////////////////////////////////////////////////

double * matrix_expm(double * X, int K){
	double * Y=(double *)malloc(K*K*sizeof(double));
	cout << "Please install package for matrix exponent, see expm_logm.cpp for detail." << endl;
	system("pause");
	//input:	X is a K*K double matrix, the entry of row m column n is X[m*K+n]
	//			K indicates the size of the matrix
	//output:	Y is a K*K double matrix, the entry of row m column n is Y[m*K+n] 
	//
	//			Please complete the function to fill Y, such that Y=matrix_exponent(X), and remove the two lines above.
	//
	//Notice:	The space of Y is created by malloc here in the first line.
	//			We will free it at other place of our program.
	//			Please don't remove the "malloc" line when you remove the "cout" and "system" lines above.
	return Y;
}

double * matrix_logm(double * X, int K){
	double * Y=(double *)malloc(K*K*sizeof(double));
	cout << "Please install package for matrix logarithm, see expm_logm.cpp for detail." << endl;
	system("pause");
	//input:	X is a K*K double matrix, the entry of row m column n is X[m*K+n]
	//			K indicates the size of the matrix
	//output:	Y is a K*K double matrix, the entry of row m column n is Y[m*K+n] 
	//
	//			Please complete the function to fill Y, such that Y=matrix_logarithm(X), and remove the two lines above.
	//
	//Notice:	The space of Y is created by malloc here in the first line.
	//			We will free it at other place of our program.
	//			Please don't remove the "malloc" line when you remove the "cout" and "system" lines above.
	return Y;
}


