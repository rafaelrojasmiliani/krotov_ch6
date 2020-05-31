/*!
 *       \file:  test.cpp
 *
 *		 \brief: DESCRIPTION
 *        Version:  1.0
 *        Created:  13/11/14 14:03:27
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */
#include<gsl/gsl_odeiv2.h>
#include<gsl/gsl_errno.h>
#include"eh1dofmb.h"
#include<assert.h>
//#define _GNU_SOURCE
#include <fenv.h>


int main(){
	cEh1DofMb test;
	double omega;
	double T;
	FILE *out=fopen("datahf.dat","w");

	test.Omega=1.0;
	test.setDelta0(0.0);
	test.setDeltat(0.2);
	test.findMin(1.0e-4,true,true);
}
