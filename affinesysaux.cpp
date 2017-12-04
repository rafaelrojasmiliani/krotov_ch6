/*!
 *       \file:  affinesysaux.cpp
 *
 *		 \brief: DESCRIPTION
 *        Version:  1.0
 *        Created:  03/08/15 20:21:29
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */



#include<eigen3/Eigen/Core>
#include"indexvector.h"
#include"affinesys.h"
#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
using namespace Eigen;

void affinesys::dynSysLinearPartFiniteDiffTest(double t,const double xin[]){
	Map<const VectorXd> x0(xin,stateDim);
	typedef Matrix<double,Dynamic,Dynamic,RowMajor> MatrixXX;
	MatrixXX hi(stateDim,affineDim);
	MatrixXX dhdxError(stateDim,affineDim);
	MatrixXX dhdxiFD(stateDim,affineDim);
	double buffdhdx[500];
	Map<MatrixXX> dhdxi(NULL,0,0);
	VectorXd x(stateDim);

	static const double coeff[9]={	 1.0/280.0	,-4.0/105.0	, 1.0/5.0,
									-4.0/5.0	, 0.0		, 4.0/5.0, 	
									-1.0/5.0	, 4.0/105.0	,-1.0/280.0};
	static size_t uici,uicj;
	const double dx=0.01;
	dynSysLinearPart_1G(t,xin,buffdhdx);


	const size_t stride=stateDim*affineDim;
	for(uici=0;uici<stateDim;uici++){
		new (&dhdxi) Map< MatrixXX > (buffdhdx+uici*stride,stateDim,affineDim);
		dhdxiFD.setZero();
		x=x0;
		for(uicj=0;uicj<9;uicj++){
			x(uici)=x0(uici)+((double)uicj-4.0)*dx;
			dynSysLinearPart(t,x.data(),hi.data());
			dhdxiFD+=coeff[uicj]*hi;
		}
//		printf ( "x(%zu) error\n",uici );
		dhdxError=dhdxi-dhdxiFD/dx;
//		std::cout << dhdxError;
//		printf ( "\n-----------------\n" );

		assert(fabs(dhdxError.maxCoeff())<zeroTol);
	}
	
}

void affinesys::runningCostLinearPartFiniteDiffTest(double t,const double xin[]){
	Map<const VectorXd> x0(xin,stateDim);
	typedef Matrix<double,Dynamic,Dynamic,RowMajor> MatrixXX;
	VectorXd h0(affineDim);
	VectorXd dh0dxiError(affineDim);
	VectorXd dh0dxiFD(affineDim);
	double buffdhdx[500];
	Map<VectorXd> dh0dxi(NULL,0);
	VectorXd x(stateDim);

	static const double coeff[9]={	 1.0/280.0	,-4.0/105.0	, 1.0/5.0,
									-4.0/5.0	, 0.0		, 4.0/5.0, 	
									-1.0/5.0	, 4.0/105.0	,-1.0/280.0};
	static size_t uici,uicj;
	const double dx=0.01;
	runningCostLinearPart_1G(t,xin,buffdhdx);


	const size_t stride=affineDim;
	for(uici=0;uici<stateDim;uici++){
		new (&dh0dxi) Map< VectorXd > (buffdhdx+uici*stride,affineDim);
		dh0dxiFD.setZero();
		x=x0;
		for(uicj=0;uicj<9;uicj++){
			x(uici)=x0(uici)+((double)uicj-4.0)*dx;
			runningCostLinearPart(t,x.data(),h0.data());
			dh0dxiFD+=coeff[uicj]*h0;
		}
//		printf ( "x(%zu) error\n",uici );
		dh0dxiError=dh0dxi-dh0dxiFD/dx;
//		std::cout << dh0dxiError;
//		printf ( "\n-----------------\n" );
		assert(fabs(dh0dxiError.maxCoeff())<zeroTol);

	}
	
}
void affinesys::sfdttest(double t,const double xin[]){
	Map<const VectorXd> x0(xin,stateDim);
	typedef Matrix<double,Dynamic,Dynamic,RowMajor> MatrixXX;
	VectorXd dsdt(affineDim);
	VectorXd si(affineDim);
	VectorXd dsdtError(affineDim);
	VectorXd dsdtFD(affineDim);
	VectorXd x(stateDim);

	static const double coeff[9]={	 1.0/280.0	,-4.0/105.0	, 1.0/5.0,
									-4.0/5.0	, 0.0		, 4.0/5.0, 	
									-1.0/5.0	, 4.0/105.0	,-1.0/280.0};
	static size_t uici,uicj;
	const double dt=0.001;


	x=x0;
	double t0=t;
	swichingVector_partial_dt(t,xin,dsdt.data());
	dsdtFD.setZero();
	for(uicj=0;uicj<9;uicj++){
		t=t0+((double)uicj-4.0)*dt;
		computeSwichingVector(t,x.data(),si.data());
		dsdtFD+=coeff[uicj]*si;
	}
		printf ( "dsdt error\n");
		dsdtError=dsdt-dsdtFD/dt;
		std::cout << dsdtError;
		printf ( "\n-----------------\n" );
		assert(fabs(dsdtError.maxCoeff())<zeroTol);

	
}

//void affinesys::taskBeforeStateStepForTestEq_10(double t){
//		double r[10];
//		const double *x0;
//		const double *u0;
//		const double *cs;
//		double a,b,S;
//
//		x0= getState0(t);
//		u0= getControl0(t);
//		computeSwichingVector(t,iterator::state,&S);
//		swichingVector_dt_LinSys(t,iterator::state,&b,&a);
//		
//		r[0]=solver::R(t,state,control);
//		r[1]=solver::R(t,state,u0);
//		r[2]=solver::R(t,x0,u0);
//		r[3]=solver::phi(t,iterator::state);
//		r[4]=solver::phi(t,x0);
//		r[5]=(double)affinesys::convergetManifolds.size();
//		r[6]=solver::phi_dt(t,state,control);
//		r[7]=*u0;
//		r[8]=S;
//		r[9]=b+a**control;
//
////		runningCostLinearPartFiniteDiffTest(t,state);
//
////		Fl(t,iterator::state,NULL,iterator::control,&r[4]);
////		Cp(t,iterator::state,NULL,iterator::control,&r[5]);
//
//
////		sfdttest(t,state);
//		extraDataCont->push(t,r);
//
//}
//void affinesys::taskBeforeStateStepForTestR_6(double t){
//		double r[6];
//		const double *x0;
//		const double *u0;
//		const double *cs;
//		double a,b;
//
//		x0= getState0(t);
//		u0= getControl0(t);
//		
//		r[0]=R(t,state,control);
//		r[1]=R(t,state,u0);
//		r[2]=R(t,x0,u0);
//		r[3]=solver::G(t,iterator::state,iterator::control);
//		r[4]=solver::G(t,x0,u0);
//		r[5]=(double)affinesys::convergetManifolds.size();
//
//		extraDataCont->push(t,r);
//
//;/	}
