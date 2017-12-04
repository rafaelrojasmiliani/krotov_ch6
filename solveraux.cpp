
#include"solver.h"

#include<gsl/gsl_errno.h>
#include<iostream>
#include<eigen3/Eigen/Core>
#include <iomanip>
using namespace Eigen;
double solver::R(double t,const double xin[],const double uin[]){

	double result=0.0;
//			Map<const VectorXd> x0(getState0(t).data(),stateDim);
	
	const double* const csptr= getCostate(t);
	const double* const x0ptr= getState0(t);
	const double* const u0ptr= getControl0(t);
	Map<const VectorXd> x0(x0ptr,stateDim);
	Map<const VectorXd> psi(csptr,stateDim);
	sym2full_cpy(sigmaMem1,csptr+stateDim,stateDim);
	dynSysCostate(t,csptr,auxMem1);	
	sym2full_cpy(sigmaMem2,auxMem1+stateDim,stateDim);
	dynSystem(t,x0ptr,u0ptr,jacMem+stateDim);
	dynSystem(t,xin,uin,jacMem);

	static Map< VectorXd >
				psi_dot(auxMem1,stateDim) ;
	static Map< VectorXd > auxv(auxMem1,stateDim) ;

	Map <const VectorXd> x(xin,stateDim);
	Map <VectorXd> x_dot(jacMem,stateDim);
	Map <VectorXd> x_dot0(jacMem+stateDim,stateDim);

//			std::cout << "x \n" << x.transpose() << "\n-----\n";
	result-= psi_dot.transpose()*x;
//			std::cout << "psidot*x \n" << result << "\n-----\n";

//			std::cout << "sigma \n" << sigma << "\n-----\n";
	auxv=(sigma*(x_dot-x_dot0)-0.5*sigma_dot*(x-x0));

//			printf ( "\nsigma*df - 0.5 sigma*dx\n" );
//			std::cout << auxv.transpose() << "\n-----\n";

	result+=(x-x0).transpose()*auxv;

	result+= psi.transpose()*x_dot;

	result-=runningCost(t,xin,uin);


//	auxv= psi +sigma*(x_dot-x_dot0);
//			std::cout << "lm in R = " << auxv.transpose() << "\n";

//			printf ( "S in R= %+14.7E\n",S );
	return result;
	
}

void solver::d2Rdxdx(double t,const double xin[],const double uin[]){
	static size_t uici,uicj;
		const double dx=2.0;
		static double xi;
		VectorXd x(stateDim);
		double aux1,aux2,aux3,aux4;
	printf ( "-------------------------------------------\n" );
	for(uici=0;uici<stateDim;uici++){
		for(uicj=0;uicj<stateDim;uicj++){
		Map<const VectorXd> x0(xin,stateDim);
			x=x0;
			x(uici)+=  dx;
			x(uicj)+=  dx;
			aux1=R(t,x.data(),uin);

			x=x0;
			x(uici)+=  dx;
			x(uicj)+= -dx;
			aux2=R(t,x.data(),uin);

			x=x0;
			x(uici)+= -dx;
			x(uicj)+=  dx;
			aux3=R(t,x.data(),uin);

			x=x0;
			x(uici)+= -dx;
			x(uicj)+= -dx;
			aux4=R(t,x.data(),uin);
			printf ( " %+14.7E ",(aux1-aux2-aux3+aux4)/dx/dx/4.0 );
		}
		printf ( "\n" );
	}
//	const double *cs=getCostate(t);
//	const double *u0ptr=getControl0(t);
//	const double *x0ptr=getState0(t);
//		Map<const VectorXd> x0(x0ptr,stateDim);
//	double csd[9];
//	dynSysCostate(t,cs,csd);
//	double A= 10+*u0ptr;
//	double d= 6.5;
//	printf ( "dr2dxdx(1,3) = %+14.7E\n",  
//					-csd[5] + 2.0*A*x0(2)*cs[8]-4*(d-x0(0))*x0(2)**u0ptr);
//
//	printf ( "-----------------------+++++++++--------------------\n" );
}

void solver::testR(double t,const double x[],const double u[],int level){
		

	if(level==0) return;
	const double *x0,*u0,*cs;
	x0=getState0(t);
	u0=getControl0(t);
	cs=getCostate(t);

	int assertion=1;

	double rtxu,rtxu0,rtx0u0;

	rtxu=R(t,x,u);

	rtxu0=R(t,x,u0);

	rtx0u0=R(t,x0,u0);

	
	if(rtxu - rtxu0 <= -5.0e-5){
		printf ( "R(t,x,u) - R(t,x,u0) is not working\n" );
		printf ( "R(t,x ,u ) = %+14.7E\n",rtxu );
		printf ( "R(t,x ,u0) = %+14.7E\n",rtxu0 );
		printf ( "dR         = %+14.7E\n",rtxu-rtxu0 );
		assert(0);
	}
	if(level==1) return;


	if(rtxu0 - rtx0u0 <= -5.0e-5){
		printf ( "R(t,x ,u0) - R(t,x0,u0) is not working\n" );
		printf ( "R(t,x ,u ) = %+14.7E\n",rtxu0 );
		printf ( "R(t,x0,u0) = %+14.7E\n",rtx0u0 );
		printf ( "dR         = %+14.7E\n",rtxu0-rtx0u0 );
			VectorXd x1=Map<const VectorXd>(x0,stateDim);
			VectorXd x2=Map<const VectorXd>(x,stateDim);
			VectorXd x3=Map<const VectorXd>(x,stateDim);
			printf ( "variation in direction from x0 to x\n" );
				for(double tt=0.0;tt<=1.1;tt+=0.1){
					x3=tt*(x2-x1) + x1;
					std::cout << x3.transpose() << " ";
					printf ( "%+14.7E  %+14.7E\n",tt,R(t,x3.data(),u0) );
				}
		assert(0);
	}

	if(level==2) return;


}



void solver::dRdx(double t,const double xin[],const double uin[],double res[]){

		static const double coeff[9]={	 1.0/280.0	,-4.0/105.0	, 1.0/5.0,
							-4.0/5.0	, 0.0		, 4.0/5.0, 	
							-1.0/5.0	, 4.0/105.0	,-1.0/280.0};
		static size_t uici,uicj;
		const double dx=1.0;
		static double xi;
		Map<const VectorXd> x0(xin,stateDim);
		VectorXd x(stateDim);

		for(uici=0;uici< stateDim;uici++){
			x=x0;
			res[uici]=0.0;
			for(uicj=0;uicj<9;uicj++){
				x(uici)=x0(uici)+((double)uicj-4.0)*dx;
				res[uici]+=coeff[uicj]*R(t,x.data(),uin);
			}
			res[uici]/=dx;
		}
}

double solver::G(double t,const double xin[],const double uin[]){
	const double* const csptr= getCostate(t);
	const double* const x0ptr= getState0(t);
	const double* const u0ptr= getControl0(t);
	Map<const VectorXd> x0(x0ptr,stateDim);
	Map<const VectorXd> psi(csptr,stateDim);
	sym2full_cpy(sigmaMem1,csptr+stateDim,stateDim);
	Map <const VectorXd> x(xin,stateDim);
	static VectorXd dx(stateDim);
	static double result;
	dx=x-x0;
	result=  finalCost(t,xin,uin);
	result+= psi.transpose()*x;
	result+=0.5*dx.transpose()*sigma.selfadjointView<Upper>()*dx;

	return result;
}
double solver::phi(double t,const double xin[]){
	const double* const csptr= getCostate(t);
	const double* const x0ptr= getState0(t);
	const double* const u0ptr= getControl0(t);
	Map<const VectorXd> x0(x0ptr,stateDim);
	Map<const VectorXd> psi(csptr,stateDim);
	sym2full_cpy(sigmaMem1,csptr+stateDim,stateDim);
	Map <const VectorXd> x(xin,stateDim);
	static VectorXd dx(stateDim);
	static double result;
	dx=x-x0;
	result= 0.0;;
	result+= psi.transpose()*x;
	result+=0.5*dx.transpose()*sigma.selfadjointView<Upper>()*dx;

	return result;
}

double solver::phi_dt(double t,const double xin[],const double uin[]){

	double result=0.0;
//			Map<const VectorXd> x0(getState0(t).data(),stateDim);
	
	const double* const csptr= getCostate(t);
	const double* const x0ptr= getState0(t);
	const double* const u0ptr= getControl0(t);
	Map<const VectorXd> x0(x0ptr,stateDim);
	Map<const VectorXd> psi(csptr,stateDim);
	sym2full_cpy(sigmaMem1,csptr+stateDim,stateDim);
	dynSysCostate(t,csptr,auxMem1);	
	sym2full_cpy(sigmaMem2,auxMem1+stateDim,stateDim);
	dynSystem(t,x0ptr,u0ptr,jacMem+stateDim);
	dynSystem(t,xin,uin,jacMem);

	static Map< VectorXd >
				psi_dot(auxMem1,stateDim) ;
	static Map< VectorXd > auxv(auxMem1,stateDim) ;

	Map <const VectorXd> x(xin,stateDim);
	Map <VectorXd> x_dot(jacMem,stateDim);
	Map <VectorXd> x_dot0(jacMem+stateDim,stateDim);

	result-= psi_dot.transpose()*x;

	auxv=(sigma*(x_dot-x_dot0)-0.5*sigma_dot*(x-x0));

	result+=(x-x0).transpose()*auxv;

	result+= psi.transpose()*x_dot;

	return result;
	
}
double solver::findMin
		(double tol,bool flagtext,bool flagfile,const char* prefix,
		 size_t nd,size_t iar[],double d0[],double dt[]){
	static double val0,val1,val2;
	size_t i,uicj;
	char str[20];
	findMinIter=0;
	val2=val1=val0=init();
		if(flagtext)
		printf ( "[%03i] %+11.4E localimp = %+11.4E totalimp = %+11.4E \n",
				0,val2,(val2-val1)/fabs(val2), (val2-val0)/fabs(val0));
		if(flagfile){
			sprintf(str,"%s%03i",prefix,0);
			writeFilesToText(str);
		}
	uicj=0;
	for(i=1; i<500;i++){
		findMinIter++;
		if(nd!=0 && i==iar[uicj]){
			if(d0!=NULL)
				setDelta0(d0[uicj]);
			if(dt!=NULL)
				setDeltat(dt[uicj]);
			uicj++;
		}
		val2=step();
		if(flagtext)
		printf ( "[%03i] %+11.4E localimp = %+11.4E totalimp = %+11.4E \n",
				i,val2,(val2-val1)/fabs(val1), (val2-val0)/fabs(val0));
		if(flagfile){
			sprintf(str,"%s%03i",prefix,i);
			writeFilesToText(str);
		}

//		if(val2 > val1 && nd==0){
//			delta+=2.0;
//			mytest.setDelta(delta);
//		}


		if( fabs((val2-val1)/val2) < tol){
			if(flagtext)
			printf ( "optimal found! at iteration %zu\n",i );
			break;
		}
		val1=val2;
	}
	return val2;
}
void solver::taskBeforeStateStep(double t){


}


int  solver::costateDynSysWrapperF
			(double t,const double y[],double yp[],void *params){
	static int res;
	static size_t uici;
   res=static_cast<solver*>(params)->dynSysCostate(t,y,yp);
	for(uici=0;uici<static_cast<solver*>(params)->costateDim;uici++)
		yp[uici]=-yp[uici];
	return res;
}
void solver::testCostateForward(){
	
	const double *cs0;
	double *sol= (double*)malloc(costateDim*sizeof(double));
	double t;
	gsl_odeiv2_driver *drivercs;
	gsl_odeiv2_system syscs;
	syscs.function	= &this->costateDynSysWrapperF;
	syscs.jacobian	= NULL;
	syscs.dimension	= costateDim;
	syscs.params		= (void*)this;
	cSolCont solcont(costateDim);


	drivercs			= gsl_odeiv2_driver_alloc_y_new(	
											&syscs,
											this->odeStepCostate,
											cTimeStep*1.0e-8,
											1.0e-9,
											1.0e-9);
	gsl_odeiv2_driver_set_hmax(drivercs,cTimeStep);
	processR->gotoBeginInTime();
	costateCont.gotoBeginInTime();

	t=0;
	wwas=1;
	size_t uici;
	cs0=costateCont.getAtTimeInterp1(t);
	dynSysCostate(t,cs0,sol);
	solcont.push(t,sol);

	do{
		t+=cTimeStep;
		cs0=costateCont.getAtTimeInterp1(t);
		dynSysCostate(t,cs0,sol);
		solcont.push(t,sol);
	}while(t<=getEndTime());
	solcont.writeascii("costate_dt.dat");
//	auxCont3.writeascii("process_1.dat");
	
	free (sol);
}
