/*!
 *       \file:  iteratordynsys.cpp
 *		 \brief: This file contains the definition of the class which represent
 *		 the process of iteration of the optimizers.
 *		 In particular here we have the classes which implment the iterations
 *		 as solution of ODES
 *        Created:  12/11/14 17:06:45
 *         Author:  Rafael A. Rojas 
 * =====================================================================================
 */

#include"iterator.h"
#include<assert.h>
#include<eigen3/Eigen/Core>
#include<iostream>
using namespace Eigen;

void iterator::firstIteration(){
	double t=0.0;
	*(state+stateDim)=0.0;
	processW->clear();
	int percentage=0,lp=-1;
	if(extraDataCont != NULL)
		extraDataCont->clear();
	taskInitStateSolution();
	do{
		controller0(t,state,control);
		processW->push(t,process);
		if(flagText){
			percentage= t/endTime*100.0;
			if(percentage % 5 == 0 && percentage!=lp){
				printf ( "solving iteration 0: %03i %% \r",percentage );
				fflush(stdout);
				lp=percentage;
			}
		}
		driverStateStepFirstIter(&t);
		taskAfterFirstIterStateStep(t);
	}while(t<=endTime);
	if(flagText)
			printf ( "solving iteration 0: 100 %% \n" );
	processW->push(t,process);
	lastEndTime=t;

}
void iterator::solveCostate(){
	wwas= 0;
	cSolCont *aux2;
	aux2=processW;
	processW=processR;
	processR=aux2;
	processR->gotoEndInTime();
	process1=processR->readBackwardInTime();
	double s=0.0;	
	double t;
	double ts;
	int percentage=0,lp=-1;
	costateCont.clear();
	do{	
		process0=processR->readBackwardInTime();
		t=*process1;
		ts=*process1 - *process0;
		if(ts < 0.0){
				printf ( "dt= %+.20E",ts );
				printf ( "t1 = %+.14E    \t t0 = %+.14E\n",
					  *process1,*process0 );
			assert(0);
		}
		if(flagText){
			percentage= s/endTime*100.0;
			if(percentage % 5 == 0 && percentage !=lp){ 
				printf ( "solving costate: %03i %% \r",percentage );
				fflush(stdout);
				lp=percentage;
			}
		}
		costateCont.push(*process1,costate);
		taskBeforeCostateStep(t);
		gsl_odeiv2_driver_apply(odeDriverCostate,&s,s+ts,costate);
		t=*process0;
		taskAfterCostateStep(t);
		process1=process0;
	}while(!processR->isAtBeginTime());
	taskBeforeCostateStep(t-2.0*ts);
	if(flagText)
	printf ( "solving costate: 100 %% \n" );
	costateCont.push(*process0,costate);
//		fprintf(costatep,"%+14.7E ",t);
//		for(size_t uici=0;uici<costateDim;uici++)
//			fprintf(costatep,"%+14.7E ",*(costate+uici));
//		fprintf(costatep,"\n");
//	costateCont.writeascii("costatexxx.dat");

}

void iterator::driverStateStep(double *ti){
		gsl_odeiv2_driver_apply(odeDriverState,ti,*ti+cTimeStep,state);
}
void iterator::driverStateStepFirstIter(double *ti){
		gsl_odeiv2_driver_apply(odeDriverState,ti,*ti+cTimeStep,state);
}
void iterator::driverCostateStep(double t1,double t0){
	assert(t1 > t0);
	static double s0,s1;
	s0=lastEndTime-t1;
	s1=lastEndTime-t0;
	gsl_odeiv2_driver_apply(odeDriverCostate,&s0,s1,costate);
}
int  iterator::costateDynSysWrapper
			(double t,const double y[],double yp[],void *params){
	t=static_cast<iterator*>(params)->lastEndTime-t;
	return static_cast<iterator*>(params)->dynSysCostate(t,y,yp);
}
void iterator::solveState(){


	wwas= 1;
	static double t;
	double dr0=0,dr1;
	static double taux;
	static unsigned int uici,uicj;
	static int counter=0;
	counter++;
	t=0.0;
	*(state+stateDim)=0.0;
	processR->gotoBeginInTime();
	costateCont.gotoBeginInTime();
	processW->clear();
	if(extraDataCont != NULL)
		extraDataCont->clear();
//	FILE *processp=fopen("processp.dat","w");
//	FILE *costatep=fopen("costatep.dat","w");
//	std::cout << Map<RowVectorXd>(state,stateDim) << "\n";
	const double* x0;
	const double* u0;
	const double* cs;
	const double* sigma0;
	int percentage=0,lp=-1;
	taskInitStateSolution();
	do{
		x0=	getState0(t);;
		u0=	getControl0(t);
		cs=	getCostate(t);
		sigma0= cs+stateDim;
		

// compute lm {{{
		for(uici=0;uici<stateDim;uici++){
			*(lm+uici)=*(cs+uici);
			for(uicj=0;uicj<stateDim;uicj++)
				*(lm+uici)+=*(sigma0+trg2lin(uici,uicj,stateDim))*
									(*(state+uicj)-*(x0+uicj));
		}// }}}
		
//		computeLm(t,state,lm);
		maxHamiltonian(t,state,lm,control);
		if(flagText){
			percentage= t/endTime*100.0;
			if(percentage % 5 == 0 && percentage !=lp){
				printf ( "solving state: %03i %% \r",percentage );
				fflush(stdout);
				lp=percentage;
			}
		}
		processW->push(t,process);
		taskBeforeStateStep(t);
		driverStateStep(&t);
		taskAfterStateStep(t);
		assert( t > taux);

	}while(t <= endTime);
		if(flagText)
			printf ( "solving state: 100 %% \n" );
	lastEndTime=t;
	taskBeforeStateStep(t);
	processW->push(t,process);
//	printf ( "\n" );
//	fclose(costatep);
//	fclose(processp);

//	printf ( "max dr = %+14.7E \n",dr0 );

}


int  iterator::stateDynSysWrapper(
				double t,const double y[],double yp[],void *params){
	*(yp+static_cast<iterator*>(params)->stateDim)=
		static_cast<iterator*>
		(params)->runningCost(t,y,static_cast<iterator*>(params)->control);

	return static_cast<iterator*>(params)->dynSysState(t,y,yp);
}
int  iterator::costateJacobianWrapper
		(double t,const double y[],double dfdx[],
		 double dfdt[],void *params){
	t=static_cast<iterator*>(params)->endTime-t;
	return static_cast<iterator*>(params)->
								jacobianCostate(t,y,dfdx,dfdt);
}

int  iterator::stateJacobianWrapper
		(double t,const double y[],double dfdx[],
		 double dfdt[],void *params){
	return static_cast<iterator*>(params)->
								jacobianState(t,y,dfdx,dfdt);
}
void iterator::taskBeforeStateStep(double t){

}
void iterator::taskAfterStateStep(double t){

}
