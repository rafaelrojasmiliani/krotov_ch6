
#include"solver.h"

#include<gsl/gsl_errno.h>
#include<iostream>
#include<eigen3/Eigen/Core>
#include <iomanip>
using namespace Eigen;


double solver::init(){
	static double trc,tfc;
	setState(x0);
	iterator::firstIteration();
	trc= *(state+stateDim);
	tfc=finalCost(0.0,state,control);
	return *(state+stateDim)+finalCost(0.0,state,control);
}
void solver::init(double *rc,double *fc){
	setState(x0);
	iterator::firstIteration();
	*rc= *(state+stateDim);
	*fc=finalCost(getEndTime(),state,control);
}
double solver::step(){
	Map<VectorXd> psiT(iterator::costate,stateDim);
	finalCost_1Gx(	getEndTime(),iterator::state,iterator::control,
					iterator::costate);
	
	psiT=-psiT;

	finalCost_2Gxx(getEndTime(),iterator::state,iterator::control,auxMem1);

	sym2full_cpy(jacMem,auxMem1,stateDim);
	jac=-jac;
	jac-=delta0.asDiagonal();



	full2sym_cpy(iterator::costate+stateDim,jacMem,stateDim);

//			sym2full_cpy(sigmaMem1,costate+stateDim,stateDim);
//			std::cout << "---------------sigma: \n" << sigma << "\n";
	solveCostate();

//	testCostateForward();
//	return 0.0;
	setState(x0);
	solveState();

//			printf ( "---------------------------\nextension data: \n" );
//			double G0=G(0.0,state0,control0,state0,control0,costate0);
//			double G =this->G(0.0,state,control,state0,control0,costate0);
//			printf ( "G  =  %lf \nG0 =  %lf\n\n\n",G,G0 );
	trc= *(state+stateDim);
	tfc=finalCost(0.0,state,control);
//			printf ( "rc= %lf    \t\t fc=  %lf \n", *(state+stateDim),finalCost(0.0,state,control));
	return *(state+stateDim)+finalCost(0.0,state,control);

}
//solver::solver(	size_t sdim, size_t cdim,
//		double cTimeStep,double endTime,const double *ic):
solver::solver(	size_t sdim, size_t cdim,
		double cTimeStep,double endTime,const double *ic,
			 double sae,double sre,
			 const gsl_odeiv2_step_type* sst,
			 double cae,double cre,
			 const gsl_odeiv2_step_type* cst
			 ):
	iterator(sdim,cdim,cTimeStep,endTime,sae,sre,sst,cae,cre,cst),
	jac(NULL,0,0),sigma(NULL,0,0),sigma_dot(NULL,0,0),psi(NULL,0),
	psi_dot(NULL,0),delta0(sdim),deltat(sdim),
	x0(ic),
	jacMem((double*)malloc(sdim*sdim*sdim*sizeof(double))),
	sigmaMem1((double*)malloc(sdim*sdim*sizeof(double))),
	sigmaMem2((double*)malloc(sdim*sdim*sizeof(double))),
auxMem1((double*)malloc(sdim*sdim*(sdim+(sdim*(sdim+1))/2)*sizeof(double)))
{
	
	
	
	new (&jac) 
		Map< Matrix<double,Dynamic,Dynamic,RowMajor> >(jacMem,
														stateDim,stateDim) ;
	new (&sigma) 
		Map< Matrix<double,Dynamic,Dynamic,RowMajor> >(sigmaMem1,
														iterator::stateDim,
														iterator::stateDim); 
	new (&sigma_dot) 
		Map< Matrix<double,Dynamic,Dynamic,RowMajor> >(sigmaMem2,
														iterator::stateDim,
														iterator::stateDim);

	delta0.setConstant(3.0);
	deltat.setConstant(3.0);

}
solver::~solver(){
	free (jacMem);
	free (sigmaMem2);
	free (sigmaMem1);
	free(auxMem1);
}
int solver::dynSysState(double time,const double x[],double xp[]){
		int res=dynSystem(time,x,iterator::control,xp);
//				*(xp+stateDim)=runningCost(time,x,iterator::control);
		return res;
}
int solver::dynSysCostate(double t,const double y[],double ydot[]){
	
	new (&psi) 
		Map< const VectorXd >(y,stateDim) ;
	new (&psi_dot) 
		Map< VectorXd >(ydot,stateDim) ;

	double *hessians;
	unsigned uich,uici,uick,uicj;
	double* const sdot_ptr=ydot+stateDim;
	const double* const s_ptr=y+stateDim;
	const double* const x0ptr=getState0(t);
	const double* const u0ptr=getControl0(t);

	#define sigmadot(i,j) (*(ydot+stateDim+trg2linsd(i,j)))
	#define sigma(i,j) (*(y+stateDim+trg2linsd(i,j)))


	runningCost_2Gxx(t,x0ptr,u0ptr,sdot_ptr);// ydot+stateDim=d2f^0dxdx
	for(uici=0;uici<sigmaDim;uici++)
		*(sdot_ptr+uici)= -*(sdot_ptr+uici);

	dynSys_2Gxx(t,x0ptr,u0ptr,jacMem);

	hessians=jacMem;

	for(uich=0;uich<stateDim;uich++){
		for(uici=0;uici<sigmaDim;uici++)
			*(sdot_ptr+uici)+= psi(uich)**(hessians+uici);
		hessians+= sigmaDim;
	}

	dynSys_1Gx(t,x0ptr,u0ptr,jacMem);


//	memset((void*)sdot_ptr,0,sigmaDim*sizeof(double));

	for(uici=0 ;uici < stateDim ;uici++){
		for(uicj=uici ;uicj < stateDim ;uicj++)
			for(uick=0 ;uick < stateDim ;uick++)
				*(sdot_ptr+trg2linsd(uici,uicj))+=
					 *(s_ptr+trg2linsd(uici,uick))*jac(uick,uicj)
					+*(s_ptr+trg2linsd(uicj,uick))*jac(uick,uici);

		*(sdot_ptr+trg2linsd(uici,uici))-= deltat(uici);

	}

//	printf ( "-----------------------------\n" );
//	for(uici=0 ;uici < stateDim ;uici++){
//		for(uicj=0 ;uicj < stateDim ;uicj++){
//				printf ( " %+14.7E ",*(sdot_ptr+trg2linsd(uici,uicj)) );
//		}
//		printf ( "\n" );
//	}
//	printf ( "-----------------------------\n" );
//		assert( 0 );
	runningCost_1Gx(t,x0ptr,u0ptr,ydot);

//	printf ( "psi dod 1 = %+14.7E\n",ydot[1] );                      
//	printf ( "psi dod 1 = %+14.7E\n",
//					2.0*u0ptr[0]*(x0ptr[1]-0.95*cos(0.95*t) ) +u0ptr[1]); 

	psi_dot= jac.transpose()*psi-psi_dot;


//	printf ( "psi dod 1 = %+14.7E\n",ydot[1] );                      
//	printf ( "psi dod 1 = %+14.7E\n",
//					-y[0]+2.0*(u0ptr[0]+0.1)*y[1]-  
//					2.0*(x0ptr[1]-0.95*cos(0.95*t) ) +u0ptr[1]); 

	#undef sigma
	#undef sigmadot

//		double p[4]={*(x0ptr+0),
//					 *(x0ptr+1),
//					 *(x0ptr+2),
//					 *(u0ptr+0)};
//	if(wwas==0){
//		auxCont.push(t,ydot);
//		auxCont2.push(t,p);
//	}else
//		auxCont3.push(t,p);
	return GSL_SUCCESS;
	
}

int solver::jacobianCostate
	(double t,const double y[],double dfdx[],double dfdt[]){

	static Map<Matrix<double,Dynamic,Dynamic,RowMajor> > 
					csJac(dfdx,costateDim,costateDim);
	static Matrix<double,Dynamic,Dynamic,RowMajor> sJac(stateDim,stateDim);
	const double* const x0=getState0(t);
	const double* const u0=getControl0(t);
	static VectorXd dfdxx(stateDim*stateDim*((stateDim*(1+stateDim))/2));
	size_t uici,uicj,uick,uicl,csJacRow,csJacCol;
	double auxMem[500];

	csJac.setZero();
	dynSys_1Gx(t,x0,u0,sJac.data());

	csJac.topLeftCorner(stateDim,stateDim)= sJac.transpose();

	dynSys_2Gxx(t,x0,u0,dfdxx.data());

	for(uici=0;uici<stateDim;uici++)
		csJac.col(uici).tail(sigmaDim)=
				dfdxx.segment(uici*sigmaDim,sigmaDim);

	for(uici= 0; uici < stateDim; uici++)
		for(uicj=uici; uicj < stateDim; uicj++)
			for(uick=0; uick < stateDim; uick++)
				for(uicl=uick; uicl < stateDim; uicl++){
					csJacRow= stateDim+trg2linsd(uici,uicj);
					csJacCol= stateDim+trg2linsd(uick,uicl);
					if(uicj == uick)
						csJac(csJacRow,csJacCol)+=sJac(uicl,uici);
					if(uici == uicj)
						csJac(csJacRow,csJacCol)+=sJac(uicl,uicj);
				}

//		static const double coeffFwr[7]={-49.0/20.0,	6.0,  	-15.0/2.0,    20.0/3.0,   
//										 -15.0/4.0,     6.0/5.0,-1.0/6.0};
//		const double dt=0.0001;
//
//		Map<VectorXd > csdV(auxMem,costateDim);
//		Map<VectorXd > csd_dt(dfdt,costateDim);
//		csd_dt.setConstant(0.0);
//
//		static double s0,si;
//		s0=s;
////
//		for(uicj=0;uicj<7;uicj++){
//			si=s0+((double)uicj)*dt;
//			dynSysCostate(si,y,auxMem);
//			csd_dt+=coeffFwr[uicj]*csdV;
//		}
//		csd_dt/=dt;
		return GSL_SUCCESS;
}
