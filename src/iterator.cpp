/*!
 *       \file:  iterator.cpp
 *		 \brief: This file contains the definition of the class which represent
 *		 the process of iteration of the optimizers
 *        Created:  12/11/14 17:06:45
 *         Author:  Rafael A. Rojas 
 * =====================================================================================
 */
#include"iterator.h"
#include<gsl/gsl_odeiv2.h>
#include<iostream>

using namespace Eigen;
#define lpB 20
//iterator::iterator(size_t s_dim,size_t c_dim,double ts,double et):
iterator::iterator
			(size_t s_dim,size_t c_dim,double ts,
			 double et,
			 double sae,double sre,
			 const gsl_odeiv2_step_type* sst,
			 double cae,double cre,
			 const gsl_odeiv2_step_type* cst
			 ):
	processCont0(s_dim+c_dim+1),
	processCont1(s_dim+c_dim+1),
	flagText(true),
	processR(NULL),
	processW(NULL),
	costateCont(s_dim+((s_dim+1)*s_dim)/2,'-'),
	extraDataCont(NULL),
	stateDim(	s_dim	),
	controlDim(	c_dim	),
	costateDim(	s_dim+((s_dim+1)*s_dim)/2	),
	sigmaDim(	((s_dim+1)*s_dim)/2			),
	process(	(double*)malloc((s_dim+1+c_dim)*sizeof(double))),
	costate(	(double*)malloc((s_dim+((s_dim+1)*s_dim)/2)*sizeof(double))	),
	SigmaPtr(costate+s_dim),
	state(	process),
	control(process+s_dim+1),
	lm(			(double*)malloc(stateDim*sizeof(double))),
	endTime(et),
	cTimeStep(ts)
{

	this->odeStepState				= sst;
	this->odeStepCostate			= cst;

	this->odeStateSys.function		= &this->stateDynSysWrapper;
	this->odeStateSys.jacobian		= &this->stateJacobianWrapper;
	this->odeStateSys.dimension		= stateDim+1;
	this->odeStateSys.params		= (void*)this;

	this->odeCostateSys.function	= &this->costateDynSysWrapper;
	this->odeCostateSys.jacobian	= &this->costateJacobianWrapper;
	this->odeCostateSys.dimension	= costateDim;
	this->odeCostateSys.params		= (void*)this;

	this->odeDriverState			= 
								gsl_odeiv2_driver_alloc_y_new(	
											&this->odeStateSys,
											this->odeStepState,
											5.0e-10,
											sae,
											sre);

	printf ( "step type? %s\n",gsl_odeiv2_step_name(odeDriverState->s));
	gsl_odeiv2_driver_set_hmax(this->odeDriverState,cTimeStep);

	this->odeDriverCostate			= 
								gsl_odeiv2_driver_alloc_yp_new(	
											&this->odeCostateSys,
											this->odeStepCostate,
											cTimeStep*1.0e-8,
											cae,
											cre);
	gsl_odeiv2_driver_set_hmax(this->odeDriverCostate,cTimeStep);
	processW=&processCont0;
	processR=&processCont1;

}
iterator::~iterator(){

	free (lm);
	free (process);
	free (costate);
	gsl_odeiv2_driver_free(odeDriverState);
	gsl_odeiv2_driver_free(odeDriverCostate);
	if(extraDataCont != NULL)
		delete extraDataCont;

}

void iterator::setState(const double* ptr){	
	if(ptr!=NULL)
		memcpy((void*)state,(const void*)ptr,stateDim*sizeof(double));
	else
		memset((void*)state,0,sizeof(double)*stateDim);
}

void iterator::setCostate(const double* ptr){	
	if(ptr!=NULL)
		memcpy((void*)costate,(const void*)ptr,costateDim*sizeof(double));
	else
		memset((void*)costate,0,sizeof(double)*stateDim);
}

void iterator::setControl(const double* ptr){	
	if(ptr!=NULL)
		memcpy((void*)control,(const void*)ptr,controlDim*sizeof(double));
	else
		memset((void*)control,0,sizeof(double)*stateDim);
}
size_t iterator::trg2lin(size_t i,size_t j,size_t dim)const {
	if (i <= j)
		return i*dim - ((i-1)*((i-1) + 1))/2 + j - i;
	return j*dim - ((j-1)*((j-1) + 1))/2 + i - j;
}
const double* iterator::getCostate(double t){
	if(wwas == 0) return costate;
	return costateCont.getAtTimeInterp1(t);
//	static VectorXd f0(costateDim);
//	static VectorXd f1(costateDim);
//	static double theta;
//	static const double *x0ptr,*x1ptr;
//	static const double *u0ptr,*u1ptr;
//	const double *b0,*b1;
//	static double t0,t1;
//	static double ts;
//	static Map<VectorXd> xres(csInterpol,costateDim);
//	costateCont.getAtTimeDataToIterp(t,&b0,&b1);
//	t0=*b0;			x0ptr=b0+1;
//	t1=*b1;			x1ptr=b1+1;
//	Map<const VectorXd> x0(x0ptr,costateDim);
//	Map<const VectorXd> x1(x1ptr,costateDim);
//	theta=(t-t0)/(t1-t0);
//if(t0<= t && t <=t1){
//	dynSysCostate(t0,x0ptr,f0.data());
//	dynSysCostate(t1,x1ptr,f1.data());
//	ts=t1-t0;
//	assert(ts>0);
//
//	xres=(1.0-theta)*x0 +theta*x1 + theta*(theta-1.0)*
//			((1.0-2.0*theta)*(x1-x0)-(theta-1.0)*ts*f0-theta*ts*f1);
//
//}else{
//	xres=	x0+theta*(x1-x0);
//}
//	return csInterpol;
}
const double* iterator::getState0(double t){
	static VectorXd f0(stateDim);
	static VectorXd f1(stateDim);
	static double theta;
	static const double *x0ptr,*x1ptr;
	static const double *u0ptr,*u1ptr;
	const double *b0,*b1;
	static double t0,t1;
	static double ts;
//	static Map<VectorXd> xres(sInterpol,stateDim);
	if(wwas == 0){
//		assert((*process0-1.0e-8<= t && t <=*process1+1.0e-8));
		Map<VectorXd>(lm,stateDim) = 
				(Map<const VectorXd>(process1+1,stateDim) -
				 Map<const VectorXd>(process0+1,stateDim))
				/(*process1-*process0) *(t - *process0) 
				+Map<const VectorXd>(process0+1,stateDim);
		return lm;
	}
	
	return processR->getAtTimeInterp1(t);
//	processR->getAtTimeDataToIterp(t,&b0,&b1);
//	t0=*b0;			x0ptr=b0+1; u0ptr=b0+1+stateDim+1;
//	t1=*b1;			x1ptr=b1+1; u1ptr=b0+1+stateDim+1;
//	Map<const VectorXd> x0(x0ptr,stateDim);
//	Map<const VectorXd> x1(x1ptr,stateDim);
//	theta=(t-t0)/(t1-t0);
//if(t0<= t && t <=t1){
//	dynSystem(t0,x0ptr,u0ptr,f0.data());
//	dynSystem(t1,x1ptr,u1ptr,f1.data());
//	ts=t1-t0;
//
//	xres=(1.0-theta)*x0 +theta*x1 + theta*(theta-1.0)*
//			((1.0-2.0*theta)*(x1-x0)+(theta-1.0)*ts*f0+theta*ts*f1);
//
//}else{
//	xres=	x0+theta*(x1-x0);
//}
//	return sInterpol;
}
const double* iterator::getControl0(double t){
	const double *b0,*b1;
	if(wwas == 0){
		return process0+stateDim+2;
	}else
		processR->getAtTimeDataToIterp(t,&b0,&b1);
	return b0+stateDim+2;
}
const double* iterator::getCostate_dt(double t){
	const double *b0,*b1;
	if(wwas == 0) assert(0);
	costateCont.getAtTimeDataToIterp(t,&b0,&b1);

	Map<VectorXd>(lm,stateDim) = 
			(Map<const VectorXd>(b1+1,stateDim) -
			 Map<const VectorXd>(b0+1,stateDim))
			/(*b1-*b0);

	return lm;
}
void iterator::writeFilesToText(const char* postfix,const char* frm){
	char str[100];
	sprintf(str,"process%s.dat",postfix);
	processW->writeascii(str,frm);
	sprintf(str,"costate%s.dat",postfix);
	costateCont.writeascii(str,frm);
	if(extraDataCont !=NULL){
		sprintf(str,"extra%s.dat",postfix);
		extraDataCont->writeascii(str,frm);
	}

}

void iterator::initExtraDataCont(size_t n){
	extraDataCont= new cSolCont(n);
}

void iterator::computeLm(double t,const double xin[],double lm[]){
	size_t uici,uicj;
	const double* x0;
	const double* cs;
	const double* sigma0;
	x0=	getState0(t);;
	cs=	getCostate(t);
	sigma0= cs+stateDim;
	
// compute lm {{{
		for(uici=0;uici<stateDim;uici++){
			*(lm+uici)=*(cs+uici);
			for(uicj=0;uicj<stateDim;uicj++)
				*(lm+uici)+=*(sigma0+trg2lin(uici,uicj,stateDim))*
									(*(xin+uicj)-*(x0+uicj));
		}// }}}
}
void iterator::computeLm(double t,const double xin[],
							const double cs[],const double x0[],
							double lm[]){
	size_t uici,uicj;
	const double* sigma0;
	sigma0= cs+stateDim;
	
// compute lm {{{
		for(uici=0;uici<stateDim;uici++){
			*(lm+uici)=*(cs+uici);
			for(uicj=0;uicj<stateDim;uicj++)
				*(lm+uici)+=*(sigma0+trg2lin(uici,uicj,stateDim))*
									(*(xin+uicj)-*(x0+uicj));
		}// }}}
}
