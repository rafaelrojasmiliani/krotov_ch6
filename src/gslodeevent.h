/*!
 *       \file:  gslodeevent.h
 *
 *		 \brief: DESCRIPTION
 *        Version:  1.0
 *        Created:  19/06/15 09:56:00
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */


#ifndef  GSLODEEVENT_INC
#define  GSLODEEVENT_INC
#include<gsl/gsl_odeiv2.h>
#include<eigen3/Eigen/Core>
using namespace Eigen;
/** This class implements the integration methods available in the GSL library to
 * perform an event detection integration of ODE's.
 * An event detection integration integrates the ODE
 * \f{eqnarray}{
 *	\dot{x} &=& f(t,x) 
 * \f} 
 * between times \f$t_0\f$ and \f$t_1\f$, but it stops at time \f$t_i\f$ with \f$t_0< t_i< t_1\f$
 * whenever some component of a vector \f$S(t,x)\f$ meets a root. To obtain this value different
 * root finding method are implement.
 * The integration algorithm is the following:
 * -# Integrate the ODE from the point \f$x_0\f$ at time \f$t_0\f$ to 
 *  time \ft_1\f$ and returns the point \f$x_1\f$. 
 * -# Test if some component of the switching vector \f$S(t,x)\f$ have
 *  changed sign between \f$(t_0,x_0)\f$ and \f$(t_1,x_1)\f$. 
 * -# If some component of \f$S(t,x)\f$ have changed its sign, we perform
 *  a root finding method to obtain the root of this component. This root
 *  is found using an Hermite interpolating polynomial between \f$(t_0,x_0)\f$ and \f$(t_1,x_1)\f$.
 * */
class cGslOdeEvent{
	private:
		double* const localFI;//!< Auxiliary value of interpolation.
		double* const localFF;//!< Auxiliary value of interpolation.
		double* const localXI;//!< Auxiliary value of interpolation.
		double* const localXF;//!< Auxiliary value of interpolation.
		double* const xAuxRoot;//!< Auxiliary value of interpolation.
		double localTI;//!< Auxiliary value of interpolation.
		double localTF;//!< Auxiliary value of interpolation.
		gsl_odeiv2_driver *gslOdeDriver;
		gsl_odeiv2_system gslOdeSys;
		static int  gslOdeFuntionWrapper(
						double t,const double y[],double yp[],void *params){
			return static_cast<cGslOdeEvent*>(params)->dynSys(t,y,yp);

		}
		static int gslOdeJacobianWrapper(
						double t, const double x[], 
						double *dfdx,double dfdt[], void *params){

			return static_cast<cGslOdeEvent*>(params)->gslOdeEventJacobian(t,x,dfdx,dfdt);
		}
	
		/** Rider Method for root finding. */
		void localRootFindRidder(double *time,double x[]);
		/** Rider Method for root finding. */
		int localRootFindBrent(double *time,double x[]);
		/** Rider Method for root finding. */
		int localRootFindBrent(double *time,double x[],size_t i);
		/** Initializes Hermite interpolator. */
		void localInitHermiteInterpolator(
						double tI,const double xI[],double tF, const double xF[]);
		/** Get Hermite Interpolation. */
		void localHermiteInterpolator(double t,double x[]);
	public:
		const double zeroTol;//!< Tolerance.
		const double timeStep; //!< Time step of the integration method.
		const size_t dim;//!< Dimension of the dynamic system.
		const size_t edim; //!< Dimension of the event vector.
		const gsl_odeiv2_step_type* const gslOdeStepType;
		cGslOdeEvent
				(size_t d1,size_t d2, double ts,double zt,
				 double errAbs=1.0e-7,double errRel=1.0e-7,
				 const gsl_odeiv2_step_type* st=gsl_odeiv2_step_rk8pd,const double *x0=NULL):
		localFI(	(double*)malloc(d1*sizeof(double))),
		localFF(	(double*)malloc(d1*sizeof(double))),
		localXI(	(double*)malloc(d1*sizeof(double))),
		localXF(	(double*)malloc(d1*sizeof(double))),
		xAuxRoot(	(double*)malloc(d1*sizeof(double))),
		zeroTol(zt),
		timeStep(ts),
		dim(d1),
		edim(d2),
		gslOdeStepType(st)
		{

			double ts1=ts/10.0;
			gslOdeSys.function		= &gslOdeFuntionWrapper;
			gslOdeSys.jacobian		= &gslOdeJacobianWrapper;
			gslOdeSys.dimension		= dim;
			gslOdeSys.params		= (void*)this;

			if(x0!=NULL)
				ts1=computeFirstStep(x0,errAbs,errRel);



			gslOdeDriver			= gsl_odeiv2_driver_alloc_standard_new(	
											&gslOdeSys,
											gslOdeStepType,
											ts1,
											errAbs,
											errRel,1.0,10.0);

			gsl_odeiv2_driver_set_hmax(gslOdeDriver,ts);

		}
		virtual ~cGslOdeEvent(){
			gsl_odeiv2_driver_free(gslOdeDriver);
			free(localFI);
			free(localFF);
			free(localXI);
			free(localXF);
			free(xAuxRoot);
		}
		/** Computes first integration step as the book of Hairer says.*/
		double computeFirstStep(const double *x0in,double ea,double er){
			static Map<const VectorXd> x0(x0in,dim);
			VectorXd f0(dim);
			VectorXd f1(dim);
			VectorXd x1(dim);

			dynSys(0.0,x0.data(),f0.data());
			double d0,d1,d2,sc,h0,h1;

			size_t uici,p;

			d0=x0.norm();
			d1=f0.norm();
			if(d0 > 1.0e-5 && d1 > 1.0e-5 )	
				h0= 0.01*d0/d1;
			else
				h0= 1.0e-6;
			x1=x0 + h0*f0;
			dynSys(h0,x1.data(),f1.data());

			d2=(f1-f0).norm()/h0;

			p=getOrder();
			if(d0 > 1.0e-15 && d1 > 1.0e-15)
				h1= pow(0.01/((d1> d0)?d1:d0),-(p+1.0) );
			else	
				h1= (1.0e-6 > h0*1.0e-3)?1.0e-6:h0*1.0e-3;

			return h1;
		}

		/** Get the order of the integration method. */
		int getOrder(){
				if( gsl_odeiv2_step_rk8pd == gslOdeStepType ) 
					return 8;
				if( gsl_odeiv2_step_rk2 == gslOdeStepType ) 
					return 2;
				if( gsl_odeiv2_step_rk2imp == gslOdeStepType ) 
					return 2;
			return 4;
		}
		/** Dynamic system which we are integrating. */
		virtual int dynSys(double t,const double x[],double xp[])=0;
//		void hermiteInterpolator
//				(double t0,const double x0[],double t1,const double x1[],double t,
//				 double x[]);
		/** Event Vector.. */
		virtual int eventFucntion(double t,const double xin[],double h[])=0;

		/** The jacobian of the dynamic system may be required by the integration method. */
		virtual int gslOdeEventJacobian
				(double t,const double x[],double dfdx[],double dfdt[])
		{}

		virtual void projectonSF0(double t,double *x){}

		/** Perform the event detection integration method. */
		void stepEventDectector(double *t0,double x0[]);


};

//#include"gslodeeventhermitetemplate.h"
//#include"gslodeeventrootriddertemplate.h"
//#include"gslodeeventsteptemplate.h"

#endif   /* ----- #ifndef GSLODEEVENTTEMPLATE_INC  ----- */
