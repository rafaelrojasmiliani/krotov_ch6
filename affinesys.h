

#ifndef  AFFINESYS_INC
#define  AFFINESYS_INC
#include<iostream>
#include<algorithm>
#include<vector>
#include"solver.h"
#include<eigen3/Eigen/Core>
#include"indexvector.h"
#include"gslodeevent.h"
using namespace Eigen;
/** This class is an abstraction of an optimal control problem
 * which is affine with respect to the control variable.
 *
 * Represents the following Optimal Control problem
 *
 * \f{eqnarray}{
 *	\min I &=& F\big(T,x(T)\big) +\int_0^T g^0(t,x) + h^0(t,x)\cdot u \ dt \cr
 *	\dot{x} &=& g(t,x) + h(t,x) u 
 * \f} 
 *	where \f$x\in\mathbb{R}^n\f$, \f$u\in\mathbb{R}^m\f$, 
 *	\f$ h^0(t,x)\in\mathbb{R}^m  \f$, \f$g(t,x)\in\mathbb{R}^n\f$,
 *  \f$h(t,x)\in\mathbb{R}^{n\times}m\f$. 
 *  The solution of the Krotov's method for this problem requires special
 *  event-detection features of the solution the arising differential equations
 *  that are implemented in the class \ref cGslOdeEvent. 
 *  This class define several pure virtual functions in that represents
 *  different components of the optimal control problem that this class
 *  intent to solve. Also this class implements different features required
 *  to integrate the arising differential equation with discontinuous right-hand
 *  side that arises during the Krotov's optimization process. 
 *  - Computation of the swishing vector \f$S(t,x)\f$.
 *  - Computation of the equivalent control in the sliding surface. 
 *  - Projection on the sliding manifold.
 *
 * */
class affinesys: public solver,public cGslOdeEvent{
	private:

//		double* const swichingVector;
//		double* const swichingVector_dt;
		affinesys& operator=(affinesys &that);
		Map< VectorXd > f0;
		Map< VectorXd > g;
		Map< Matrix<double,Dynamic,Dynamic,RowMajor> > A;
		Map< Matrix<double,Dynamic,Dynamic,RowMajor> > dhdxi;
		Map< Matrix<double,Dynamic,Dynamic,RowMajor> > dhdt;
		Map< Matrix<double,Dynamic,Dynamic,RowMajor> > m1;
		Map< Matrix<double,Dynamic,Dynamic,RowMajor> > skxl;

		
		/// Implements the dynamic system of \ref cGslOdeEvent
		int dynSys(double t,const double xin[],double xdot[]){
			dynSystem(t,xin,iterator::control,xdot);
			xdot[stateDim]= runningCost(t,xin,iterator::control);
			return GSL_SUCCESS;
		}
		/// Implements the event Function of \ref cGslOdeEvent
		int eventFucntion(double t,const double xin[],double h[]){
			computeSwichingVector(t,xin,h);
			return GSL_SUCCESS;
		}

		double proyTolSfVal;
		double proyTolDx;
		unsigned long proyMaxIter;
	protected:
		/** Is a vector containing the values of the swishing vector
		 * which are zero and define the sliding manifold. */
		tIndexvector convergetManifolds;
		double DEBUG_b1;
		double DEBUG_b2a;
		double DEBUG_b2b;
		double DEBUG_b3a;
		double DEBUG_b3b;
		double DEBUG_a1;
		double DEBUG_a2;
		double DEBUG_a3;
		double localProyGain;
	public:
		const size_t affineDim;//!< The numbers of controls which are affine.

		affinesys(	size_t sdim,size_t cdim,size_t afdim,
					double ts,double et,double zt,const double *ic,
					double ae=1.0e-7,double re=1.0e-7,
					const gsl_odeiv2_step_type* st=gsl_odeiv2_step_rk8pd,
					double csae=1.0e-7,double csre=1.0e-7,
					const gsl_odeiv2_step_type* csst=gsl_odeiv2_step_rk8pd
					):
		solver(sdim,cdim,ts,et,ic,ae,re,st,csae,csre,csst),
		cGslOdeEvent(sdim+1,afdim,ts,zt,ae,re,st),
		f0(NULL,0),g(NULL,0),A(NULL,0,0),dhdxi(NULL,0,0),
		dhdt(NULL,0,0),m1(NULL,0,0),skxl(NULL,0,0),
		affineDim(afdim)
		{
			localProyGain=0.01;
			proyTolSfVal=zeroTol;
			proyTolDx=zeroTol;
			proyMaxIter=100;
		}

		/** Set a gain for the projection process. */
		void setProyLocalGain(double p){
			localProyGain=p;
		}
		virtual ~affinesys(){}
		/** Returns \f$g(t,x)\f$. */
		virtual void dynSysAffinePart(double t,const double *x,const double u2[], double *g)=0;
		/**Returns \f$h(t,x)\f$.
		 * This method return a matrix stateDim*affineDim, 
		 * that is the linear 
		 * part of the affine system. */
		virtual void dynSysLinearPart(double t,const double *x,double *h)=0;
		/**Returns the first gradient of  \f$h(t,x)\f$.
		 * This method returns a third order tensor in a linear array.
		 * this linear array is composed by blocks of stateDim*controDim,
		 * matrices given in row major form. each matrix represents the
		 * linear part derivative with respect to each state component.*/
		virtual void dynSysLinearPart_1G(double t,const double *x,double *dhdx)=0;
		/** Returns the time derivative of \f$h(t,x)\f$.
		 * It returns a linear array with the components in row major format. */
		virtual void dynSysLinearPart_dt(double t,const double *x,double *dhdx)=0;
		/** Returns \f$g^0(t,x)\f$. */
		virtual double runningCostAffinePart(double t,const double *x,const double *u2)=0;
		/** Returns \f$h^0(t,x)\f$.
		 * Returns a linear array of \ref affineDim elements. */
		virtual void runningCostLinearPart(double t,const double *x,double *r)=0;
		/** Returns the first gradient of \f$h^0(t,x)\f$.
		 * Returns a linear array of \f$n\times m\f$ elements. Each block
		 * of \f$m\f$ elements if the partial derivative of \f$h^0(t,x)\f$ 
		 * with respect to a state component \f$x\f$.*/
		virtual void runningCostLinearPart_1G(double t,const double *x,double *r)=0;
		/** Returns the time derivative of \f$h^0(t,x)\f$.
		 * Returns a linear array of \f$m\f$ elements.*/
		virtual void runningCostLinearPart_dt(double t,const double *x,double *r)=0;


		/** Computes the swishing vector.
		 *  Computes the swishing vector using the following formula
 		 * \f{eqnarray}{
			S^k &=& \bigg( \psi_n + \sigma_{nj} \Delta x_j \bigg)	h_{nk} - h^0_k
		 * \f}	
		 * In this function the costate values are computed.
		 * */
		virtual void computeSwichingVector(double t,
									const double xin[],double *res,int flag=0);
		/** Computes the swishing vector.
		 *  Computes the swishing vector using the following formula
 		 * \f{eqnarray}{
			S^k &=& \bigg( \psi_n + \sigma_{nj} \Delta x_j \bigg)	h_{nk} - h^0_k
		 * \f}	
		 * In this function the costate values must be furnished.
		 * */
		virtual void computeSwichingVector(double t,
									const double xin[],const double cs[],
									const double x0in[],double *res,int flag=0);

		/** Computes the derivative of the swishing vector.
		 * Computes the derivative of the swishing vector using the following 
		 * formula:
 		 * \f{eqnarray}{
		\dot{S}^k &=&  A^k_i(t,x) u^i + b^k(t,x)
		 * \f}	
		 * where
 		 * \f{eqnarray}{
				 b^k(t,x) &=&
		\left[ \dot{\psi}_n(t) + \dot{\sigma}_{nj}(t)  g^j(t,x) - \sigma_{nj}f_i^j  \right] h^n_k 
			+\left( \psi_n + \sigma_{nj} \Delta x^j \right) 
			\left[ \frac{\partial h^n_k}{\partial x^l} g^l + \frac{\partial h^n_k}{\partial t} \right] \cr
			&& + \sigma_{nj} g^j h^i_k - \frac{\partial h^0_k}{\partial x^l} g^l - \frac{\partial h^0_k}{\partial t}
		 * \f}	
		 * and
 		 * \f{eqnarray}{
				A^k_r &=&  \bigg( \sigma_{nj} h^j_r h^n_k + \big(\psi_n +\sigma_{nj} \Delta x^j\big)
			\frac{\partial h^n_k}{\partial x^l} h^l_r -  \frac{\partial h^0_k}{\partial x^l} h^l_r\bigg) 
		 * \f}	
		 *
		 * */
		virtual void swichingVector_dt_LinSys(double t,
									const double xin[],
									double *bout,double *Aout);
		/** Returns the partial derivative of the switching vector with respect to time. */
		virtual void swichingVector_partial_dt(double t,
									const double xin[],double *res);
		/** Returns the first gradient of the switching vector.
		 * This values is returned as a linear array \ref res of \ref affineDim \f$\times\f$ \ref stateDime elements
		 * which represents a row major stored matrix. */
		virtual void swichingVector_G1(double t,
									const double xin[],double *res);
		/** Returns the first gradient of the switching vector.
		 * This values is returned as a linear array \ref res of \ref affineDim \f$\times\f$ \ref stateDime elements
		 * which represents a row major stored matrix. It needs the value of the costate and the previous process.*/ 
		virtual void swichingVector_G1(double t,
									const double xin[],const double cs[],
									const double x0in[],double *res);
		
		/** Define the integration step. Integrates as \ref cGslOdeEvent and project on 
		 * the sliding manifold when it is necessary. */
		virtual void driverStateStep(double *ti){
			stepEventDectector(ti,iterator::state);
			projectonSF0(*ti,iterator::state);
		}


		/** Computes the projection of a point \f$x\in \mathbb{R}^n\f$ into the sliding manifold. */
		virtual void projectonSF0(double t,double *x);
		/** Computes the projection of a vector \f$v\in \mathbb{R}^n\f$ into the tangent space to the sliding manifold. */
		virtual void projectonSF0Tangent(double t,const double xin[],double vin[]);
		virtual void projectonSF0Time(double *t,double xin[]);
		/** Debug function. */
		virtual void testR(double t,const double x[],const double u[],int level=1);

		/** Implements \ref  gslOdeEventJacobian of \ref cGslOdeEvent.*/
		int jacobianState(
					double t,const double y[],double dfdx[],double dfdt[]){ 
				return gslOdeEventJacobian(t,y,dfdx,dfdt);
		}
		/** Sets the parameter \ref proyTolSfVal of the projection algorithm. */
		void setProyTolSfVal(double t){
			proyTolSfVal=t;
		}
		/** Sets the parameter \ref proyTolDx of the projection algorithm. */
		void setProyTolDx(double t){
			proyTolDx=t;
		}
		/** Sets the parameter \ref proyMaxIter  of the projection algorithm. */
		void setProyMaxIter(unsigned long Nin){
			proyMaxIter=Nin;
		}

		/** Computes the equivalent control on the sliding manifold. */
		void computeControl(double t,const double xin[],
				const double umin[],const double umax[],double uout[],double zt);

		/** Computes the equivalent control on the sliding manifold when only one control is needed. */
		void compute_1Control(double t,const double xin[],
				const double umin[],const double umax[],double uout[],
				double zt);

		/** Computes the equivalent control on the sliding manifold when only two  control is needed. */
		void compute_2Control(double t,const double xin[],
				const double umin[],const double umax[],double uout[],
				double zt);

		/** A debug function to test the derivatives of \f$ h(t,x)\f$. */
		void dynSysLinearPartFiniteDiffTest(double t,const double xin[]);
		/** A debug function to test the derivatives of \f$ h^0(t,x)\f$. */
		void runningCostLinearPartFiniteDiffTest(double t,const double xin[]);
		/** A debug function to test the derivatives of \f$ S(t,x)\f$. */
		void sfdttest(double t,const double xin[]);
//		void taskBeforeStateStepForTestEq_10(double t);
};
#endif   /* ----- #ifndef AFFINESYS_INC  ----- */
