#ifndef  SOLVER_INC
#define  SOLVER_INC
#include"iterator.h"
#include<gsl/gsl_errno.h>
#include<iostream>
#include <iomanip>
#include<eigen3/Eigen/Core>
using namespace Eigen;
/** This class contains the dynamics of the costate, 
 * assigns the initial values of the Cauchy problems
 * and implements the computation of different values 
 * to test the behaviour of the Krotov's method.*/
class solver: public iterator
{
	private:
		unsigned int findMinIter;
	protected:
		Map< Matrix<double,Dynamic,Dynamic,RowMajor> > jac;      //!< Auxiliary matrix 
		Map< Matrix<double,Dynamic,Dynamic,RowMajor> > sigma;    //!< Auxiliary matrix
		Map< Matrix<double,Dynamic,Dynamic,RowMajor> > sigma_dot;//!< Auxiliary matrix
		Map< const VectorXd > psi;//!< Auxiliary vector
		Map< VectorXd > psi_dot;  //!< Auxiliary vector
		VectorXd delta0;//!< Contains the diagonal values of the matrix \f$\Delta_0\f$
		VectorXd deltat;//!< Contains the diagonal values of the matrix \f$\Delta_t\f$
		
		double* const jacMem; //!< Address to auxiliary memory.
		double* const sigmaMem1;//!< Address to auxiliary memory.
		double* const sigmaMem2;//!< Address to auxiliary memory.
		double* const auxMem1;  //!< Address to auxiliary memory.
		const double* const  x0;//!< Address to the initial condition of the Optimal Control Problem.

		void sym2full_cpy(double* des,const double *src,const size_t dim){
				unsigned int uici,uicj;
				for(uici=0;uici<dim;uici++)
					for(uicj=0;uicj<dim;uicj++)
						*(des+dim*uici+uicj)=*(src+trg2lin(uici,uicj,dim));
		}
		void full2sym_cpy(double* des,const double *src,size_t dim){
			unsigned int uici;
			unsigned int offset=0;
			for(uici=0;uici<dim;uici++){
				memcpy((void*)(des+offset),(const void*)(src+dim*uici+uici),(dim-uici)*sizeof(double));
				offset+=dim-uici;
			}
		}
	public:
		/** Se the value of the diagonal elements of the matrix \f$\Delta_t\f$. */
		void setDeltat(double d){
			deltat.setConstant(d);
		}
		/** Se the value of the diagonal elements of the matrix \f$\Delta_0\f$. */
		void setDelta0(double d){
			delta0.setConstant(d);
		}
		/** Defines the dynamics of the controlled system if the Optimal Control Problem. */
		virtual int dynSystem(	double t,const double y[],const double u[],
								double yp[])=0;
		/** Defines the first gradient of the dynamics of the controlled system of the Optimal Control Problem. */
		virtual int dynSys_1Gx(	double t,const double y[],const double u[],
								double *dfdy)=0;
		/** Defines the second gradient of the dynamics of the controlled system of the Optimal Control Problem. */
		virtual void dynSys_2Gxx(	double t,const double y[],const double u[],
								double *dfidydy)=0;

		/** Defines the first gradient of the running cost of the Optimal Control Problem. */
		virtual void runningCost_1Gx(double,const double y[],const double u[],
									double *dfdydy)=0;
		/** Defines the second gradient of the running cost of the Optimal Control Problem. */
		virtual void runningCost_2Gxx(double t,const double y[],const double u[],
									double *dfdydy)=0;
		/** Returns the final cost of the Optimal Control Problem. */
		virtual double finalCost(	double t,const double y[],const double u[])=0;
		/** Returns the first gradient of the final cost of the Optimal Control Problem. */
		virtual void finalCost_1Gx(	double t,const double y[],const double u[],
							double *dfdy)=0;
		/** Returns the second gradient of the final cost of the Optimal Control Problem. */
		virtual void finalCost_2Gxx(double t,const double y[],const double u[],
							double *dfdydy)=0;

		double tfc;
		double trc;
		/** Initialized the optimization process. */
		double 	init();
		/** Initialized the optimization process. This function returns the 
		 * value of the running cost, \ref rc, and the final cost \ref fc. */
		void 	init(double *rc,double *fc);
		/** Performs an step of the optimization process. */
		double step();
//		solver(	size_t sdim, size_t cdim,
//				double cTimeStep,double endTime,const double *ic);
		solver(	size_t sdim, size_t cdim,
				double cTimeStep,double endTime,const double *ic,
			 double sae=1.0e-7,double sre=1.0e-7,
			 const gsl_odeiv2_step_type* sst=gsl_odeiv2_step_rk8pd,
			 double cae=1.0e-7,double cre=1.0e-7,
			 const gsl_odeiv2_step_type* cst=gsl_odeiv2_step_rk8pd
			 );
		virtual ~solver();
		/** Returns the dynamic equation of the state in a way to be implemented
		 * in the GSL functions. */
		int dynSysState(double time,const double x[],double xp[]);
		/** Returns the dynamic equations of the costate. */
		int dynSysCostate(double t,const double y[],double ydot[]);
		/** Returns the equivalent integrant defined by the Krotov's method. */
		virtual double R(double t,const double xin[],const double uin[]);
		/** A debug function. */
		virtual void testR(double t,const double x[],const double u[],int level=1);
		/** A debug function to numerically compute the second derivative of \f$R\f$. */
		void d2Rdxdx(double t,const double xin[],const double uin[]);
		/** A debug function to numerically compute the first derivative of \f$R\f$. */
		void dRdx(double t,const double xin[],const double uin[],double res[]);
		/** Returns the equivalent final defined by the Krotov's method. */
		double G(double t,const double xin[],const double uin[]);
		/** Returns the value of the bounding function. */
		double phi(double t,const double xin[]);
		/** Returns the value of the derivative bounding function. */
		double phi_dt(double t,const double xin[],const double uin[]);
		/** Performs the optimization procedure with a tolerance of \ref tol.. */
		double findMin(double tol=1.0e-4,bool flagtext=false,bool flagfile=false,
						const char* prefix="",
						size_t dn=0,size_t *iar=NULL,double *d0=NULL,double *dt=NULL);
		virtual void taskBeforeStateStep(double t);
		virtual int jacobianState(
					double t,const double y[],double dfdx[],double dfdt[]){}
		/** Computes the jacobian of the costate dynamic equations. */
		virtual int jacobianCostate(
					double t,const double y[],double dfdx[],double dfdt[]);
		/** A debug function. */
		void testCostateForward();
		/** A debug function. */
		static int  costateDynSysWrapperF
			(double t,const double y[],double yp[],void *params);
		/** A debug function. */
		unsigned int getFindMinIter(){
			return findMinIter;
		}
};
#endif   /* ----- #ifndef SOLVER_INC  ----- */
