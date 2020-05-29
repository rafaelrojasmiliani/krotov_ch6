/*!
 *       \file:  iterator.h
 *
 *		 \brief: DESCRIPTION
 *        Version:  1.0
 *        Created:  12/11/14 16:38:24
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  ITERATOR_INC
#define  ITERATOR_INC
#include<gsl/gsl_odeiv2.h>
#include<eigen3/Eigen/Core>
#include<vector>
#include"csolcont.h"
/** This class is intended to be the abstraction of
 * the iteration process. It only make the io-task of
 * writing the distinct solutions in files an read of
 * the iterations and calls the GSL function to solve
 * different ODEs's required by the Krotov's method. */
class iterator{
	private:
		cSolCont processCont0;//!< Contains a process.
		cSolCont processCont1;//!< Contains a process.
		bool flagText;
		double endTime;

	protected:
	gsl_odeiv2_driver *odeDriverState,*odeDriverCostate;
	gsl_odeiv2_system odeStateSys,odeCostateSys;
	const gsl_odeiv2_step_type *odeStepState,*odeStepCostate;
	double lastEndTime;
	cSolCont *processR;//!< A pointer to the container where we are reading the last solution.
	cSolCont *processW;//!< A pointer to the container where we are storing the current solution.
	cSolCont costateCont;//!< The container which stores the costate solution.
	cSolCont *extraDataCont;//!< A pointer to a container which may store different data.
	const size_t stateDim;//!< Dimension of the dynamic system of the problem.
	const size_t controlDim;//!< Dimension of the control variable.
	const size_t costateDim;//!< Dimension of the costate variable.
	const size_t sigmaDim;//!< Auxiliary variable.
	double* const process;//!< Memory address where is stored the state, KPI, and control at one instant.
	double* const costate;//!<  Memory address where is stored the costate at one instant.
	double* const SigmaPtr;//!< Memory address where is stored the the matrix \f$\sigma\f$ of the costate.
	double* const state;//!< Memory address where is stored the state.
	double* const control;//!< Memory address there is stored the control.
	double* const lm;//!< Auxiliary memory to compute the approximate lagrange multipliers.
//	double csInterpol[100];//
//	double sInterpol[100];
	const double* process1;//!< A pointer used during the solution of the costate.
	const double* process0;//!< A pointer used during the solution of the costate.
	const double cTimeStep;//!< The time step for solving the ODEs. 
	int wwas;

//		iterator
//			(size_t s_dim,size_t c_dim,
//			 double cTimeStep,double EndTime);

		iterator
			(size_t s_dim,size_t c_dim,double cTimeStep,
			 double EndTime,
			 double sae=1.0e-7,double sre=1.0e-7,
			 const gsl_odeiv2_step_type* sst=gsl_odeiv2_step_rk8pd,
			 double cae=1.0e-7,double cre=1.0e-7,
			 const gsl_odeiv2_step_type* cst=gsl_odeiv2_step_rk8pd
			 );

/** A wrapper for the dynamic system used by GSL. */
		static int  stateDynSysWrapper(
						double t,const double y[],double yp[],void *params);
/** A wrapper for the jacobian used by GSL. */
		static int stateJacobianWrapper (
						double t, const double x[], 
						double *dfdy,double dfdt[], void *params);
/** A wrapper for the jacobian system used by GSL. */
		static int costateJacobianWrapper (
						double t, const double x[], 
						double *dfdy,double dfdt[], void *params);
/** A wrapper for the dynamic system used by GSL. */
		static int  costateDynSysWrapper(double t,const double y[],double yp[],void *params);

		/** Solves the controlled system and stores the solution in the a \ref cSolCont. */
		void solveState();
		virtual void taskInitStateSolution(){}
		/** A virtual function which allows to implement a task to
		 * perform after the integration step of the controlled system. */
		virtual void taskAfterStateStep(double t);
		/** A virtual function which allows to implement a task to
		 * perform after the integration step of the controlled system
		 * at the first iteration. */
		virtual void taskAfterFirstIterStateStep(double t){}
		/** A virtual function which allows to implement a task to
		 * perform before the integration step of the controlled system. */
		virtual void taskBeforeStateStep(double t);
		/** A virtual function which allows to implement a task to
		 * perform after the integration step of the costate. */
		virtual void taskAfterCostateStep(double t){}
		/** A virtual function which allows to implement a task to
		 * perform before the integration step of the costate. */
		virtual void taskBeforeCostateStep(double t){}
		/** Solve the costate equations and stores the solution in a \ref cSolCont. */
		void solveCostate();


		/** Perform the first integration to obtain the first guess
		 * of the optimization process. */
		virtual void  firstIteration();
		/** Get the value of the costate at time \ref t. */
		const double* getCostate(double t);
		/** Get the value of the last state at time \ref t. */
		const double* getState0(double t);
		/** Get the value of the last control at time \ref t. */
		const double* getControl0(double t);
		/** Compute an approximate derivative of the costate at time \ref t. */
		const double* getCostate_dt(double t);
	public:
		size_t getStateDim(){return stateDim;}
		size_t getControlDim(){return controlDim;}
		size_t getCostateDim(){return costateDim;}
		/** Return the final time \f$T\f$ of the Optimal Control Problem.  */
		double getEndTime(){return endTime;}
		/** Set The final time \f$T\f$ of the Optimal Control Problem. */
		void setEndTime(double t){endTime=t;}
		/** A function to get the \ref i, \ref j element of a square matrix
		 * of dimension \ref dim when this matrix is stored in a form
		 * which does not repeats elements for a symmetric matrix.  */
		size_t trg2lin(size_t i,size_t j,size_t dim) const;
		/** A function to get the \ref i, \ref j element of a square matrix
		 * of dimension \ref stateDim when this matrix is stored in a form
		 * which does not repeats elements for a symmetric matrix.  */
		size_t trg2linsd(size_t i,size_t j){                      
			if (i <= j)                                         
				return i*stateDim- ((i-1)*(i))/2 + j - i;   
			return j*stateDim - ((j-1)*(j))/2 + i - j;      
		}                                                       
//		double (controller1*)(double t,double x[]);

		/** Dynamic system of the costate. */
		virtual int dynSysCostate(double t,const double y[],double dydt[]) 	=0;
		/** Jacobian of the dynamic system of the costate. */
		virtual int jacobianCostate
				(double t,const double y[],
				 double dfdx[],double dfdt[]) 	=0;

		/** Jacobian of the controlled system.  */
		virtual int jacobianState
				(double t,const double y[],
				 double dfdx[],double dfdt[]) 	=0;

		/** Dynamic system of the state in a way which can be interfaced with GSL.. */
		virtual int dynSysState(double t,const double y[],double dydt[])		=0;

		/** Dynamic system of the state. */
		virtual int dynSystem(	double t,const double y[],const double u[],
								double yp[])=0;

		/** Returns the value of the control which maximizes the Hamiltonian of the system. */
		virtual void maxHamiltonian(double t,const double y[],const double *lambda,
									double *u)=0;
		/** Returns the running cost. */
		virtual double runningCost(	double t,const double y[],const double u[] )=0;

		virtual ~iterator();


		/** Set the value of the state variable. */
		void setState(const double*);
		/** Set the value of the costate variable. */
		void setCostate(const double*);
		/** Set the value of the control variable. */
		void setControl(const double*);


		/** Returns the control of the first iteration. */
		virtual void controller0(double t,const double x[],double u[])=0;

		/** Returns the equivalent integral defined by the Krotov's method. */
		virtual double R(double t,const double xin[],const double uin[])=0;
		/** Returns the equivalent final defined by the Krotov's method. */
		virtual double G(double t,const double xin[],const double uin[])=0;
		/** A debug function. */
		virtual void testR(double t,const double x[],const double u[],int level=1){}

		/** A debug function to numerically compute the second derivative of \f$R\f$. */
		virtual void d2Rdxdx(double t,const double xin[],const double uin[])=0;

		/** A debug function to numerically compute the first derivative of \f$R\f$. */
		virtual void dRdx(double t,const double xin[],const double uin[],double res[])=0;
		/** Updates the value of the variable \ref state for an integration step.*/	
		virtual void driverStateStep(double *ti);
		/** Updates the value of the variable \ref state for an integration step for the first iteration.*/	
		virtual void driverStateStepFirstIter(double *ti);
		/** Updates the value of the variable \ref costate for an integration step
		 * from time \ref ti to \ref tf.*/	
		virtual void driverCostateStep(double ti,double tf);
		/** Writes the stored solutions of the costate and the state to an ascii file
		 * with format \ref frm. */
		void writeFilesToText(const char* postfix="",const char* frm="%+14.7E ");
		/** Initialize the container for extra data. */
		void initExtraDataCont(size_t n);
		/** Returns a container which contains the last computed solution of the costate. */
		cSolCont* getLastCostate(){
			return costateCont.clone();
		}
		void setNoText(){
			flagText=false;
		}
		/** Compute the approximation of the lagrange multipliers. */
		virtual void computeLm(double t,const double xin[],double lm[]);
		/** Compute the approximation of the lagrange multipliers. */
		virtual void computeLm(double t,const double xin[],
							const double cs[],const double x0[],
							double lm[]);


};
#endif   /* ----- #ifndef ITERATOR_INC  ----- */
