
#ifndef EH1DOFMB_H
#define EH1DOFMB_H

#include<affinesys.h>
#include<math.h>
class cEh1DofMb:public affinesys{
    public:
    double x0[2];
    double u0[1];
    double umin[1];
    double umax[1];
	double zeta0;
	double Omega;

    cEh1DofMb():
	affinesys(2,1,1,1.0e-3,20.0,1.0e-7,x0),
		zeta0 ( zeta0 ) ,
		Omega ( Omega ) 
	{
		x0[0]= 0.0;
		x0[1]= 0.0;
		u0[0]= 0.5;
		umax[0]= 10.0;
		umin[0]= 0.0;

	}

    virtual void controller0(double t,const double x[],double u[]){
		u[0] = 0.5;
    }

	int dynSystem(double t,const double x[], const double u[], double xd[]){
		xd[0]= x[1];
		xd[1]= // Not supported in C:
// etadd
-x[0] - x[1]*(2.0*u[0] + 2.0*zeta0) + etadd(t);
		return GSL_SUCCESS;
	}
	int dynSys_1Gx(double t,const double x[],const double u[],double dfdx[]){
		dfdx[0]= 0;
		dfdx[1]= 1.0;
		dfdx[2]= -1.0;
		dfdx[3]= -2.0*u[0] - 2.0*zeta0;
		return GSL_SUCCESS;
	}
    void dynSys_2Gxx(double t,const double x[],const double u[],
                                                             double dfdxx[]){
		 dfdxx[0]= 0;
		 dfdxx[1]= 0;
		 dfdxx[2]= 0;

		 dfdxx[3]= 0;
		 dfdxx[4]= 0;
		 dfdxx[5]= 0;

	}

    double runningCost(double t,const double x[],const double u[]){
        return -u[0]*pow(x[1], 2);
    }
    
    void runningCost_1Gx(double t,const double x[],const double u[],
                                double df0dx[]){
		df0dx[0]= 0;
		df0dx[1]= -2.0*u[0]*x[1];
	}

    void runningCost_2Gxx(double t,const double x[],const double u[],
                                double df0dxx[]){
		df0dxx[0]= 0;
		df0dxx[1]= 0;
		df0dxx[2]= -2.0*u[0];
	}

    double finalCost(double t,const double x[],const double u[]){
        return 0;
    }
    
    void finalCost_1Gx(double t,const double x[],const double u[],
                                                                    double *dFdx){
		dFdx[0]= 0;
		dFdx[1]= 0;
	}

    void finalCost_2Gxx(double t,const double x[],const double u[],
                               double *dFdxx){
		dFdxx[0]= 0;
		dFdxx[1]= 0;
		dFdxx[2]= 0;
	}

    double hamiltonian(double t,const double x[], const double lambda[],const double u[]){
        return // Not supported in C:
// etadd
lambda[0]*x[1] + lambda[1]*(-x[0] - x[1]*(2*u[0] + 2*zeta0) + etadd(t)) + u[0]*pow(x[1], 2);
    }

    void maxHamiltonian(double t,const double x[],
        const double lambda[],double u[]){

            compute_1Control(t,x,umin,umax,u,1.0e-6);
    }

    double etadd(double t){
        return sin(Omega*t);
            
	}

    void dynSysAffinePart
        (double t,const double x[],const double u2[],double g[]){
		g[0]= x[1];
		g[1]= // Not supported in C:
// etadd
-x[0] - 2*x[1]*zeta0 + etadd(t);
	}

        void dynSysLinearPart(double t,const double x[],double h[]){
		h[0]= 0;
		h[1]= -2.0*x[1];
	}

        void dynSysLinearPart_1G(double t,const double *x,double dhdx[]){
		dhdx[0]= 0;
		dhdx[1]= 0;
//--end of matrix

		dhdx[2]= 0;
		dhdx[3]= -2.0;
//--end of matrix

	}

        void dynSysLinearPart_dt(double t,const double *x,double dhdt[]){
		dhdt[0]= 0;
		dhdt[1]= 0;
	}

        double runningCostAffinePart(double t,const double *x,const double u2[]){
               return  0;
        }

    void runningCostLinearPart(double t,const double x[],double h0[]){
		h0[0]= -pow(x[1], 2);
	}

    void runningCostLinearPart_1G(double t,const double x[],double dh0dx[]){

		dh0dx[0]= 0;
		dh0dx[1]= -2.0*x[1];
	}

    void runningCostLinearPart_dt(double t,const double *x,double h0t[]){
		h0t[0]= 0;
	}

};
#endif
    