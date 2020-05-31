/*!
 *       \file:  test.cpp
 *
 *		 \brief: DESCRIPTION
 *        Version:  1.0
 *        Created:  30/06/15 14:09:56
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */
#include<iterator.h>
#include<solver.h>
#include<affinesys.h>
#include"ehcapacitor.h"

#include <fenv.h>
#include<string.h>
#include<stdio.h>
#include<stdlib.h>

#define CLASS cEhCapacitor

class myTest:public CLASS{
    public:
        double y(double t){
            return 0.015*sin(t);
        }
        void controller0(double t,const double x[],double u[]){
            double xd[3];
            double uCharging[2]={umax[0],umin[1]};
            double uDischarging[2]={umin[0],umax[1]};

            if(x[1]<= 0){
                u[0]=umax[0]; u[1]=umin[1];
            }else{
                u[0]=umin[0]; u[1]=umax[1];
            }
        }
};



int main(){
    double val0,val1,val2,aux;
    int i;

    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

    myTest mytest;

    mytest.setDelta0(0.0);
    mytest.setDeltat(0.02);

    mytest.findMin(1.0e-3,true,true);
    mytest.writeFilesToText("000");
}




