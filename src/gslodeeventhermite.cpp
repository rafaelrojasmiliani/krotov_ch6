/*!
 *       \file:  gslodeeventhermite.cpp
 *
 *		 \brief: DESCRIPTION
 *        Version:  1.0
 *        Created:  19/06/15 09:58:16
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include"gslodeevent.h"
#include<eigen3/Eigen/Core>

void cGslOdeEvent::localInitHermiteInterpolator
					(double tI,const double xI[],double tF, const double xF[]){
	dynSys(tI,xI,localFI);
	dynSys(tF,xF,localFF);
	localTI= tI;
	localTF= tF;

}
void cGslOdeEvent::localHermiteInterpolator(double t,double xin[]){
	static Eigen::Map<const Eigen::VectorXd > f0(localFI,dim);
	static Eigen::Map<const Eigen::VectorXd > f1(localFF,dim);
	static Eigen::Map<const Eigen::VectorXd > x0(localXI,dim);
	static Eigen::Map<const Eigen::VectorXd > x1(localXF,dim);
	Eigen::Map< Eigen::VectorXd > x(xin,dim);
	static double theta;
	assert(localTI <= t && t <= localTF);
	theta= (t-localTI)/(localTF-localTI);

		x.noalias()=(1.0-theta)*x0 + theta*x1 +
				theta*(theta-1.0)*(
							(1.0-2.0*theta)*(x1-x0) +
							(theta-1.0)*timeStep*f0   +
							theta*timeStep*f1 );
}
