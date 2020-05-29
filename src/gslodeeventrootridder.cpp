/*!
 *       \file:  gslodeeventrootridder.h
 *
 *		 \brief: DESCRIPTION
 *        Version:  1.0
 *        Created:  19/06/15 10:00:40
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include"gslodeevent.h"
#include<assert.h>
#include<math.h>

void cGslOdeEvent::localRootFindRidder(double *time,double x[]){
	
	double ef0,ef1,efm,ef2;
	double eta0,eta1,etam,eta2;
	double t0,t1,t2,tm;

	eta0=0.0;
	etam=0.5;
	eta1=1.0;
	int counter=0;
	double efh[200];

	do{
		t0= localTI + eta0*timeStep;
		t1= localTI + eta1*timeStep;
		tm= localTI + etam*timeStep;


		localHermiteInterpolator(t0,xAuxRoot);
		eventFucntion(t0,xAuxRoot,&ef0);

		localHermiteInterpolator(t1,xAuxRoot);
		eventFucntion(t1,xAuxRoot,&ef1);

		localHermiteInterpolator(tm,xAuxRoot);
		eventFucntion(tm,xAuxRoot,&efm);

		assert(ef0*ef1 < 0.0);

		if(ef0 >= ef1){
			eta2=  etam + (etam - eta0) * efm/sqrt(efm*efm - ef0*ef1);
		}else
			eta2=  etam - (etam - eta0) * efm/sqrt(efm*efm - ef0*ef1);


		t2= localTI + eta2*timeStep;

		localHermiteInterpolator(t2,x);
		eventFucntion(t2,x,&ef2);

//		printf ( "eta2= %+14.7E   ef= %+14.7E counter = %i\n",eta2,ef2,counter );
		
		assert(eta0 < eta2 && eta2 < eta1);
//		efh[counter]=ef2;

		assert(counter < 10000);
		if(fabs(ef2) < zeroTol){
			*time= localTI+eta2*timeStep;
//			printf ( "\ncouter = %i\n",counter );
			return;
		}

		if(ef0*ef2 <  0.0)
			eta1= eta2;
		else
			eta0= eta2;
		etam=(eta0 + eta1)*0.5;

		counter++;
		
//		if(counter > 100){
//			for(int i=0;i<100;i++)
//				printf ( "ef[%03i] = %+14.7E\n",i,efh[i] );
//			assert(0);
//		}

	}while(1);

}
