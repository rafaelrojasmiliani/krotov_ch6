/*!
 *       \file:  gslodeeventstep.cpp
 *
 *		 \brief: DESCRIPTION
 *        Version:  1.0
 *        Created:  19/06/15 10:04:14
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include"gslodeevent.h"
#include"indexvector.h"
#include<math.h>
#include<string.h> //memcpy 
void cGslOdeEvent::stepEventDectector(double *time,double x0in[]){
	double tAux;
	static double ev0[10],ev1[10];
	static double xAux[30];
	static int rres;
	tIndexvector status0,status1;


	localTI=*time;

	
	memcpy((void*)localXI,(void*)x0in,dim*sizeof(double));
	memcpy((void*)&localTI,(void*)time,1*sizeof(double));
	eventFucntion(*time,x0in,ev0);
//	printf ( "ev t0= %+14.7E  ev0 = %+14.7E  %+14.7E\n",*time,ev0[0],ev0[1] );
//	status0=

//	printf ( "entering step\n" );
	gsl_odeiv2_driver_apply(gslOdeDriver,time,*time+timeStep,x0in);
//	printf ( "step done\n" );

	eventFucntion(*time,x0in,ev1);
//	printf ( "ev t1= %+14.7E  ev1 = %+14.7E  %+14.7E\n",*time,ev1[0],ev1[1] );
	
	memcpy((void*)localXF,(void*)x0in,dim*sizeof(double));
	memcpy((void*)&localTF,(void*)time,1*sizeof(double));
	localInitHermiteInterpolator(localTI,localXI,localTF,localXF);

	for(size_t uici=0;uici < edim;uici++){
		if(fabs(ev0[uici]) > zeroTol && ev0[uici]*ev1[uici] < 0.0){
	//		printf ( "\tenterin root finder\n" );
			rres=localRootFindBrent(&tAux,xAux,uici);
			if(rres==1 &&  tAux< *time){
				*time=tAux;
				memcpy((void*)x0in,(void*)xAux,dim*sizeof(double));
			}
	//		printf ( "\t done\n" );
		}
	}
}
