/*!
 *       \file:  gslodeeventrootbrent.cpp
 *
 *		 \brief: DESCRIPTION
 *        Version:  1.0
 *        Created:  21/06/15 12:53:58
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

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define ITMAX 100
#define EPS (zeroTol)

double tHist[ITMAX];
double fHist[ITMAX];

int cGslOdeEvent::localRootFindBrent(double *time,double x[]){
	
	int iter;
	double a=localTI,b=localTF,c=localTF,d,e=100.0,min1,min2;
	double fc,p,q,r,s,tol1,xm;
	double fa,fb;
	eventFucntion(localTI,localXI,&fa);
	eventFucntion(localTF,localXF,&fb);

	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
		assert(0);
	fc=fb;
	for (iter=0;iter<ITMAX;iter++) {
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
			c=a;
	//		Rename a, b, c and adjust bounding intervald.
			fc=fa;
			e=d=b-a;
		}
		if (fabs(fc) < fabs(fb)) {
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}
		tol1=1.0e-23; //Convergence check.
		xm=0.5*(c-b);
		if (fabs(fb) < zeroTol ){
//			printf (" the interval is less than"
//					" zero tolerance or one of its"
//					" liimits is exactly zero\n" );
			*time=b;
//			printf ( "couter = %i\n",iter );
			return 1;
		}
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
			s=fb/fa; //Attempt inverse quadratic interpolation.
			if (a == c) {
				p=2.0*xm*s;
				q=1.0-s;
			} else {
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q; //Check whether in bounds.
			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e=d;	//Accept interpolation.
				d=p/q;
			} else {
				d=xm;//Interpolation failed, use bisection.
				e=d;
			}
		} else{ //Bounds decreasing too slowly, use bisection.
			d=xm;
			e=d;
		}
		a=b;
		fa=fb;
		if (fabs(d) > tol1)
			b += d;
		else
			b += SIGN(fabs(xm)*1.0e-4,xm);


		localHermiteInterpolator(b,x);
//		printf ( "computing wf:\n" );
		eventFucntion(b,x,&fb);
//		printf ( "done.\n" );
		tHist[iter]= b;
		fHist[iter]= fb;
	}
	printf ( "[gslEvent:root] Brent Method failed\n" );
//	for(iter=0;iter < ITMAX; iter++)
//		printf ( "t = %+14.7E \t f = %+14.7E\n",tHist[iter],fHist[iter] );
//	assert(0);
	return 0;
//	printf ( "interval = [ %14.7E  %+14.7E ]\n",localTI,localTF );
//
////	nrerror("Maximum number of iterations exceeded in zbrent");

}
int cGslOdeEvent::localRootFindBrent(double *time,double x[],size_t i){
	
	int iter;
	double a=localTI,b=localTF,c=localTF,d,e=100.0,min1,min2;
	static double fc,p,q,r,s,tol1,xm;
	static double fa,fb;
	static double aux[10];
	eventFucntion(localTI,localXI,aux);
	fa=aux[i];
//	printf ( "...................................\n" );
//	printf ( "ev t0= %+14.7E  ev0 = %+14.7E  %+14.7E   (i  = %zu )\n",localTI,aux[0],aux[1],i );
	eventFucntion(localTF,localXF,aux);
	fb=aux[i];
//	printf ( "ev t1= %+14.7E  ev1 = %+14.7E  %+14.7E   (i  = %zu )\n",localTF,aux[0],aux[1],i );

	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)){
			printf ( "[][gslEvent:root] fatal error, "
							"there is not a root\n" );
		return 0;
	}
	fc=fb;
	for (iter=0;iter< ITMAX;iter++) {
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
			c=a;
	//		Rename a, b, c and adjust bounding intervald.
			fc=fa;
			e=d=b-a;
		}
		if (fabs(fc) < fabs(fb)) {
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}
		tol1=1.0e-13; //Convergence check.
		xm=0.5*(c-b);
		if (fabs(fb) < zeroTol && fabs(c-b) < zeroTol){
//			printf (" the interval is less than"
//					" zero tolerance or one of its"
//					" liimits is exactly zero\n" );
			*time=b;
//			printf ( "couter = %i\n",iter );
			return 1;
		}
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
			s=fb/fa; //Attempt inverse quadratic interpolation.
			if (a == c) {
				p=2.0*xm*s;
				q=1.0-s;
			} else {
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q; //Check whether in bounds.
			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e=d;	//Accept interpolation.
				d=p/q;
			} else {
				d=xm;//Interpolation failed, use bisection.
				e=d;
			}
		} else{ //Bounds decreasing too slowly, use bisection.
			d=xm;
			e=d;
		}
		a=b;
		fa=fb;
		if (fabs(d) > tol1)
			b += d;
		else
			b += SIGN(fabs(xm)*1.0e-10,xm);


		localHermiteInterpolator(b,x);
		eventFucntion(b,x,aux);
		fb=aux[i];
		tHist[iter]= b;
		fHist[iter]= fb;
	}
	printf ( "[gslEvent:root] Brent Method failed in the interval %+14.7E  %+14.7E\n",localTI,localTF );
//	for(iter=0;iter < ITMAX; iter++)
//		printf ( "t = %.20E \t f = %+.20E\n",tHist[iter],fHist[iter] );
//	assert(0);
	return 0;
//	printf ( "interval = [ %14.7E  %+14.7E ]\n",localTI,localTF );
//
//	for(iter=0;iter < ITMAX; iter++)
//		printf ( "t = %+14.7E \t f = %+14.7E\n",tHist[iter],fHist[iter] );
////	nrerror("Maximum number of iterations exceeded in zbrent");
//	assert(0);

}
//int cGslOdeEvent::localBrentStep(double ts[3],double fts[3]){
//	
//	static double d,e=100.0,min1,min2;
//	static double p,q,r,s,tol1,xm;
//
////	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
////		return -1
////	fc=fb;
////	for (iter=1;iter<=ITMAX;iter++) {
//		if ((fts[1] > 0.0 && fts[2] > 0.0) || (fts[1] < 0.0 && fts[2] < 0.0)) {
//			ts[2]=ts[0];
//	//		Rename a, b, c and adjust bounding intervald.
//			fts[2]=fts[0];
//			e=d=ts[1]-ts[0];
//		}
//		if (fabs(fts[2]) < fabs(fts[1])) {
//			ts[0]=ts[1];
//			ts[1]=ts[2];
//			ts[2]=ts[0];
//			fts[0]=fts[1];
//			fts[1]=fts[2];
//			fts[2]=fts[0];
//		}
//		tol1=1.0e-13; //Convergence check.
//		xm=0.5*(ts[2]-ts[1]);
//		if (fabs(fts[1]) < zeroTol ){
////			printf (" the interval is less than"
////					" zero tolerance or one of its"
////					" liimits is exactly zero\n" );
//			*time=ts[1];
////			printf ( "couter = %i\n",iter );
//			return 1;
//		}
//		if (fabs(e) >= tol1 && fabs(fts[0]) > fabs(fts[1])) {
//			s=fts[1]/fts[0]; //Attempt inverse quadratic interpolation.
//			if (ts[0] == ts[2]) {
//				p=2.0*xm*s;
//				q=1.0-s;
//			} else {
//				q=fts[0]/fts[2];
//				r=fts[1]/fts[2];
//				p=s*(2.0*xm*q*(q-r)-(ts[1]-ts[0])*(r-1.0));
//				q=(q-1.0)*(r-1.0)*(s-1.0);
//			}
//			if (p > 0.0) q = -q; //Check whether in bounds.
//			p=fabs(p);
//			min1=3.0*xm*q-fabs(tol1*q);
//			min2=fabs(e*q);
//			if (2.0*p < (min1 < min2 ? min1 : min2)) {
//				e=d;	//Accept interpolation.
//				d=p/q;
//			} else {
//				d=xm;//Interpolation failed, use bisection.
//				e=d;
//			}
//		} else{ //Bounds decreasing too slowly, use bisection.
//			d=xm;
//			e=d;
//		}
//		ts[0]=ts[1];
//		fts[0]=fts[1];
//		if (fabs(d) > tol1)
//			ts[1] += d;
//		else
//			ts[1] += SIGN(1.0e-20,xm);
//
//	return 0;
//
////	}
////	printf ( "interval = [ %14.7E  %+14.7E ]\n",localTI,localTF );
////
////	for(iter=0;iter < ITMAX; iter++)
////		printf ( "t = %+14.7E \t f = %+14.7E\n",tHist[iter],fHist[iter] );
//////	nrerror("Maximum number of iterations exceeded in zbrent");
////	assert(0);
//
//}
