/*!
 *       \file:  affinesys.cpp
 *
 *		 \brief: Base class to represent affine systems. 
 *        Created:  19/06/15 09:58:16
 *
 *         Author:  Rafael A. Rojas 
 *
 * =====================================================================================
 */

#include<eigen3/Eigen/Core>
#include"indexvector.h"
#include"affinesys.h"
#define MatrixXdRowMajor Matrix<double,Dynamic,Dynamic,RowMajor>
using namespace Eigen;
void affinesys::computeSwichingVector(double t,
									const double xin[],double *res,int flag){
// b->res, l->lambda, h->sigmaMem2 {{{
			const double* const csptr= getCostate(t);
			const double* const x0ptr= getState0(t);
			Map<const VectorXd> x(xin,stateDim) ;
			Map<const VectorXd> x0(x0ptr,stateDim) ;
			Map<const VectorXd> psi(csptr,stateDim) ;

			Map<VectorXd> b(res,affineDim) ;
			static Map<VectorXd> l(auxMem1,stateDim) ;

 			static Map< MatrixXdRowMajor > h(sigmaMem2,stateDim,affineDim);
 			static Map< RowVectorXd > h0(sigmaMem2,affineDim);

			sym2full_cpy(sigmaMem1,csptr+stateDim,stateDim);
			// }}}
			
//			if(flag ==1 ){
//					printf ( "[csf] t = %+14.7E\n",t );
//				iterator::processR->printActualInterval("[csf] x0 ");
//				iterator::costateCont.printActualInterval("[csf] cs ");
//				std::cout << "[csf] cs: "<<Map<const RowVectorXd> (csptr,costateDim) << "\n";
//				std::cout << "[csf] x0: "<<Map<const RowVectorXd> (x0ptr,stateDim) << "\n";
//			}
			l = psi+ sigma*(x - x0 );
//			projectonSF0Tangent(t,xin,l.data());
//			printf ( "\t\t\t\t\t\t\t(%+.5E) l = %+14.7E ="
//					": %+14.7E"
//					" + %+14.7E  +%+14.7E + %+14.7E+ %+14.7E \n",
//					t,l(3),psi(3),
//					sigma(3,0)*(x(0)-x0(0)),sigma(3,1)*(x(1)-x0(1)),
//					sigma(3,2)*(x(2)-x0(2)),sigma(3,3)*(x(3)-x0(3)));
			dynSysLinearPart(t,xin,h.data());
			b= h.transpose()*l;
			if(flag ==1 ){
					printf ( "[csf] h^T*lambda = %+14.7E\n",b(0) );
			}


			runningCostLinearPart(t,xin,h0.data());
			if(flag ==1 ){
					printf ( "[csf] h0         = %+14.7E\n",h0(0) );
			}

			b-= h0.transpose();
		}

void affinesys::computeSwichingVector(double t,
									const double xin[],const double csptr[],
									const double x0ptr[],double *res,int flag){
// b->res, l->lambda, h->sigmaMem2 {{{
			Map<const VectorXd> x(xin,stateDim) ;
			Map<const VectorXd> x0(x0ptr,stateDim) ;
			Map<const VectorXd> psi(csptr,stateDim) ;

			Map<VectorXd> b(res,affineDim) ;
			static Map<VectorXd> l(auxMem1,stateDim) ;

 			static Map< MatrixXdRowMajor > h(sigmaMem2,stateDim,affineDim);
 			static Map< RowVectorXd > h0(sigmaMem2,affineDim);

			sym2full_cpy(sigmaMem1,csptr+stateDim,stateDim);
			// }}}
			
//			if(flag ==1 ){
//					printf ( "[csf] t = %+14.7E\n",t );
//				iterator::processR->printActualInterval("[csf] x0 ");
//				iterator::costateCont.printActualInterval("[csf] cs ");
//				std::cout << "[csf] cs: "<<Map<const RowVectorXd> (csptr,costateDim) << "\n";
//				std::cout << "[csf] x0: "<<Map<const RowVectorXd> (x0ptr,stateDim) << "\n";
//			}
			l = psi+ sigma*(x - x0 );
//			projectonSF0Tangent(t,xin,l.data());
//			printf ( "\t\t\t\t\t\t\t(%+.5E) l = %+14.7E ="
//					": %+14.7E"
//					" + %+14.7E  +%+14.7E + %+14.7E+ %+14.7E \n",
//					t,l(3),psi(3),
//					sigma(3,0)*(x(0)-x0(0)),sigma(3,1)*(x(1)-x0(1)),
//					sigma(3,2)*(x(2)-x0(2)),sigma(3,3)*(x(3)-x0(3)));
			dynSysLinearPart(t,xin,h.data());
			b= h.transpose()*l;
			if(flag ==1 ){
					printf ( "[csf] h^T*lambda = %+14.7E\n",b(0) );
			}


			runningCostLinearPart(t,xin,h0.data());
			if(flag ==1 ){
					printf ( "[csf] h0         = %+14.7E\n",h0(0) );
			}

			b-= h0.transpose();
		}
/*! Computes the derivative of the linear system corresponding to the switching vector w.r.t. the time.  
 * \param t, time
 * \param xin state
 * \param *bout, derivate of the rhs of the linear system
 * \param *Aout, derivate of the matrix in the lhs of the linear system
*/  
void affinesys::swichingVector_dt_LinSys(const double t,
									const double xin[],
									double *bout,double *Aout){
			
// initialization b,A,x,x0,psi,sigma{{{
	Map<const VectorXd > x(xin,stateDim);
	Map< MatrixXdRowMajor > A(Aout,affineDim,affineDim);
	Map<VectorXd> b(bout,affineDim) ;

	const double* const csptr=getCostate(t);
	Map<const VectorXd> psi(csptr,stateDim) ;
	sym2full_cpy(sigmaMem1,csptr+stateDim,stateDim);
	const double* const x0ptr=getState0(t);
	Map<const VectorXd> x0(x0ptr,stateDim);
	const double* const u0ptr=getControl0(t);

	static Map<VectorXd> l(auxMem1,stateDim) ;
		static Map< MatrixXdRowMajor > h(sigmaMem2,stateDim,affineDim);
	static VectorXd g(stateDim);
	static VectorXd f0(stateDim);
	static size_t uici,uicj,uick;
// }}}

//std::cout <<  t << " " << Map<const RowVectorXd>(csptr,costateDim) << "\n";

	double b1,b2,b3,a1,a2,a3;
	dynSysCostate(t,csptr,auxMem1);
	Map<const VectorXd> psi_dot(auxMem1,stateDim) ;
	sym2full_cpy(sigmaMem2,auxMem1+stateDim,stateDim);

	dynSysAffinePart(t,xin,iterator::control,g.data());
	dynSystem(t,x0ptr,u0ptr,f0.data());

	l=-psi_dot-sigma_dot*(x-x0)+sigma*(g-f0);

	

	dynSysLinearPart(t,xin,h.data());

	b= h.transpose()*l; // --- s dot operation (b_1)
//			DEBUG_b1=b(0);

			A= h.transpose()*sigma*h; // --- s dot operation (A_1)
//			DEBUG_a1=A(0,0);

			l= psi+sigma*(x-x0);

			dynSysLinearPart_1G(t,xin,jacMem);// ->jacMem

//  m1 to sigmaMem2 as affineDim x stateDim and {{{
// computation of the contraction in i of
//		m1^k_l=l_i* {dh^i_k \over x^l}
	new (&m1) Map< MatrixXdRowMajor >(sigmaMem1,affineDim,stateDim);

static const size_t stride=affineDim*stateDim;
#define dhikdxl(i,k,l) (*(jacMem+stride*(l)+affineDim*(i)+(k)) )
	for(int uick=0;uick<affineDim;uick++){
		for(int uicl=0;uicl<stateDim;uicl++){
			m1(uick,uicl)=0.0;
			for(int uici=0;uici<stateDim;uici++)
				m1(uick,uicl)+=dhikdxl(uici,uick,uicl)*l(uici);
		}
	}
#undef dhikdxl
//		}}}
		
	A+= m1*h; // s dot operation  A2
	b+= m1*g; 
//			DEBUG_b2a=b(0)-DEBUG_b1;
//			DEBUG_a2=A(0,0)-DEBUG_a1;
	new (&m1) 
		Map< Matrix<double,Dynamic,Dynamic,RowMajor> >
												(jacMem,stateDim,affineDim);
	dynSysLinearPart_dt(t,xin,m1.data());
	b+= m1.transpose()*l ;
//			DEBUG_b2b=b(0)-DEBUG_b1-DEBUG_b2a;


	new (&m1) 
				Map< Matrix<double,Dynamic,Dynamic,RowMajor> >
												(sigmaMem1,stateDim,affineDim);
	runningCostLinearPart_1G(t,xin,m1.data());

	A-= m1.transpose()*h; // s dot operation (A3)

	b-= m1.transpose()*g; // s dot operation (b3.1)

	new (&m1) 
		Map< Matrix<double,Dynamic,Dynamic,RowMajor> >
												(sigmaMem1,1,affineDim);
	runningCostLinearPart_dt(t,xin,m1.data());

	b-=  m1.transpose();	// s dot operation (b3.2)
}


	// This method retur a matrix m(i,j)= d S_i / d x^j
void affinesys::swichingVector_G1(double t,
									const double xin[],double *res){

// initialization skxl,x,x0,psi,sigma,h{{{
	static int flag=1;
			Map< MatrixXdRowMajor> skxl(res,affineDim,stateDim);
			Map< const VectorXd > x(xin,stateDim) ;

			const double* const csptr=getCostate(t);
			const double* const x0ptr=getState0(t);

			Map<const VectorXd> x0(x0ptr,stateDim);
			Map<const VectorXd> psi(csptr,stateDim) ;
			sym2full_cpy(sigmaMem1,csptr+stateDim,stateDim);

			static Map<VectorXd> l(auxMem1,stateDim) ;
 			static Map< MatrixXdRowMajor > h(sigmaMem2,stateDim,affineDim);

			new (&m1) 
				Map<  Matrix<double,Dynamic,Dynamic,RowMajor> >
												(sigmaMem1,affineDim,stateDim);
//// }}}
		
			dynSysLinearPart(t,xin,h.data());
			skxl=h.transpose()*sigma;

			l=psi+sigma*(x-x0);

//	if( t > 0.149 && t < 0.153){
//			printf ( "[dsdx] h sima (0,0)  %+14.7E\n",skxl(0,0) );
//	}
			dynSysLinearPart_1G(t,xin,jacMem);
			static const int stride=stateDim*affineDim;
/** m1(k,l)= l_i * ( dh^i_k / dx^l) *  */
#define dhikdxl(i,k,l) (*(jacMem+stride*(l)+affineDim*(i)+(k)) )
			for(int uicl=0;uicl<stateDim;uicl++){
				for(int uick=0;uick<affineDim;uick++){
					m1(uick,uicl)=0.0;
					for(int uici=0;uici<stateDim;uici++)
						m1(uick,uicl)+=dhikdxl(uici,uick,uicl)*l(uici);
				}
			}
#undef dhikdxl
			skxl+=m1;

			runningCostLinearPart_1G(t,xin,h.data());
			skxl-=h.transpose();

}
void affinesys::swichingVector_G1(double t,
									const double xin[],
									const double csptr[],
									const double x0ptr[],double *res){

// initialization skxl,x,x0,psi,sigma,h{{{
	static int flag=1;
			Map< MatrixXdRowMajor> skxl(res,affineDim,stateDim);
			Map< const VectorXd > x(xin,stateDim) ;

			Map<const VectorXd> x0(x0ptr,stateDim);
			Map<const VectorXd> psi(csptr,stateDim) ;
			sym2full_cpy(sigmaMem1,csptr+stateDim,stateDim);

			static Map<VectorXd> l(auxMem1,stateDim) ;
 			static Map< MatrixXdRowMajor > h(sigmaMem2,stateDim,affineDim);

			new (&m1) 
				Map<  Matrix<double,Dynamic,Dynamic,RowMajor> >
												(sigmaMem1,affineDim,stateDim);
//// }}}
		
			dynSysLinearPart(t,xin,h.data());
			skxl=h.transpose()*sigma;

			l=psi+sigma*(x-x0);

//	if( t > 0.149 && t < 0.153){
//			printf ( "[dsdx] h sima (0,0)  %+14.7E\n",skxl(0,0) );
//	}
			dynSysLinearPart_1G(t,xin,jacMem);
			static const int stride=stateDim*affineDim;
/** m1(k,l)= l_i * ( dh^i_k / dx^l) *  */
#define dhikdxl(i,k,l) (*(jacMem+stride*(l)+affineDim*(i)+(k)) )
			for(int uicl=0;uicl<stateDim;uicl++){
				for(int uick=0;uick<affineDim;uick++){
					m1(uick,uicl)=0.0;
					for(int uici=0;uici<stateDim;uici++)
						m1(uick,uicl)+=dhikdxl(uici,uick,uicl)*l(uici);
				}
			}
#undef dhikdxl
//			std::cout << "m1= [" << m1 << "]\n";
//
//	if( t > 0.149 && t < 0.153){
//			printf ( "[dsdx] contrac (0,0) %+14.7E\n",m1(0,0) );
//	}
			skxl+=m1;
//		std::cout << "2dsdx= [" << skxl << " ]\n";

			runningCostLinearPart_1G(t,xin,h.data());
//			std::cout << "dh0dx \n" << h << "\n" ;
				
			skxl-=h.transpose();
//		std::cout << "3dsdx= [" << skxl << " ]\n";


		}

/*! Computes the partial derivative of the switching vector w.r.t. the time.  
 * \param t, time
 * \param xin state
 * \param res, partial derivative of the switching vector
*/  
void affinesys::swichingVector_partial_dt(double t,
									const double xin[],double *res){

// initialization skxl,x,x0,psi,sigma,h{{{
	typedef Matrix<double,Dynamic,Dynamic,RowMajor> RowMajorMatrix;
	static int flag=1;
	Map< VectorXd > dsdt(res,affineDim);
	Map< const VectorXd > x(xin,stateDim) ;

	const double* const csptr=getCostate(t);
	const double* const x0ptr=getState0(t);
	const double* const u0ptr=getControl0(t);

	Map<const VectorXd> x0(x0ptr,stateDim);
	Map<const VectorXd> psi(csptr,stateDim) ;
	sym2full_cpy(sigmaMem1,csptr+stateDim,stateDim);

	static RowMajorMatrix h(stateDim,affineDim);
	static RowMajorMatrix dhdt(stateDim,affineDim);
	static MatrixXd dh0dt(1,affineDim);

	static VectorXd f0(stateDim);

//// }}}
		
	dynSysCostate(t,csptr,auxMem1);
	Map<const VectorXd> psi_dot(auxMem1,stateDim) ;
	sym2full_cpy(sigmaMem2,auxMem1+stateDim,stateDim);

	dynSystem(t,x0ptr,u0ptr,f0.data());
	dynSysLinearPart(t,xin,h.data());

	dsdt=h.transpose()*(-psi_dot-sigma_dot*(x-x0)-sigma*f0);
	dynSysLinearPart_dt(t,xin,dhdt.data());
	runningCostLinearPart_dt(t,xin,dh0dt.data());
	dsdt+=dhdt.transpose()*(psi+sigma*(x-x0))-dh0dt.transpose();

}

/*! Test and print to stdout the Krotov's R function  
 * \param t, time
 * \param x state
 * \param u control
 * \param level, level of verbosity
*/  
void affinesys::testR(double t,const double x[],const double u[],int level){
		

	if(level==0) return;
	const double *x0,*u0,*cs;
	x0=getState0(t);
	u0=getControl0(t);
	cs=getCostate(t);

	int assertion=1;

	double rtxu,rtxu0,rtx0u0;

	rtxu=R(t,x,u);

	rtxu0=R(t,x,u0);

	rtx0u0=R(t,x0,u0);

	
	if(rtxu - rtxu0 <= -5.0e-5){
		printf ( "R(t,x,u) - R(t,x,u0) is not working\n" );
		printf ( "R(t,x ,u ) = %+14.7E\n",rtxu );
		printf ( "R(t,x ,u0) = %+14.7E\n",rtxu0 );
		printf ( "dR         = %+14.7E\n",rtxu-rtxu0 );
		double S;
		computeSwichingVector(t,x,&S);
		printf ( "sf         = %+14.7E\n",S );
		assert(0);
	}
	if(level==1) return;


	if(rtxu0 - rtx0u0 <= -5.0e-5){
		printf ( "R(t,x ,u0) - R(t,x0,u0) is not working\n" );
		printf ( "R(t,x ,u ) = %+14.7E\n",rtxu0 );
		printf ( "R(t,x0,u0) = %+14.7E\n",rtx0u0 );
		printf ( "dR         = %+14.7E\n",rtxu0-rtx0u0 );
			VectorXd x1=Map<const VectorXd>(x0,stateDim);
			VectorXd x2=Map<const VectorXd>(x,stateDim);
			VectorXd x3=Map<const VectorXd>(x,stateDim);
			printf ( "variation in direction from x0 to x\n" );
				for(double tt=0.0;tt<=1.1;tt+=0.1){
					x3=tt*(x2-x1) + x1;
					std::cout << x3.transpose() << " ";
					printf ( "%+14.7E  %+14.7E\n",tt,R(t,x3.data(),u0) );
				}
		assert(0);
	}

	if(level==2) return;


}

//double affinesys::R(double t,const double xin[],const double uin[]){
//
//	double result=0.0;
////			Map<const VectorXd> x0(getState0(t).data(),stateDim);
//	
//	const double* const csptr= getCostate(t);
//	const double* const x0ptr= getState0(t);
//	const double* const u0ptr= getControl0(t);
//	Map<const VectorXd> x0(x0ptr,stateDim);
//	Map<const VectorXd> psi(csptr,stateDim);
//	sym2full_cpy(sigmaMem1,csptr+stateDim,stateDim);
//	dynSysCostate(t,csptr,auxMem1);	
//	sym2full_cpy(sigmaMem2,auxMem1+stateDim,stateDim);
//	dynSystem(t,x0ptr,u0ptr,jacMem+stateDim);
//	dynSystem(t,xin,uin,jacMem);
//
//	static Map< VectorXd >
//				psi_dot(auxMem1,stateDim) ;
//	static Map< VectorXd > auxv(auxMem1,stateDim) ;
//
//	Map <const VectorXd> x(xin,stateDim);
//	Map <VectorXd> x_dot(jacMem,stateDim);
//	Map <VectorXd> x_dot0(jacMem+stateDim,stateDim);
//
////			std::cout << "x \n" << x.transpose() << "\n-----\n";
//	result-= psi_dot.transpose()*x;
////			std::cout << "psidot*x \n" << result << "\n-----\n";
//
////			std::cout << "sigma \n" << sigma << "\n-----\n";
//	auxv=(sigma*(x_dot-x_dot0)-0.5*sigma_dot*(x-x0));
//
////			printf ( "\nsigma*df - 0.5 sigma*dx\n" );
////			std::cout << auxv.transpose() << "\n-----\n";
//
//	result+=(x-x0).transpose()*auxv;
//
//	result+= psi.transpose()*x_dot;
//
//	result-=runningCost(t,xin,uin);
//
//
////	auxv= psi +sigma*(x_dot-x_dot0);
////			std::cout << "lm in R = " << auxv.transpose() << "\n";
//
////			printf ( "S in R= %+14.7E\n",S );
//	return result;
//	
//}
#undef MatrixXdRowMajor 
