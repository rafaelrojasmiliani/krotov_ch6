/*!
 *       \file:  affinesyscontrol.cpp
 *
 *		 \brief: DESCRIPTION
 *        Version:  1.0
 *        Created:  23/07/15 21:44:04
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include"affinesys.h"
#include<eigen3/Eigen/Core>
#include<eigen3/Eigen/Eigenvalues>
#include<complex>
#include<iostream>
using namespace Eigen;

void affinesys::computeControl(double t,const double xin[],
				const double umin[],const double umax[],double uout[],double zt){
	typedef Matrix<double,Dynamic,Dynamic,RowMajor> tMdXX;
	typedef Matrix<std::complex<double>,Dynamic,Dynamic,RowMajor> tMcdXX;
	static VectorXd Sv(affineDim);
	static tMdXX  A0	(affineDim,affineDim);
	static tMdXX  Aaux	(affineDim,affineDim);
	static tMdXX  A1	(affineDim,affineDim);
	static tMcdXX A2	(affineDim,affineDim);
	static tMcdXX Ident	(affineDim,affineDim);
	static tMdXX  A3	(affineDim,affineDim);
	static tMdXX  Ab	(affineDim,affineDim);
	static tMcdXX  Vminus(affineDim,affineDim);
	static tMcdXX  P(affineDim,affineDim);
	static ComplexEigenSolver<tMcdXX> eigenSolver;
	static VectorXd b0(affineDim);
	static VectorXd b1(affineDim);
	static VectorXd b3(affineDim);
	static VectorXd ueq(affineDim);

	static tIndexvector zeroSv;
	static tIndexvector nonZeroSv;
	static tIndexvector negEv;
	static tIndexvector nessearyCondU;
	static tIndexvector slidingControls;
	size_t uici,uicj,uick;

	convergetManifolds.clear();
	zeroSv.clear();
	negEv.clear();

	computeSwichingVector(t,xin,Sv.data());

	ueq.resize(affineDim);
	// store the indexes of Sv wich are zero in zeroSv
	for(uici=0;uici<affineDim;uici++)
		if(fabs(Sv(uici)) < zt){
			zeroSv.push_back(uici);
			ueq(uici)=0.0;
		}else{
			uout[uici]=(Sv(uici)>0.0)?umax[uici]:umin[uici];
			ueq(uici)=uout[uici];
			nonZeroSv.push_back(uici);
		}

	if(zeroSv.size() == 0) return;

	swichingVector_dt_LinSys(t,xin,b0.data(),A0.data());

//	std::cout << "---matrix A \n"<< A0 << "\n";
//	std::cout << "---vector b \n"<< b0 << "\n";
	b0+=A0*ueq; // compute the part of b associated to
				// the values of Sv that are not zero 
	
//	zeroSv.print("zeroSV = ");
	zeroSv.removeRows(Aaux,A0);
	zeroSv.removeRows(b1,b0);
	zeroSv.removeCols(A1,Aaux);
	
//	std::cout << "---matrix A for zeros\n"<< A1 << "\n";
//	std::cout << "---vector b for zeros\n"<< b0 << "\n";

// obtain the matrix corresponding to the
// values of Sv wich are zero.
	A2=A1.cast<std::complex<double> >();


	eigenSolver.compute(A2);

// copute wich eigven values have negative real part
	for(uici=0;uici<zeroSv.size();uici++)
		if(eigenSolver.eigenvalues()[uici].real()<0.0)
			negEv.push_back(uici);

//	negEv.print("negative ev in A ");
// compute the eigen space of eigenvalues with negative real part
	Vminus.resize(zeroSv.size(),negEv.size());
	for(uici=0;uici<negEv.size();uici++)
		Vminus.col(uici)=eigenSolver.eigenvectors().col(negEv[uici]);
//	std::cout << "matrix Vminus\n" << Vminus << "\n";
	P.resize(zeroSv.size(),zeroSv.size());
	P.setIdentity();
// the rank of Vmin is the space where the systems Sdot = A u  + b may be asymptotically stable
// and P is the projector in the orthogonal complement of Vminus.
	P-=Vminus*Vminus.adjoint();
//	std::cout << "Vminus projector orthogonal\n" << P << "\n";
	Ident.resize(zeroSv.size(),zeroSv.size());
	Ident.setIdentity();

	nessearyCondU.clear();

	for(uici=0;uici<negEv.size();uici++){
		static double val;
		val=P.col(uici).norm();
//		printf ( "val = %+14.7E\n",val );
		if(fabs(val)<1.0e-9)
			nessearyCondU.push_back(negEv[uici]);
	}
//	nessearyCondU.print("nc : ");

	if(nessearyCondU.size()== 0) return;


	nessearyCondU.removeRows(Aaux,A1);
	nessearyCondU.removeCols(A3,Aaux);
	nessearyCondU.removeRows(b3,b1);

//	std::cout << "Aeq \n "<< A3 << "\n";
//	std::cout << "beq \n "<< b3 << "\n";
//	nessearyCondU.print("nc: ","      \n");
	ueq.resize(nessearyCondU.size());
	ueq=-A3.partialPivLu().solve(b3);
//	std::cout << "ueq \n "<< ueq << "\n";
	for(uici=0;uici<affineDim;uici++){
		if(	zeroSv.isIn(uici,&uicj) && 
		negEv.isIn(uicj,&uick)  &&
		ueq(uick) <= umax[uici] &&
		ueq(uick) >= umin[uici])			
			uout[uici]=ueq(uick);
	}
//	for(uici=0;uici < nessearyCondU.size();uici++)
//		if(umin[nessearyCondU[uici]] <= ueq(uici) &&
//			ueq(uici) <= umax[nessearyCondU[uici]]){
//				convergetManifolds.push_back(nessearyCondU[uici]);
//				uout[nessearyCondU[uici]]=ueq(uici);
//		}
				
}
void affinesys::compute_1Control(double t,const double xin[],
				const double umin[],const double umax[],double uout[],
				double zt){
	double S; 
	double b,a;
	double g1,g2;
			
	convergetManifolds.clear();


	computeSwichingVector(t,xin,&S);

	if(fabs(S) > zt){
		uout[0]=(S > 0)? umax[0] : umin[0];
		return;
	}
	swichingVector_dt_LinSys(t,xin,&b,&a);

	g1= b+a*umin[0];
	g2= b+a*umax[0];

	if(g1 > 0.0 && g2 > 0.0) uout[0]=umax[0];
	else
	if(g1 < 0.0 && g2 < 0.0) uout[0]=umin[0];
	else
	if(g1 > 0.0 && g2 < 0.0){
		uout[0]= -b/a;
		convergetManifolds.push_back(0);
	}
}

void affinesys::compute_2Control(double t,const double xin[],
				const double umin[],const double umax[],double uout[],
				double zt){
		static double S[2];
		static double Aa[4];
		static double ba[2];
		static double ueqa[2];
		static double eigval1,eigval2;
		static double det,tr;
		static double g1,g2;
		static Map<Matrix<double,2,2,RowMajor> > A(Aa);
		static Matrix2cd Ac;
		static ComplexEigenSolver<Matrix2cd> ces;
		static Map<Matrix<double,2,1> > b(ba);
		static Map<Matrix<double,2,1> > ueq(ueqa);

        convergetManifolds.clear();

		computeSwichingVector(t,xin,S);

		if(fabs(S[0]) > zt && fabs(S[1]) > zt ){
			uout[0] = (S[0] > 0.0)? umax[0] : umin[0];
			uout[1] = (S[1] > 0.0)? umax[1] : umin[1];
			return;
		}


		swichingVector_dt_LinSys(t,xin,ba,Aa);


		if(fabs(S[0]) > zt && fabs(S[1]) < zt){
			uout[0] = (S[0] > 0.0)? umax[0] : umin[0];

			g1= b(1)+A(1,1)*umin[1] + A(1,0)*uout[0];		
			g2= b(1)+A(1,1)*umax[1] + A(1,0)*uout[0];		
			if(g1 > 0.0 && g2 > 0.0) uout[1]=umax[1];
			else
			if(g1 < 0.0 && g2 < 0.0) uout[1]=umin[1];
			else
			if(g1 > 0.0 && g2 < 0.0){
//				printf ( "b  = %14.7E \n", b(1));
//				printf ( "beq= %14.7E \n", (b(1)+A(1,0)*uout[0]));
//				printf ( "Aeq= %14.7E \n", A(1,1));
				uout[1]= -(b(1)+A(1,0)*uout[0])/A(1,1);
				convergetManifolds.push_back(1);
			}
			return;
		}else
		if(fabs(S[1]) > zt && fabs(S[0]) < zt){
//				printf ( "Sv(0)=0 Sv(1) != 0\n" );
			uout[1] = (S[1] > 0.0)? umax[1] : umin[1];

			g1= b(0)+A(0,0)*umin[0] + A(0,1)*uout[1];		
			g2= b(0)+A(0,0)*umax[0] + A(0,1)*uout[1];		
//	printf ( "g1  = %14.7E g2  = %14.7E \n", g1,g2);
			if(g1 > 0.0 && g2 > 0.0) uout[0]=umax[0];
			else
			if(g1 < 0.0 && g2 < 0.0) uout[0]=umin[0];
			else
			if(g1 > 0.0 && g2 < 0.0){
//				printf ( "b  = %14.7E \n", b(0));
//				printf ( "beq= %14.7E \n", (b(0)+A(0,1)*uout[1]));
//				printf ( "Aeq= %14.7E \n", A(0,0));
				uout[0]= -(b(0)+A(0,1)*uout[1])/A(0,0);
				convergetManifolds.push_back(0);
//				printf ( "c2: A(0,0) = %+14.7E\n",A(0,0) );
			}
			return;
		}
//		printf ( "%+14.7E    %+14.7E\n",S[0],S[1] );

		Ac=A.cast<std::complex<double> >();

		ces.compute(Ac,false);

		eigval1 = ces.eigenvalues()[0].real();
		eigval2 = ces.eigenvalues()[1].real();
//		printf ( "real part of eigen values %+14.7E  %+14.7E\n",eigval1,eigval2 );

		if(eigval1 < -1.0e-18 && eigval2 < -1.0e-18){
//				printf ( "solving\n" );
			ueq=-A.partialPivLu().solve(b);

//			printf ( " t= %+14.7E  we get ueq = %+14.7E  %+14.7E\n",t,ueq(0),ueq(1) );

			if(	umin[0] <= ueq(0) && ueq(0) <= umax[0] &&
				umin[1] <= ueq(1) && ueq(1) <= umax[1] ){
				uout[0]=ueq(0);
				uout[1]=ueq(1);
				convergetManifolds.push_back(0);
				convergetManifolds.push_back(1);
//				printf ( "the two!!!\n" );
				return;
			}else
			if(	umin[0] <= ueq(0) && ueq(0) <= umax[0] && fabs(A(0,1)) < 1.0e-18 ){
				uout[0]= ueq(0);
				g1= b(1)+A(1,1)*umin[1] + A(1,0)*uout[0];		
				g2= b(1)+A(1,1)*umax[1] + A(1,0)*uout[0];		
				if(g1 > 0.0 && g2 > 0.0) uout[1]=umax[1];
				else
				if(g1 < 0.0 && g2 < 0.0) uout[1]=umin[1];
//				printf ( "just 0\n" );
				return;
			}else
			if(	umin[1] <= ueq(1) && ueq(1) <= umax[1] && fabs(A(1,0)) < 1.0e-18 ){
				uout[1]= ueq(1);
				g1= b(0)+A(0,0)*umin[0] + A(0,1)*uout[1];		
				g2= b(0)+A(0,0)*umax[0] + A(0,1)*uout[1];		
				if(g1 > 0.0 && g2 > 0.0) uout[0]=umax[0];
				else
				if(g1 < 0.0 && g2 < 0.0) uout[0]=umin[0];
//				printf ( "just 1\n" );
				return;
			}else{

			}
		}
}
