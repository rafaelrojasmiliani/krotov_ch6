#include"affinesys.h"


#include<eigen3/Eigen/Core>
#include<eigen3/Eigen/LU>
#include<eigen3/Eigen/Cholesky>
using namespace Eigen;


void affinesys::projectonSF0(double t,double xin[]){
	static VectorXd sfv(affineDim),sfvred(affineDim);
	static VectorXd x0(stateDim),dx(stateDim),dlambda(affineDim),
								 bx(stateDim),blambda(affineDim),x0_ant(stateDim);
	static VectorXd lambda0(affineDim),lambda0_ant(affineDim);
	static Matrix<double,Dynamic,Dynamic,RowMajor> dsdx(affineDim,stateDim);
	static Matrix<double,Dynamic,Dynamic,RowMajor> dsdxred(affineDim,stateDim);
	static Matrix<double,Dynamic,Dynamic,RowMajor> A1(affineDim,affineDim);
	Map< VectorXd > xhat(xin,stateDim);
	static size_t subManifoldDim;
	tIndexvector svStatus1;
//	convergetManifolds.print( "----------------------computinf projection on ",
//				   "----------------------------------\n"	);
//std::cout << "--x0 : " << xhat.transpose() << "\n";
	if(convergetManifolds.empty()) return;
	computeSwichingVector(t,xhat.data(),sfv.data());//,( t > 0.149 && t < 0.153)?1:0);
//std::cout << "--sfv0 : " << sfv.transpose() << "\n";
//	svStatus1.print("actual zero of sfv:");
	convergetManifolds.removeRows(sfvred,sfv);
//	convergetManifolds.print("projecting into :");

//	if( t > 0.149 && t < 0.153){
//			printf ( "[proy] sf = %+14.7E\n",sfv(0) );
//	}
	if(sfvred.norm() < zeroTol) return;

	subManifoldDim=convergetManifolds.size();
	dlambda.resize(subManifoldDim);
	blambda.resize(subManifoldDim);
	dsdxred.resize(subManifoldDim,stateDim);
	A1.resize(subManifoldDim,subManifoldDim);
	lambda0.resize(subManifoldDim);

	x0=xhat;



	lambda0.setOnes();
	lambda0*=0.0005;


	static double sfhis[1000];
	static int counter;
	counter= 0;
	double localGain=localProyGain;
//	if( t > 0.15 && t < 0.153){
//		printf ( "= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = \n" );
//		printf ( "= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = \n" );
//		printf ( "= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = \n" );
//	}
//  std::cout << "----sfv = [ " << sfv.transpose() << "] \n";
//	printf ( "iterating ------------------------------\n" );
	while(1){
		x0_ant=x0;
		lambda0_ant=lambda0;

		swichingVector_G1(t,x0.data(),dsdx.data());
		convergetManifolds.removeRows(dsdxred,dsdx);

//	if( t > 0.149 && t < 0.153){
//std::cout << "-------x0= [ " << x0.transpose() << "\n";
//std::cout << "--lambda0= [ " << lambda0.transpose() << "]\n ";
//	}
		bx		= - localGain*((xhat-x0)+dsdxred.transpose()*lambda0);
		blambda	= -sfvred-dsdxred*bx;

		A1.noalias()= dsdxred*dsdxred.transpose();

		dlambda=- A1.llt().solve(blambda);
		dx= bx - dsdxred.transpose()*dlambda;

//std::cout << "  x - xhat= [ " << (x0-xhat).transpose() << "]\n ";
		x0+=dx;
		lambda0+= dlambda;

//	if( t > 0.149 && t < 0.153){
////{{{ 
//std::cout << "-----sfv = [ " << sfv.transpose() << "] \n";
//std::cout << "----xhat = [" << xhat.transpose() << "\n";
//std::cout << "----dsdx = [ " << dsdx << "\n-----\n";
//std::cout << "----  bx = [ " << bx.transpose()	<< "]\n";
//std::cout << "dSdx^t dl = [ " << (dsdxred.transpose()*dlambda).transpose()<< "]\n";
//std::cout << "--blambda= [ " << blambda.transpose() << "]\n";
//std::cout << "--localgain =" << localGain<< "\n";
//printf ( "    --      A= %+14.7E\n",A1(0,0) );																					
//std::cout << "\t\t\t\t dsdxred= [ \n" << dsdxred << "\n-----\n";
//std::cout << "\t\t\t\t      dx= [ " << dx.transpose() << "]\n ";
//std::cout << "\t\t\t\t      x0= [ " << x0.transpose() << "\n";
//std::cout << "\t\t\t\t dlambda= [ " << dlambda.transpose() << "]\n ";
//std::cout << "\t\t\t\t lambda0= [ " << lambda0.transpose() << "]\n ";
////}}}
//	}
		computeSwichingVector(t,x0.data(),sfv.data());

		convergetManifolds.removeRows(sfvred,sfv);
//		printf ( "time = %+14.7E\n",t );
//std::cout << "\n---sfv = [ " << sfv.transpose() << "] \n";
//		svStatus1.print("are zero: ",".------------------------\n");
//	convergetManifolds.print("convergetManifolds:");
		
		if(sfvred.norm() < proyTolSfVal && dx.norm() < proyTolDx)
			break;

		sfhis[counter]= sfv(0);	
		if(counter > 1 && 
			(fabs(sfhis[counter]) >= fabs(sfhis[counter-1]))){
			localGain/= 1.2;
			x0=x0_ant;
			lambda0=lambda0_ant;
		}
		counter++;
		if(counter>proyMaxIter){
			printf ( "[warning affine Proy]: "
			"filed on projector at time %+14.7E "
			"wth an error of %+14.7E\n",t,sfv(0) );
			convergetManifolds.print("[warning affine Proy]: "
							"failed projectin on :");
			std::cout << "[warning affine Proy]: xhat= ["
					<< xhat.transpose() << "\n";
			std::cout << "[warning affine Proy]: x0  = ["
					<< x0.transpose() << "\n";
			printf ( "[warning affine Proy] proyTolSfVal = %+14.7E \n",proyTolSfVal );
			printf ( "[warning affine Proy] proyTolDx = %+14.7E \n",proyTolDx );
				xhat=x0;
			return;
			for(int uici=0;uici < 100; uici++)
					printf ( "sf %i  = %+14.7E\n",uici,sfhis[uici] );
			assert(0);

		}
//		printf ( "sf_i= %+14.7E \n",sfv(0) );
	};
//	printf ( "counter projector  %i\n",counter );
//	printf ( "-------------------------------------------\n" );
//	printf ( "sf= %+14.7E \n",sf );

//	if( t > 0.149 && t < 0.153){
//		printf ( "[proy] xhat ="
//						"[ %+14.7E %+14.7E %+14.7E %+14.7E ]\n",
//					xhat(0),xhat(1),xhat(2),xhat(3) );
//		printf ( "[proy] x    ="
//						"[ %+14.7E %+14.7E %+14.7E %+14.7E ]\n",
//					x0(0),x0(1),x0(2),x0(3) );
//		printf ( "================================\n" );
//		printf ( "================================\n" );
//		printf ( "================================\n" );
//	}
	xhat=x0;
}

void affinesys::projectonSF0Time(double *t,double xin[]){
	typedef Matrix<double,Dynamic,Dynamic,RowMajor> RowMajorMatrix;
	static VectorXd sfv(affineDim),sfvred(affineDim);
	static VectorXd x0(stateDim+1),dx(stateDim+1),dlambda(affineDim),
								 bx(stateDim+1),blambda(affineDim),x0_ant(stateDim+1);
	static VectorXd lambda0(affineDim),lambda0_ant(affineDim);
	static RowMajorMatrix dsdx(affineDim,stateDim);
	static MatrixXd dsdt(affineDim,1);
	static RowMajorMatrix dsdxred(affineDim+1,stateDim+1);
	static RowMajorMatrix dsdxt(affineDim,stateDim+1);
	static RowMajorMatrix A1(affineDim,affineDim);

	double *xin2=(double*)malloc((stateDim+1)*sizeof(double));

	memcpy((void*)xin2,(void*)xin,(stateDim)*sizeof(double));
	*(xin2+stateDim)=*t;
	Map< VectorXd > xhat(xin2,stateDim+1);
	static size_t subManifoldDim;
	tIndexvector svStatus1;
//	convergetManifolds.print( "----------------------computinf projection on ",
//				   "----------------------------------\n"	);
//std::cout << "--x0 : " << xhat.transpose() << "\n";
	if(convergetManifolds.empty()) return;
	computeSwichingVector(*t,xhat.data(),sfv.data());//,( t > 0.149 && t < 0.153)?1:0);
//std::cout << "--sfv0 : " << sfv.transpose() << "\n";
//	svStatus1.print("actual zero of sfv:");
	convergetManifolds.removeRows(sfvred,sfv);
//	convergetManifolds.print("projecting into :");

//	if( t > 0.149 && t < 0.153){
//			printf ( "[proy] sf = %+14.7E\n",sfv(0) );
//	}
	if(sfvred.norm() < zeroTol) return;

	subManifoldDim=convergetManifolds.size();
	dlambda.resize(subManifoldDim);
	blambda.resize(subManifoldDim);
	dsdxred.resize(subManifoldDim,stateDim);
	A1.resize(subManifoldDim,subManifoldDim);
	lambda0.resize(subManifoldDim);

	x0=xhat;



	lambda0.setOnes();
	lambda0*=0.0005;


	static double sfhis[1000];
	static int counter;
	counter= 0;
	double localGain=localProyGain;

	while(1){
		x0_ant=x0;
		lambda0_ant=lambda0;


		swichingVector_G1(x0(stateDim),x0.data(),dsdx.data());
		swichingVector_partial_dt(x0(stateDim),x0.data(),dsdt.data());
		dsdxt.leftCols(stateDim)=dsdx;
		dsdxt.rightCols(1)=dsdt;
		convergetManifolds.removeRows(dsdxred,dsdxt);

//	if( t > 0.149 && t < 0.153){
//std::cout << "-------x0= [ " << x0.transpose() << "\n";
//std::cout << "--lambda0= [ " << lambda0.transpose() << "]\n ";
//	}
		bx		= - localGain*((xhat-x0)+dsdxred.transpose()*lambda0);
		blambda	= -sfvred-dsdxred*bx;

		A1.noalias()= dsdxred*dsdxred.transpose();

		dlambda=- A1.llt().solve(blambda);
		dx= bx - dsdxred.transpose()*dlambda;

//std::cout << "  x - xhat= [ " << (x0-xhat).transpose() << "]\n ";
		x0+=dx;
		lambda0+= dlambda;

//	if( t > 0.149 && t < 0.153){
////{{{ 
//std::cout << "-----sfv = [ " << sfv.transpose() << "] \n";
//std::cout << "----xhat = [" << xhat.transpose() << "\n";
//std::cout << "----dsdx = [ " << dsdx << "\n-----\n";
//std::cout << "----  bx = [ " << bx.transpose()	<< "]\n";
//std::cout << "dSdx^t dl = [ " << (dsdxred.transpose()*dlambda).transpose()<< "]\n";
//std::cout << "--blambda= [ " << blambda.transpose() << "]\n";
//std::cout << "--localgain =" << localGain<< "\n";
//printf ( "    --      A= %+14.7E\n",A1(0,0) );																					
//std::cout << "\t\t\t\t dsdxred= [ \n" << dsdxred << "\n-----\n";
//std::cout << "\t\t\t\t      dx= [ " << dx.transpose() << "]\n ";
//std::cout << "\t\t\t\t      x0= [ " << x0.transpose() << "\n";
//std::cout << "\t\t\t\t dlambda= [ " << dlambda.transpose() << "]\n ";
//std::cout << "\t\t\t\t lambda0= [ " << lambda0.transpose() << "]\n ";
////}}}
//	}
		computeSwichingVector(x0(stateDim),x0.data(),sfv.data());

		convergetManifolds.removeRows(sfvred,sfv);
//		printf ( "time = %+14.7E\n",t );
//std::cout << "\n---sfv = [ " << sfv.transpose() << "] \n";
//		svStatus1.print("are zero: ",".------------------------\n");
//	convergetManifolds.print("convergetManifolds:");
		
		if(sfvred.norm() < proyTolSfVal && dx.norm() < proyTolDx)
			break;

		sfhis[counter]= sfv(0);	
		if(counter > 1 && 
			(fabs(sfhis[counter]) >= fabs(sfhis[counter-1]))){
			localGain/= 1.2;
			x0=x0_ant;
			lambda0=lambda0_ant;
		}
		counter++;
		if(counter>proyMaxIter){
			printf ( "[warning affine Proy]: "
			"filed on projector at time %+14.7E "
			"wth an error of %+14.7E\n",*t,sfv(0) );
			convergetManifolds.print("[warning affine Proy]: "
							"failed projectin on :");
			std::cout << "[warning affine Proy]: xhat= ["
					<< xhat.transpose() << "\n";
			std::cout << "[warning affine Proy]: x0  = ["
					<< x0.transpose() << "\n";
			printf ( "[warning affine Proy] proyTolSfVal = %+14.7E \n",proyTolSfVal );
			printf ( "[warning affine Proy] proyTolDx = %+14.7E \n",proyTolDx );
				xhat=x0;
			return;
			for(int uici=0;uici < 100; uici++)
					printf ( "sf %i  = %+14.7E\n",uici,sfhis[uici] );
			assert(0);

		}
	};
	
	memcpy((void*)xin,(void*)x0.data(),(stateDim)*sizeof(double));
	*t=x0(stateDim);

	free(xin2);
}
void affinesys::projectonSF0Tangent(double t,const double xin[],double vin[]){
	static Matrix<double,Dynamic,Dynamic,RowMajor> dsdx(affineDim,stateDim);
	static Matrix<double,Dynamic,Dynamic,RowMajor> dsdxred(affineDim,stateDim);
	static size_t subManifoldDim;
	static size_t uici;
	static VectorXd vInNorm(stateDim);
	static double proyscalar;
	Map<VectorXd> v(vin,stateDim);
	if(convergetManifolds.empty()) return;

	subManifoldDim=convergetManifolds.size();
	dsdxred.resize(subManifoldDim,stateDim);

	swichingVector_G1(t,xin,dsdx.data());
	convergetManifolds.removeRows(dsdxred,dsdx);

	vInNorm.setZero();
	for(uici=0;uici<subManifoldDim;uici++){
//		printf ( "dsdxred cols=%zu rows= %zu\n",dsdxred.cols(),dsdxred.rows() );
//		printf ( "v       cols=%zu rows= %zu\n",v.cols(),v.rows());
//		printf ( "dsdxred.row(uici)  cols=%zu rows= %zu\n",
//						dsdxred.row(uici).cols(),dsdxred.row(uici).rows());
		proyscalar= (dsdxred.row(uici)*v);
		proyscalar/= dsdxred.row(uici).norm();
		vInNorm+=proyscalar*dsdx.row(uici).transpose();
	}
	v-=vInNorm;

}

//void affinesys::computeLm(double t,const double xin[],double lm[]){
//	iterator::computeLm(t,xin,lm);
////	if(convergetManifolds.size()!=0)
////	std::cout << "lm proposed  " << Map<RowVectorXd>(lm,stateDim) << "\n";
////	projectonSF0Tangent(t,xin,lm);
////	if(convergetManifolds.size()!=0)
////	std::cout << "lm effective " << Map<RowVectorXd>(lm,stateDim) << "\n";
//}

