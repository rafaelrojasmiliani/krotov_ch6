#include"csolcont.h"
#include<math.h>

cSolCont::cSolCont(size_t sdim,char ts):
	fileActualPosBytes(0),
	fileActualPosElements(0),
	fileLengthBytes(0),
	fileLengthElements(0),
	file(tmpfile()),
	fileStrideBytes(sizeof(double)*(sdim+1)),
	buffer((double*)malloc(BUFFER_NUMBER*(sdim+1)*sizeof(double))),
	bufferStride(sdim+1),
	dim(sdim),
	sense(ts),
	bufferPos(0),
	readForwardInTimePtr((ts=='+')?		&cSolCont::readForwardInFile:
										&cSolCont::readBackwardInFile),
	readBackwardInTimePtr((ts=='+')?	&cSolCont::readBackwardInFile:
										&cSolCont::readForwardInFile),
	gotoBeginInTimePtr((ts=='+')? 		&cSolCont::fileGotoBegin:
										&cSolCont::fileGotoEnd),
	gotoEndInTimePtr((ts=='+')?			&cSolCont::fileGotoEnd:
										&cSolCont::fileGotoBegin),
	getCurrentTimeInLinesPtr((ts=='+')?	&cSolCont::getFileActualPosInLines:
										&cSolCont::getFileCompltPosInLines),
	isAtEndTimePtr((ts=='+')?			&cSolCont::isAtFileEnd:
										&cSolCont::isAtFileBegin),
	isAtBeginTimePtr((ts=='+')?			&cSolCont::isAtFileBegin:
										&cSolCont::isAtFileEnd),
	pushAndWriteBuff((ts=='+')?			&cSolCont::pushForwardInTime:
										&cSolCont::pushBackwardInTime),
	interpolBuffer((double*)malloc((sdim)*sizeof(double)))
{
	assert(ts=='+' || ts =='-')	;
	size_t uici;
	for(uici=0;uici<BUFFER_NUMBER;uici++){
		buffer_t[uici]=NULL;
		buffer_sol[uici]=NULL;
	}
//	printf ( "file = %p\n",file );

}

cSolCont::cSolCont(const cSolCont& that):
	fileActualPosBytes(0),
	fileActualPosElements(0),
	fileLengthBytes(0),
	fileLengthElements(0),
	file(tmpfile()),
	fileStrideBytes(sizeof(double)*(that.dim+1)),
	buffer((double*)malloc(BUFFER_NUMBER*(that.dim+1)*sizeof(double))),
	bufferStride(that.dim+1),
	dim(that.dim),
	sense(that.sense),
	bufferPos(0),
	readForwardInTimePtr((that.sense=='+')?		&cSolCont::readForwardInFile:
												&cSolCont::readBackwardInFile),
	readBackwardInTimePtr((that.sense=='+')?	&cSolCont::readBackwardInFile:
												&cSolCont::readForwardInFile),
	gotoBeginInTimePtr((that.sense=='+')? 		&cSolCont::fileGotoBegin:
												&cSolCont::fileGotoEnd),
	gotoEndInTimePtr((that.sense=='+')?			&cSolCont::fileGotoEnd:
												&cSolCont::fileGotoBegin),
	getCurrentTimeInLinesPtr((that.sense=='+')?	&cSolCont::getFileActualPosInLines:
												&cSolCont::getFileCompltPosInLines),
	isAtEndTimePtr((that.sense=='+')?			&cSolCont::isAtFileEnd:
												&cSolCont::isAtFileBegin),
	isAtBeginTimePtr((that.sense=='+')?			&cSolCont::isAtFileBegin:
												&cSolCont::isAtFileEnd),
	pushAndWriteBuff((that.sense=='+')?			&cSolCont::pushForwardInTime:
												&cSolCont::pushBackwardInTime),
	interpolBuffer((double*)malloc((that.dim)*sizeof(double)))
{
			
	size_t uici;
	for(uici=0;uici<BUFFER_NUMBER;uici++){
		buffer_t[uici]=NULL;
		buffer_sol[uici]=NULL;
	}
	*this= that;
}


cSolCont* cSolCont::clone(){
	cSolCont* res;
	res= new cSolCont(*this);
	*res = *this;
	return res;
}

cSolCont::~cSolCont(){
	fclose(file);
	free(buffer);
	free(interpolBuffer);
}

const double* cSolCont::interp1
				(double t,double t0,const double x0[],double t1,const double x1[]){
	static size_t uici;
	static double dt;
	dt=t1-t0;
	if(dt==0) dt=1e-4;
	for(uici=0;uici<dim;uici++)
		*(interpolBuffer+uici) =
				(x1[uici]-x0[uici])/dt*(t-t0) + x0[uici];
	return interpolBuffer;

}

//const double* cSolCont::getAtTime_dt(double t){
//	assert(*buffer_t[0] <= t && t <=*buffer_t[BUFFER_NUMBER-1] );
//	static size_t uici;
//
//	for(uici=0;uici<BUFFER_NUMBER;uici++){
//		if(t == *buffer_t[uici])
//	}
//}


