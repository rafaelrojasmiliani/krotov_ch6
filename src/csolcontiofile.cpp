/*!
 *       \file:  csolcontiofile.cpp
 *
 *		 \brief: Here we hace direct interaction with the file
 *        Version:  1.0
 *        Created:  03/06/15 11:25:05
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */
#include"csolcont.h"
#include<stdlib.h>
#include<string.h>

/** Here we write in the file, but first we add this to the buffer,
 * Then we actualize the possition indications of the buffer,
 * and the file, the first assertion is to assure that
 * we are only writing in a file that as been just created
 * or have been cleaded. */
const double* cSolCont::pushInFile(double t, const double *x){
	assert(fileActualPosBytes==fileLengthBytes);

	if(bufferPos == BUFFER_NUMBER) bufferPos=0;
	double *cbp= buffer+bufferPos*bufferStride;
	bufferPos++;

	memcpy((void*)cbp,(const void*)&t,1*sizeof(double));
	memcpy((void*)(cbp+1),(const void*)x,dim*sizeof(double));

	fwrite((const void*)cbp,sizeof(double),bufferStride,file);

	fileActualPosBytes= ftell(file);
	fileLengthBytes= fileActualPosBytes;
	fileActualPosElements=	fileActualPosBytes/(sizeof(double)*bufferStride);
	fileLengthElements= fileActualPosElements;
	return cbp;
}


/** Here we read in the forward sense of the file, then 
 * we actualize the actual position variables. The data is
 * placed on the buffer */
const double* cSolCont::readForwardInFile(){
	if(bufferPos == BUFFER_NUMBER) bufferPos=0;
	double *cbp= buffer+bufferPos*bufferStride;
	bufferPos++;

	fread((void*)cbp,sizeof(double),bufferStride,file);

	fileActualPosBytes=		ftell(file);
	fileActualPosElements=	fileActualPosBytes/(sizeof(double)*bufferStride);

	return cbp;
}

/** Here we read in the backward sense of the file using
 * the fseek function, then 
 * we actualize the actual position variables. The data is
 * placed on the buffer */
const double* cSolCont::readBackwardInFile(){
	if(bufferPos == BUFFER_NUMBER) bufferPos=0;
	double *cbp= buffer+bufferPos*bufferStride;
	bufferPos++;

	fseek(file,-fileStrideBytes,SEEK_CUR);
	fread((void*)cbp,sizeof(double),bufferStride,file);
	fseek(file,-fileStrideBytes,SEEK_CUR);

	fileActualPosBytes=		ftell(file);
	fileActualPosElements=	fileActualPosBytes/(sizeof(double)*bufferStride);

	return cbp;
}

/** Here we place the file at the end of the last written line, because
 * pushInFile is the only function that modify fileLengthBytes, and
 * actualize the current position variables.*/
void cSolCont::fileGotoEnd(){
	fseek(file,fileLengthBytes,SEEK_SET);
	fileActualPosBytes=fileLengthBytes;
	fileActualPosElements=fileActualPosBytes/(sizeof(double)*bufferStride);
}

/** Here we place the file at its begin */
void cSolCont::fileGotoBegin(){
	fseek(file,0L,SEEK_SET);
	fileActualPosBytes=0;
	fileActualPosElements=0;
	bufferPos=0;
}

/** Where we place the file at its begin and set to zero
 * the variables of the lenght of the file and its actual position */
void cSolCont::clear(){
	fseek(file,0L,SEEK_SET);
	fileActualPosBytes=0;
	fileActualPosElements=0;
	fileLengthBytes=0;
	fileLengthElements=0;
	bufferPos=0;
	size_t uici;
	for(uici=0;uici<BUFFER_NUMBER;uici++){
		buffer_t[uici]=NULL;
		buffer_sol[uici]=NULL;
	}
}

/** The asignation procedure consinst in copying the
 * file of the rhs variable.  This is done dong this:
 * 1) The object is cleaded, i.e., we delete its
 * 		actual data.
 * 2) we rewind the rhs file
 * 3) we start to read the rhs and copy its information
 * 		intp te lhs file until the end of file
 *4) then the replace the rhs file into its original position.
 * 		*/

cSolCont& cSolCont::operator =(const cSolCont& that){
	clear();
	double *ab;
	long thatActualPos;
	size_t rs;

	assert(this->dim == that.dim);
	assert(this->sense == that.sense);


	
	ab= (double*)malloc(bufferStride*sizeof(double));
	
	thatActualPos= ftell(that.file);
	fseek(that.file,0,SEEK_SET);

	while((rs = fread(ab,sizeof(double),bufferStride,that.file))){
		push(*ab,(ab+1));
	}
	fseek(that.file,thatActualPos,SEEK_SET);


	free (ab);

	fileLengthBytes		= that.fileLengthBytes;
	fileLengthElements 	= that.fileLengthElements; 
	gotoBeginInTime();

	return *this;
}


void cSolCont::writeascii(const char* name,const char* fmt){
	FILE* output=NULL;
	const double *cbp;	
	size_t uici;
	output=fopen(name,"w");

	if(fileLengthBytes==0) return;
	(this->*cSolCont::gotoBeginInTimePtr)();
	while(!(this->*cSolCont::isAtEndTimePtr)()){
		cbp=(this->*cSolCont::readForwardInTimePtr)();
		for(uici=0;uici<dim+1;uici++)
			fprintf(output,fmt,*(cbp+uici));
		fprintf(output,"\n");
	}
	fclose(output);
	output=NULL;
}

void cSolCont::load(const char* name){

	double *ab= (double*)malloc(bufferStride*sizeof(double));
	size_t rs;
	FILE* in=fopen(name,"rb");
	clear();

	while((rs = fread(ab,sizeof(double),bufferStride,in)))
		push(*ab,(ab+1));

	free (ab);
	fclose(in);
	gotoBeginInTime();
}

void cSolCont::save(const char* name){

	double *ab= (double*)malloc(bufferStride*sizeof(double));
	size_t rs;
	FILE* out=fopen(name,"wb");

	fseek(file,0,SEEK_SET);
	while((rs = fread(ab,sizeof(double),bufferStride,file)))
		fwrite((void*)ab,sizeof(double),bufferStride,out);

	fclose(out);
	free (ab);
	gotoBeginInTime();

}
/** We place the file at the position of the begin
 * of the place element. and actualize the position variables.
 * The we have to fill the buffer to contain that position.  */
//void cSolCont::operator[](const long place){
//	assert(place >= 0 && place <= fileLengthElements);
//	long placeInBytes= place*bufferStride*sizeof(double);
//	int aux= fileLengthElements-place;
//	static const double* cbp;
//	if(place < BUFFER_NUMBER){
//		goto
//	}
//	if(fileLengthElements-place < BUFFER_NUMBER){
//		/** retyrb the place in buffer_t, bu the rigth */
//	}
//
//
//	fseek(file,placeInBytes-BUFFER_NUMBER/2,SEEK_SET);
//	for(uici=0;uici<BUFFER_NUMBER;uici++){
//		cbp=readForwardInFile();
//		if(sense=='+'){
//			buffer_t[uici]=cbp;
//			buffer_sol[uici]=cbp+1;
//		}else{
//			buffer_t[BUFFER_NUMBER-1-uici]=cbp;
//			buffer_sol[BUFFER_NUMBER-1-uici]=cbp+1;
//		}
//
//	}
//
//	
//}

//tPtrPair FinalTimePair(){
//	tPtrPair res;
//	gotoEndInTime();
//	res.second= buffer_[BUFFER_NUMBER-1];
//	res.first= buffer_[BUFFER_NUMBER-2];
//}
//
//
//
//tPtrPair beginTimePair(){
//	
//}
//
//tPtrPair cSolCont::operator++(){
//
//}
