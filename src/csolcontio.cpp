
#include"csolcont.h"
#include<string.h>
#include<math.h>
#include<stdlib.h>

/** Write to the file and actualize the buffer direct access
 * variables buffer_t and buffer_sol, The last write is
 * place and the end of the buffer direct access varaibles, */
void cSolCont::pushForwardInTime(double t,const double* x){
	const double* cbp;
	size_t uici;
	cbp= pushInFile(t,x);

	for(uici=0;uici <BUFFER_NUMBER-1;uici++){
		buffer_t[uici]=buffer_t[uici+1];
		buffer_sol[uici]=buffer_sol[uici+1];
	}
	buffer_t[uici]=cbp;
	buffer_sol[uici]=cbp+1;

}

/** Write to the file and actualize the buffer direct access
 * variables buffer_t and buffer_sol, The last write is
 * place and the begin of the buffer direct access varaibles, */
void cSolCont::pushBackwardInTime(double t,const double* x){
	const double* cbp;
	size_t uici;
	cbp= pushInFile(t,x);

	for(uici=BUFFER_NUMBER-1;uici > 0;uici--){
		buffer_t[uici]=buffer_t[uici-1];
		buffer_sol[uici]=buffer_sol[uici-1];
	}
	buffer_t[0]=cbp;
	buffer_sol[0]=cbp+1;
}

void cSolCont::gotoBeginInTime(){
	(this->*cSolCont::gotoBeginInTimePtr)();
	for(size_t uici=0;uici<BUFFER_NUMBER;uici++)
		readForwardInTime();
	(this->*cSolCont::gotoBeginInTimePtr)();
}

const double* cSolCont::readForwardInTime(){
	static const double* cbp;
	static size_t uici;
	cbp=(this->*cSolCont::readForwardInTimePtr)();
	long place= getActualTimeInLines();
	if(place >= BUFFER_NUMBER){
		for(uici=0;uici <BUFFER_NUMBER-1;uici++){
			buffer_t[uici]=buffer_t[uici+1];
			buffer_sol[uici]=buffer_sol[uici+1];
		}
		buffer_t[uici]=cbp;
		buffer_sol[uici]=cbp+1;
	}else{
		buffer_t[place]=cbp;
		buffer_sol[place]=cbp+1;
	}
	return cbp;
}
void cSolCont::gotoEndInTime(){
	(this->*cSolCont::gotoEndInTimePtr)();
	for(size_t uici=0;uici<BUFFER_NUMBER;uici++)
		readBackwardInTime();
	(this->*cSolCont::gotoEndInTimePtr)();
}
const double* cSolCont::readBackwardInTime(){
	static const double* cbp;
	static size_t uici;
	long place= getActualTimeInLines();
	long cmplplace=size()-place;
	cbp=(this->*cSolCont::readBackwardInTimePtr)();
	if(cmplplace >=BUFFER_NUMBER){
//			printf ( "entramos\n" );
		for(uici=BUFFER_NUMBER-1;uici > 0;uici--){
//				printf ( "uici= %zu\n",uici );
			buffer_t[uici]=buffer_t[uici-1];
			buffer_sol[uici]=buffer_sol[uici-1];
		}
//			printf ( "salimos\n" );
		buffer_t[0]=cbp;
		buffer_sol[0]=cbp+1;
	}else{
		buffer_t[BUFFER_NUMBER-1-cmplplace]=cbp;
		buffer_sol[BUFFER_NUMBER-1-cmplplace]=cbp+1;
	}
	return cbp;
}


tPtrPair cSolCont::operator()(double t){
	tPtrPair res(NULL,NULL);
	if(t < 0.0){
		gotoBeginInTime();
		res.first= buffer_t[0];
		res.second= buffer_t[1];
		return res;
	}else
	if(t < *buffer_t[0]){
		res.first =buffer_t[0];
		res.second=buffer_t[1];
		if(isAtBeginTime()) return res;
		do{
			res.second=res.first;
			res.first=readBackwardInTime();
		}while(t < *res.first && !isAtBeginTime());
		
		return res;
		
	}else
	if(*buffer_t[0] <= t && t <= *buffer_t[BUFFER_NUMBER-1] ){
		size_t uici;
		for(uici=0;uici < BUFFER_NUMBER;uici++){
			if(t == *buffer_t[uici]){
				if(uici != BUFFER_NUMBER-1){
					res.first=buffer_t[uici];
					res.second=buffer_t[uici+1];
				}else{
					res.first=buffer_t[uici-1];
					res.second=buffer_t[uici];
				}
				return res;
			}
			if(t<*buffer_t[uici]){
				assert( uici != 0);
				res.first=buffer_t[uici-1];
				res.second=buffer_t[uici];
				return res;
			}
		}
	}else
	if(*buffer_t[BUFFER_NUMBER-1] < t){
		res.first=buffer_t[BUFFER_NUMBER-2];
		res.second=buffer_t[BUFFER_NUMBER-1];
		if(isAtEndTime()) return res;
		do{
			res.first=res.second;
			res.second=readForwardInTime();
		}while(*res.second < t && !isAtEndTime());
		return res;
		
	}

	assert(0);

}
	
/** We should use this function after the buffer is full.
 * This function returns a local memory buffer where
 * there is the firsth degree aproximation to the
 * solution at time t. The procedure is the following:
 * 1) if t < 0 we return an interpolation
 * 		using the first and second element of te file
 * 	*/
const double* cSolCont::getAtTimeInterp1(double t){
	tPtrPair ap;
	ap= (*this)(t);

	return interp1(t,*ap.first,(ap.first+1),*ap.second,(ap.second+1));
}

void cSolCont::getAtTimeDataToIterp(double t,const double** b0, const double** b1){
	tPtrPair ap;
	ap= (*this)(t);
	*b0= ap.first;
	*b1= ap.second;
	
}

//tPtrPair stepBackInTime(){
//	long place=getTimePosElement();	
//	const double *cbp;
//	cbp= readBackwardInTime();
//	
//
//}
