/*!
 *       \file:  indexvector.cpp
 *
 *		 \brief: DESCRIPTION
 *        Version:  1.0
 *        Created:  14/03/15 15:13:15
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include"indexvector.h"

tIndexvector operator-(tIndexvector lhs,tIndexvector rhs){
	static int uici;
	static size_t sizeL,sizeR;
	tIndexvector::iterator iter,iter1,iter2;
	tIndexvector result;
	for(iter1=lhs.begin();iter1!=lhs.end();iter1++){
		iter=std::find(rhs.begin(), rhs.end(),*iter1);
		if(iter==rhs.end())
			result.push_back(*iter1);
	}
	for(iter1=rhs.begin();iter1!=rhs.end();iter1++){
		iter=std::find(lhs.begin(), lhs.end(),*iter1);
		if(iter==lhs.end()){
			iter=std::find(result.begin(), result.end(),*iter1);
			if(iter==result.end())
				result.push_back(*iter1);
		}
	}
	std::sort(result.begin(),result.end());

}
bool operator==(tIndexvector &lhs,tIndexvector &rhs){
	if(lhs.size()!=rhs.size()) return 0;
	static size_t uici;
	for(uici=0;uici<lhs.size();uici++)
		if(lhs[uici]!=rhs[uici])
			return 0;
	return 1;
}
bool operator!=(tIndexvector lhs,tIndexvector rhs){
	return !(lhs==rhs);
}
bool operator< (const tIndexvector &lhs, const tIndexvector& rhs){
	tIndexvector::const_iterator iter1;
	tIndexvector::const_iterator iter;

	for(iter1=lhs.begin();iter1!=lhs.end();iter1++){
		iter=std::find(rhs.begin(), rhs.end(),*iter1);
		if(iter==rhs.end()){
			return 0;
		}
	}

	return 1;
}
