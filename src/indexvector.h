
#ifndef  INDEXVECTOR_INC
#define  INDEXVECTOR_INC
#include<vector>
#include<iterator>
#include<algorithm>
#include<eigen3/Eigen/Core>
#include<stdio.h>
#include<complex>
using namespace Eigen;
/** This class is used to represent a set of values in an array
 * by storing a set of indexes. 
 * It also provides functions to compare with other objects of
 * the same class and remove rows or cols from matrix and vectors.  */
class tIndexvector: public std::vector<size_t>{
	public:
	tIndexvector():std::vector<size_t>(0){}
	/** Assigns all elements of the rhs. */
	tIndexvector& operator=(tIndexvector that){
		std::vector<size_t>::operator=(that);
		return *this;
	}
	/** Remove all elements contained in the rhs object. */
	tIndexvector& operator-=(tIndexvector& that){
		tIndexvector::iterator iter1,iter;
		for(iter1=that.begin();iter1!=that.end();iter1++){
			iter=std::find(this->begin(), this->end(), *iter1);
			if(iter!=this->end())
				this->erase(iter);
		}
		return *this;

	}
	/** Test is the index \ref el is stored. In addition it
	 * can return a pointer to it. */
	bool isIn(size_t el,size_t *pos=NULL) {
		tIndexvector::iterator iter;
		iter=std::find(this->begin(), this->end(), el);
		if(pos!=NULL)
			*pos= std::distance( this->begin(),iter);
		return iter!=this->end();
	}
	/** Remove the rows which indexes are contained in the object from matrix \ref src 
	 * and store the result in \ref dest. */
	void removeRows(Matrix<double,Dynamic,Dynamic,RowMajor> &dest,
					Matrix<double,Dynamic,Dynamic,RowMajor> &src) const{
		dest.resize(this->size(),src.cols());
		tIndexvector::iterator iter1;
		static size_t uici;
		for(uici=0;uici<this->size();uici++)
			dest.row(uici)= src.row((*this)[uici]);
	}
	/** Remove the rows which indexes are contained in the object from matrix with c
	 * complex vales\ref src 
	 * and store the result in \ref dest. */
	void removeRows(Matrix<std::complex<double>,Dynamic,Dynamic,RowMajor> &dest,
					Matrix<std::complex<double>,Dynamic,Dynamic,RowMajor> &src) 
		const{
		dest.resize(this->size(),dest.cols());
		tIndexvector::iterator iter1;
		static size_t uici;
		for(uici=0;uici<this->size();uici++)
			dest.row(uici)= src.row((*this)[uici]);
	}
	/** Remove the rows or cols which indexes are contained in the object from vector \ref src 
	 * and store the result in \ref dest. */
	void removeRows(VectorXd &dest,VectorXd &src) const{
		dest.resize(this->size());
		tIndexvector::iterator iter1;
		static size_t uici;
		for(uici=0;uici<this->size();uici++)
			dest(uici)= src((*this)[uici]);
		

	}
	/** Remove the rows or cols which indexes are contained in the object from array \ref src 
	 * and store the result in \ref dest. */
	void removeRows(ArrayXd &dest,ArrayXd &src) const{
		dest.resize(this->size());
		tIndexvector::iterator iter1;
		static size_t uici;
		for(uici=0;uici<this->size();uici++)
			dest(uici)= src((*this)[uici]);
	}
	/** Remove cols which indexes are contained in the object from matrix \ref src 
	 * and store the result in \ref dest. */
	void removeCols(Matrix<double,Dynamic,Dynamic,RowMajor> &dest,
					Matrix<double,Dynamic,Dynamic,RowMajor> &src) const{
		dest.resize(src.rows(),this->size());
		tIndexvector::iterator iter1;
		static size_t uici;
		for(uici=0;uici<this->size();uici++){
			dest.col(uici)= src.col((*this)[uici]);
		}
	}
	/** Remove cols which indexes are contained in the object from complex matrix \ref src 
	 * and store the result in \ref dest. */
	void removeCols(Matrix<std::complex<double>,Dynamic,Dynamic,RowMajor> &dest,
					Matrix<std::complex<double>,Dynamic,Dynamic,RowMajor> &src) 
		const{
		dest.resize(dest.rows(),this->size());
		tIndexvector::iterator iter1;
		static size_t uici;
		for(uici=0;uici<this->size();uici++)
			dest.col(uici)= src.col((*this)[uici]);
	}
	/** Print on stdout a human readable description of the object.*/
	void print(const char *prefix="",const char *postfix="\n"){
		printf ( "%s",prefix );
		printf ( " size: %zu | ",this->size() );
		tIndexvector::iterator iter1;
		for(iter1=this->begin();iter1!=this->end();iter1++)
			printf ( " %zu",*iter1 );
		printf ( "%s",postfix );

	}
};

/** Returns a \ref tIndexvector which contains the elements of \ref left 
 * without the elements of \ref right. Makes a subtraction of sets. */
tIndexvector operator-(tIndexvector left,tIndexvector right);
/** Returns true if all elements of \ref left are contained in \ref right. */
bool operator<(const tIndexvector &left,const tIndexvector &right);
/** Returns true if the \ref left and \ref right contains the same elements. */
bool operator==(tIndexvector &left,tIndexvector &right);
/** Returns true if the \ref left and \ref right contains at least one different element. */
bool operator!=(tIndexvector left,tIndexvector right);


#endif   /* ----- #ifndef INDEXVECTOR_INC  ----- */
