
#ifndef  CSOLCONT_INC
#define  CSOLCONT_INC
#include<stdlib.h>
#include<stdio.h>
#include<assert.h>
#include<utility>

#define BUFFER_NUMBER  5
typedef std::pair<const double*,const double*> tPtrPair;

/** A Solution Container. This class is intended to represent the time history of a vector-valued 
 * function which is stored as a matrix in a file.
 * The time story of this vector-valued function is stored in a binary file composed
 * of blocs of \ref dim \f$+1\f$  doubles. The first element of this blocs represents the time,
 * the other values are the vector-valued function's components.  
 * This class furnish methods to access to the value of the stored vector-valued function
 * to each instant in time through an interpolation method.
 * This class may be created to store a vector-valued function which is created adding values in
 * time-increasing order or in time-decreasing order.
 * */
class cSolCont{
	private:
		long fileActualPosBytes;//!< The    actual position in the file in bytes.
		long fileActualPosElements;//!< The actual position in the file in elements, blocks of doubles.
		long fileLengthBytes; //!< The length of the file in bytes.
		long fileLengthElements;//!< The length of the file in elements, blocks of doubles.
		FILE* const file;//!< The file where the information is stored.
		const long fileStrideBytes; //!< The distance in bytes between each block which represents an instant.
		double* const buffer;//!< Memory buffer to be used by the interpolation process.
		const size_t bufferStride;//!< the distance in bytes between each block which represents an instant in the buffer.
		const size_t dim;//!< the dimension of the vector valued function.
		const char sense;//!< the sense in which the vector valued function is stored, positive or negative in time.
		size_t bufferPos;//!< actual position in the buffer.
		const double* readForwardInFile();//!<  reads a block in the file, in forward direction.
		const double* readBackwardInFile();//!< reads a block in the file, in backwards direction.
		void fileGotoBegin();//!< place the \ref file stream at the beginning of the file.
		void fileGotoEnd();//!< place   the \ref file stream at the end of the file.
	protected:
		bool isAtFileBegin()/*< returns true if the stream \ref file is at the beginning.*/ const{
			return (fileActualPosElements==0);
		}
		bool isAtFileEnd()/*< returns true if the stream \ref file is at the end.*/ const{
			return (fileActualPosElements==fileLengthElements);
		}
		/** Write \ref dim \f$+1\f$ doubles in the file in increasing time.*/
		void pushForwardInTime(double t,const double* x) ;
		/** Write \ref dim \f$+1\f$ doubles in the file in decreasing time.*/
		void pushBackwardInTime(double t,const double* x) ;

		const double* (cSolCont::*const readForwardInTimePtr)();//!< A pointer to the function which reads in increasing time.
		const double* (cSolCont::*const readBackwardInTimePtr)();//!< A pointer to the function which reads in decreasing time.
		void (cSolCont::*const gotoBeginInTimePtr)();//!< A pointer to the func. which places the stream \ref file where time begins.
		void (cSolCont::*const gotoEndInTimePtr)();//!<   A pointer to the func. which places the stream \ref file where time ends.
		long (cSolCont::*const getCurrentTimeInLinesPtr)() const;//!< Get the time in number of stored values.
		bool (cSolCont::*const isAtEndTimePtr)()const; //!< A pointer the func. which says is the stream\ref file is where time begins.
		bool (cSolCont::*const isAtBeginTimePtr)()const;//!< A pointer the func. which says is the stream\ref file is where time ends.
		void (cSolCont::*const pushAndWriteBuff)(double,const double*);//!< A pointer to the function which reads on the file.
		double* const interpolBuffer; //!< The buffer for the interpolation computation process.
	public:
		const double* buffer_t[BUFFER_NUMBER];//!<   An array elements of \ref buffer with the time.
		const double* buffer_sol[BUFFER_NUMBER];//!< An array elements of \ref buffer without the time.
		const double* readForwardInTime(); //!< Reads the stream in direction where time increases. Stores the values in the buffer.
		const double* readBackwardInTime();//!< Reads the stream in direction where time decreases. Stores the values in the buffer.
		tPtrPair operator()(double );      //!< Returns a pair of pointers to two consecutive blocks in the buffer.
		/** Returns the actual position of the stream \ref file in number of blocks of \ref dim \f$+1\f$ doubles. */
		long getFileActualPosInLines() const{
			return fileActualPosElements;
		}
		/** Returns number of blocks of \ref dim \f$+1\f$ doubles until the end of the file from the actual position of the stream
		 * \ref file. */
		long getFileCompltPosInLines() const{
			return fileLengthElements-fileActualPosElements;
		}
		/** Returns the number of blocks of \ref dim \f$+1\f$ doubles in the file.*/
		long size() const{
			return fileLengthElements;
		}
		/** Returns the actual position in the stream \ref file in number of blocks of \ref dim \f$+1\f$ doubles.*/
		long getActualTimeInLines() const{
			return (this->*cSolCont::getCurrentTimeInLinesPtr)();
		}
		/** Returns true is the stream \ref file is located in the block that represent the last value in time. */
		bool isAtEndTime() const{
			return (this->*cSolCont::isAtEndTimePtr)();		
		}
		/** Returns true is the stream \ref file is located in the block that represent the first value in time. */
		bool isAtBeginTime() const{
			return (this->*cSolCont::isAtBeginTimePtr)();
		}
		/** Returns true if time \ref t is located inside the values used for the interpolation process. */
		bool timeIsInBuffer(double t) const{
			return (*buffer_t[0] <= t &&
				t <=*buffer_t[BUFFER_NUMBER-1]);
		}
		/** Same as \ref size. */
		size_t sizeElements() const {return fileLengthElements;}
		/** Same as \ref size but in bytes.*/
		size_t sizeBytes() const {return fileLengthBytes;}
		/** Place the stream \ref file in the position that represents the first value in time of the stored vector-valued function. */
		void gotoBeginInTime();
		/** Place the stream \ref file in the position that represents the last value in time of the stored vector-valued function. */
		void gotoEndInTime();

		cSolCont(size_t sdim,char ts='+');
		cSolCont(const cSolCont& that);
		cSolCont& operator=(const cSolCont& that);
		virtual ~cSolCont();

		cSolCont* clone();

		/** Add a new value to the time history of the vector valued function. */
		void push(double t,const double* x){
			(this->*cSolCont::pushAndWriteBuff)(t,x);
		}
		/** A debug function. */
		void printActualInterval(const char* prefix=""){
				printf ( "%s   actual interval: %+14.7E    %+14.7E\n",
					  prefix,*buffer_t[0],*buffer_t[BUFFER_NUMBER-1] );

		}
		/** Writes \ref dim \f$\times\f$ sizeof(double) bytes of \ref date in the stream \ref file. */
		const double* pushInFile(double t,const double* data) ;
		/** Makes a linear interpolation. */
		const double* interp1
			(double t,double t0,const double x0[],double t1,const double x1[]);


		void update(double t);
		/** Writes the data of the file \ref file into an ascii file with a given format \ref fmt. */
		void writeascii(const char* name,const char* fmt="%+14.7E \0");
		/** Clears all the memory. */
		void clear();
		/** Computes the approximate value of the stored vector valued function at time \ref t. */
		const double* getAtTimeInterp1(double t);

		/** Get two pointer to the buffer blocks used for an interpolation at time \ref t. */
		void getAtTimeDataToIterp(double t,const double** b0, const double** b1);

		/** Load a file. */
		void load(const char*);
		/** Save to a file. */
		void save(const char*);
};
#endif   /* ----- #ifndef CSOLCONT_INC  ----- */
