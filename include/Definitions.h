/*
 *  Definitions.h
 *
 *  Created by Gert Thijs on 27/03/06.
 *  Copyright 2006 Silicos NV. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *     + Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     + Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY SILICOS NV AND CONTRIBUTORS ``AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL SILICOS NV AND CONTRIBUTORS BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
		* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
 *
 */

#ifndef __LIB_SIMATH_DEFINITIONS__
#define __LIB_SIMATH_DEFINITIONS__

#include <stdexcept>
#include <cmath>

#define SIMATH_EXCEPTION __attribute__ ((visibility("default")))

//-------------------------------------------------------------------------------------------------------------------------
/**
	\mainpage SiMath : C++ Library for Modeling Data

 	\section _intro_ Introduction
	SiMath is a C++ library of mathematical modeling tools that work on real-valued data matrices. 

	One might wonder why SiMath was conceived since there are already several good libraries available.
  SiMath was originally created within Silicos to have access in stand-alone programs to a collection 
	of available and newly implemented tools used to work with data matrices within the context of chemoinformatics.
	The libray was designed to cover the steps when building a classifier or a predictive model. Such a procedure
	typically starts by preprocessing the data matrix, next some form of feature selection is applied and then a
	model is trained and finally the model should be applied to new data. 
	Mostly, these steps are not all present in one specific library. For instance, libsvm is an 
	excellent library to train a support vector machine but there is no functionality to do a principal component
	analysis. PCA can be done from a SVD which is avalaible in JAMA but there it uses a template-based matrix 
	representation while libsvm works with special sparse matrix representation. Combining both libraries directly
	into one program is as such a non-trivial task.
	
	The goal of SiMath was to create a simple and consistent interface to existing tools which could be easily
	added into C++ applications. So, each of the included models has more or less the same set of methods. 
	For instance, each model is initiated with a set of parameters and the model has a method to train the 
	model directly from a data matrix or vector. There are also 
	
  Since SiMath is only intended to work with real-valued matrices and vectors. Therefore, the Matrix and Vector
	class only support double values. All included tools are adapted in such a fashion that they can work directly 
	with Matrix and Vector objects. By not using a template-based implementation the design and implementation of 
	the interface of library was kept simple and rather straightforward. 
	
	Where possible, it was also decided to work with classes and algorithms from the std namespace. A typical example
	is the Vector classes which holds a std::vector<double> as internal data structure and not an array of double. 
	As such, operations on arrays can be efficiently passed to the standard algorithms 
	which decreases the probability of errors and memory leaks.
 
	There is no class available to do the IO of matrices and vectors. This is done to give the user of the
	library as much freedom as possible in the way the data are represented. Again, libsvm, for instance, imposes 
	a strict representation of a sparse matrix while TNT uses a completely different scheme to write a matrix. 
 
	\section _usage_ Usage
	
	Below is an example of a (incomplete) code snippet to do a PCA analysis on a calibrated input matrix and 
	then cluster the data with k-means based on the truncated data matrix. 
	\code 
	// example.cpp
	#include "SiMath/SiMath.h"
 
	SiMath::Matrix A(n,m,0.0); // initialise an [nxm] matrix with zeros
	// add some code to fill the matrix
	for ( int i=0; i<n; ++i )
		for ( int j=0; j<m; ++j )
			A[i][j] = some_input_function(i,j);
	
	// calibrate matrix by setting each column to have a zero mean and unit variance
	SiMath::Vector m = SiMath::meanOfAllColumns(A);
	SiMath::Vector s = SiMath::stDevOfAllColumns(A,m);
	columnNormalise(A,m,s);
 
	// to do a PCA analysis, do a SVD and store the V matrix
	SiMath::SVD svd(A,false,true); 
	
	// get the V matrix
	SiMath::Matrix V = svd.getV();

	// get the singular values and normalise them to get variances
	SiMath::Vector sv = svd.getSingularValues();
	sv /= sqrt(n-1);
	
	// SV's are sorted, find those that are above a given threshold
  int c = 0;
	for ( ; c<sv.size(); ++c ) {
		if ( sv[c] < 0.1 )
			break;
	}
	// c now holds index of first SV lower than threshold
  // if all are above treshold c is equal to sv.size()
  
	//approximate V with c best singular values
	SiMath::Matrix Vc(m,c,0);
	for ( int i=0; i<m; ++i )
		for ( int j=0; j<c; ++j )
			Vc[i][j] = V[i][j];
 
	// do the rotation of A with app
	SiMath::Matrix Arot = product(A,Vc);
	
	// cluster -n 10 clusters the data with k-means clustering
	SiMath::KMeansParameters p;
	p.dimension = c;     // dimension of the data is c
	p.nbrClusters = 10;  // nbr of clusters

	// initialise the model
	SiMath::KMeans model(p);

	// cluster the data and return the cluster labels for each data point
	std::vector<unsigned int> labels = model.cluster(Arot);
			
	// print out the labels
	for ( int i=0; i<n; ++i )
		std::cout << i << "\t" << labels[i] << std::endl;
 
	\endcode

	\section _examples_ Examples
	In the directory <b>examples</b> of the source distribution, there are currently two example programs and 
	one data file <b>data.tab</b>. The data matrix comes from the benchmark iris data set (used for instance
	in R). It has 150 data points with 5 variables. The first column gives the class label. There are three
	different classes and each class has 50 members.
 
	The first example reads in the data, normalises the columns to have zero mean and unit variance. 
	This normalised data matrix is then used to train an SVM classificator. The predicted labels are
	printed.
 
	The second example reads in the data and the number of clusters requested.
	Then it normalises the columns to have zero mean and unit variance  and transforms 
	the data matrix using PCA. Finally a k-means clustering is done using the predefined number 
	of clusters. The cluster numbers are printed in 3 columns corresponding to the three classes.
 
	To compile the code against the compiled version of libSiMath, you should use a command like 
	(if SiMath is installed in the default location)
	\code
		# first example
		$ g++ -I/usr/local/include -L/usr/local/lib -lSiMath -o example1 example1.cpp
		$ ./example1 data.tab
		
		# second example 
		$ g++ -I/usr/local/include -L/usr/local/lib -lSiMath -o example2 example2.cpp
		$ ./example2 data.tab 3
	\endcode 
	
 
	\section _installation_ Installation
	SiMath uses the GNU autoconf, so you can build it with the standard commands:
	\code
		$ ./configure
		$ make
		$ sudo make install
	\endcode
	To install SiMath in another directory than the default /usr/local, you can use the '--prefix'
	option of the configure script
	\code 
		$ ./configure --prefix=/my/install/dir
	\endcode
	
	To use SiMath in your application add the following include statement to your code
	\code
	#include "SiMath/SiMath.h"
	\endcode
 
	\section _links_ Links
	Here is a list of links to resource that have been used in creating SiMath: 
		- libsvm : http://www.csie.ntu.edu.tw/~cjlin/libsvm/
		- tnt/jama : http://math.nist.gov/tnt/index.html
		- gallery of distribution: http://www.itl.nist.gov/div898/handbook/eda/section3/eda366.htm
     
	\section _copy_ COPYRIGHT
	Copyright 2006 Silicos NV. All rights reserved.
 
  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:
		- Redistributions of source code must retain the above copyright
			notice, this list of conditions and the following disclaimer.
		-	Redistributions in binary form must reproduce the above copyright
			notice, this list of conditions and the following disclaimer in the
			documentation and/or other materials provided with the distribution.
 
  THIS SOFTWARE IS PROVIDED BY SILICOS NV AND CONTRIBUTORS ``AS IS'' AND ANY
  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  DISCLAIMED. IN NO EVENT SHALL SILICOS NV AND CONTRIBUTORS BE LIABLE FOR ANY
	DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
	(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
	LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 		
 */


//-------------------------------------------------------------------------------------------------------------


/**
	\namespace SiMath
	\brief Main namespace to hold all classses and functions part of the SiMath toolkit.
 
	\todo Further update and refine the API documentation
	\todo Include logistic regression models
	\todo Implement new performance criteria for the different modeling tools 
	\todo Add other model classes (eg. bayesian clustering, ...)
	
 */
namespace SiMath {

	// Define groups to appear as such in the documentation
	/** 
		\defgroup _error_handling_ Error Handling Classes 
		\brief All Error handling classes
		
		All error handling classes are derived from std::runtime_error. 
		To make them accessible outside the library, a visibility attribute is added in the class definition. 
	 */
	/** 
		\defgroup _model_parameters_ Model Parameters
		\brief Overview of Model parameter structure
	
		Each model present in SiMath should be initialised with a set of parameters.
	 */
	/** 
		\defgroup _clustering_ Clustering Algorithms 
		\brief Overview of Clustering Algorithms
	 */
	/** 
		\defgroup _classifiers_ Classifier Algorithms 
		\brief Overview of Classifier Algorithms
	 */
	
	
/**
	\class  SiMathError
	\brief  Exception class
 
	This exception class is derived from the std::runtime_error
	class and can be used by the different classes to 
	indicate that something went wrong during the calculations.
 
	\param m error message
	\ingroup _error_handling_
 */
	class SIMATH_EXCEPTION SiMathError: public std::runtime_error
	{
		public:
			/** Constructor
				\param     message Additional information about the error that occured.
			 */
			SiMathError(const std::string & message): runtime_error("ERROR: SiMathError : " + message) {};
	};

#define INF HUGE_VAL
#define PI 3.14159265358979323846
#define TAU 1e-12
	
#ifndef min
	template <class T> inline T min(T x,T y) { return (x<y)?x:y; }
#endif

#ifndef max
	template <class T> inline T max(T x,T y) { return (x>y)?x:y; }
#endif

#ifndef sign
	template <class T> inline T sign(const T & a, const T & b ) {return (b >= 0.0) ? ( a>=0 ? a : -a) : (a>=0 ? -a : a);}
#endif sign
	
//#ifndef swap
//	template <class T> inline void swap(T& x, T& y) { T t=x; x=y; y=t; }
//#endif swap

	
};


#endif __LIB_SIMATH_DEFINITIONS__

