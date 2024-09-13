/*
 *  Kernel.h
 *
 *  Created by Gert Thijs on 04/04/06.
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

#ifndef _SIMATH_SVM_KERNEL__
#define _SIMATH_SVM_KERNEL__

#include <cmath>
#include <string>
#include <vector>
#include "Matrix.h"
#include "Vector.h"

#include "Cache.h"
#include "SVM.h"

namespace SiMath
{
	/**
		\class Kernel 
	  \brief Base class to represent a kernel matrix
	 
    This is a pure virtual class which should not be used directly.
    Use one of the derived classes which contain kernel matrix representations
    especially tuned for specific problems formulations like classification and regression.   
   */
	class Kernel
	{
		public:
			/**
				\brief Constructor of the kernel based on the data and the parameters
				\param prob Reference to the trainings data of the SVM
				\param param Reference to the set of SVM Parameters
			 */
			Kernel(SVMProblem & prob, const SVMParameters & param);
	
	/**
	\brief Destructor
	 */
	virtual ~Kernel();
	
	/**
		\brief Pure virtual function to get a pointer to the i-th column of the kernel matrix
		
	 */
	virtual double * getQ(unsigned int i, unsigned int len) const = 0;
	
	/**
		\brief Pure virtual function to get the diagonal of the kernel matrix
	 */
	virtual double * getQD() const = 0;
	
	/**
		\brief Pure virtual function to swap columns i and j in the kernel matrix
	 */
	virtual void swapIndex(unsigned int i, unsigned int j) = 0;
	virtual	void printIndex() = 0;
	
	/**
		\brief definition of method for single kernel evaluation
	 */
	static double function(unsigned int l, const double * x, const double * y, const SVMParameters & param);
	
	protected:
		double (Kernel::*_kernelFunction)(unsigned int i, unsigned int j) const;
		SVMProblem * _dataPointer;               ///< holds pointer to original data
		double *  _xSquare;                      ///< stores x[i]^2, if needed 
	
	private:			
    //@{
		
		const KernelType _kernelType;  ///< Local stroage of kerneltype
		const unsigned int _degree;    ///< Degree of the polynomial kernel
		const double _gamma;           ///< Gamma parameter
		const double _coef0;           ///< Offset coefficient
	
	/**
		\brief Local helper function to calculate dot-product
	 \param n Length of the vector
	 \param px Pointer to x data array
	 \param py Pointer to y data array
	 */
	static double _dot(unsigned int n, const double *px, const double *py);
	
	/**
		\brief Implementation of linear kernel
	 \param i Index of first data point
	 \param j Index of second data point
			*/
	double _linearKernel(unsigned int i, unsigned  int j) const;
	
	/**
		\brief Implementation of polynomial kernel
	 \param i Index of first data point
	 \param j Index of second data point
			*/
	double _polyKernel(unsigned int i, unsigned int j) const;
	
	/**
		\brief Implementation of RBF kernel
	 \param i Index of first data point
	 \param j Index of second data point
			*/
	double _rbfKernel(unsigned int i, unsigned int j) const;
	
	/**
		\brief implementation of sigmoidal kernel
	 \param i Index of first data point
	 \param j Index of second data point
			*/
	double _sigmoidKernel(unsigned int i, unsigned int j) const;
	
	/**
		\brief implementation of spline kernel
	 \param i Index of first data point
	 \param j Index of second data point
			*/
	double _splineKernel(unsigned int i, unsigned int j) const;
	//@}
};


	/**
		\brief Kernel matrix for support vector classification problem
 
		Kernel data is internally stored in a cache, such that not the
		full matrix should be stored.
	 */
	class svcQMatrix: public Kernel
	{ 
	public:
		/**
			\brief Constructor of the kernel based on the data and the parameters
			\param prob Reference to the training data of the SVM
			\param param Reference to the set of SVM parameters
		 */
		svcQMatrix(SVMProblem & prob, const SVMParameters & param);
	
		/**
			\brief Destructor
		 */
		~svcQMatrix();
	
		/**
			\brief Get the 
		 \param i index of the row
		 \param len 
		 */
		double * getQ(unsigned int i, unsigned int len) const;
		
		/**
			\brief Get the diagonal of the kernel matrix
		 */
		double * getQD() const { return _QD; };
		
		/** 
			\brief Method to swap two columns in the kernel matrix
			\param i Index of the first column
			\param j Index of the second column
			*/	
		void swapIndex(unsigned int i, unsigned int j);
		void printIndex();

	private:
		std::vector<unsigned int> _index;   ///< Local storage of the column indices
		Cache * _cache;                     ///< Pointer to the cache where the actual data is stored
		double * _QD;                       ///< Local array to hold the diagonal
};


class oneClassQMatrix: public Kernel
{
public:
	/**
		\brief Constructor of the kernel based on the data and the parameters
		\param prob Reference to the training data of the SVM
		\param param Reference to the set of SVM parameters
	 */
	oneClassQMatrix(SVMProblem & prob, const SVMParameters & param);

	/**
		\brief Destructor
	 */
	~oneClassQMatrix();
	
	/**
		\brief Get the 
	 \param i index of the row
	 \param len 
	 */
	double * getQ(unsigned int i, unsigned int len) const;
	
	/**
		\brief Get the diagonal of the kernel matrix
	 */
	double * getQD() const { return _QD; };
	
	/** 
		\brief Method to swap two columns in the kernel matrix
		\param i Index of the first column
		\param j Index of the second column
		*/	
	void swapIndex(unsigned int i, unsigned int j);
	void printIndex();
private:	
		std::vector<unsigned int> _index;   ///< Local storage of the column indices
		Cache * _cache;                     ///< Pointer to the cache where the actual data is stored
		double * _QD;                       ///< Local array to hold the diagonal
};


/**
\class svrQMatrix
 \brief Kernel matrix for a support vector regression problem
 
 Kernel data is internally stored in a cache, such that not the
 full matrix should be stored.
 */		
class svrQMatrix: public Kernel
{ 
public:
		/**
	\brief Constructor of the kernel based on the data and the parameters
		 \param prob Reference to the training data of the SVM
		 \param param Reference to the set of SVM parameters
		 */
	svrQMatrix(SVMProblem & prob, const SVMParameters & param);
	
	/**
		\brief Destructor
	 */
	~svrQMatrix();
	
	/**
		\brief Get the 
		\param i index of the row
		\param len 
	 */
	double * getQ(unsigned int i, unsigned int len) const;

	/**
		\brief Get the diagonal of the kernel matrix
	 */
	double * getQD() const { return _QD; };
	
	/** 
		\brief Method to swap two columns in the kernel matrix
		\param i Index of the first column
		\param j Index of the second column
	 */	
	void swapIndex(unsigned int i, unsigned int j);
	void printIndex();
	
private:
	int _l;                             ///< Number of data  points in problem
	std::vector<unsigned int> _index;   ///< Local storage of the column indices
	Cache * _cache;                     ///< Pointer to the cache where the actual data is stored
	double * _QD;                       ///< Local array to hold the diagonal
	
	signed char * _sign;
	mutable int _nextBuffer;
	double * _buffer[2];
};


};

#endif  __SIMATH_SVM_KERNEL__
