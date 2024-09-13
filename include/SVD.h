/*
 *  SVD.h
 *
 *  Created by Gert Thijs on 19/06/06.
 *  Copyright 2006 Silicos NV. All rights reserved.
 *
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

#ifndef __SIMATH_SVD__
#define __SIMATH_SVD__

#include "Matrix.h"
#include "Vector.h"


namespace SiMath
{
	/**
		\class  SVDError
		\brief  Exception class
	 
		This exception class is derived from the std::runtime_error class and is used to 
		indicate that something went wrong during the SVD computation.
	 \ingroup _error_handling_
	 */	
	class SIMATH_EXCEPTION SVDError : public SiMathError
	{
		public:
		SVDError(const std::string & e) : SiMathError("SVDError : " + e) {};
	};

	/**
		\class SVD
		\brief Class to compute the Singular Value Decomposition of a matrix.
	 
	  Based on the C++ code from TNT and JAMA (more info: http://math.nist.gov/tnt/index.html)
	 */
	class SVD
	{
		public: 
			/**
				\brief Constructor processes matrix and computes SVD ( A=USV')
				
				By initialiation of the SVD object the original matrix is copied 
			  and this copy is used for further computations. By default both 
			  U and V are computed and stored internally. The computation of U or V 
			  can be switched off by using the appropriate status flag.
			 
				\param A data matrix 
				\param bU defines whether or not the U matrix should be computed
				\param bV defines whether or not the V matrix should be computed
			 */
			SVD(const Matrix & A, bool bU = true, bool bV = true);
		
			/**
				\brief Get the singular values stored in a Vector.
			 */
			Vector getSingularValues() {return _S;};
			
			/** 
				\brief Get the singular values as a matrix 
			 */
			Matrix getSingularMatrix();
			
			/**
				\brief Get the left singular vectors 
				\return Matrix of dimension [mxm]
			 */
			Matrix getU() { return _U;};
			
			/**
				\brief Get the right singular vectors 
				\return Matrix of dimension [nxn], with n the number of singular values
			 */
			Matrix getV() { return _V;};

			/**
				\brief Get the greatest singular value
			 */
			double norm2() { return _S[0];}; 
			
			/** 
				\brief Two norm of condition number (max(S)/min(S))
				*/
			double cond() { return _S[0]/_S[_S.size()-1]; };

			/**
				\brief Get rank of the problem
			 */
			int rank();

		private:
			int _m;      ///< number of rows
			int _n;      ///< number of columns
			Matrix _U;   ///< Left singular vectors
			Matrix _V;   ///< Right singular vectors
			Vector _S;   ///< Singular values
			
			bool _computeV;   ///< Check if V should be computed
			bool _computeU;   ///< Check if U should be computed
			
	};
};


#endif __SIMATH_SVD__
