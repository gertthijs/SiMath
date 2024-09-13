/*
 *  Matrix.h
 *
 *  Created by Gert Thijs on 24/03/06.
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


#ifndef  __SIMATH_MATRIX__
#define  __SIMATH_MATRIX__

#include <iostream>
#include "Definitions.h"
#include "Vector.h"

namespace SiMath
{
	
	/**
		\class MatrixError
	 \brief  Exception class
	 
	 This exception class is derived from the std::runtime_error class and is used to 
	 indicate that something went wrong when manipulating a Matrix.
	 \ingroup _error_handling_
	 */	
	class SIMATH_EXCEPTION MatrixError : public SiMathError
	{
		public:
			MatrixError(const std::string & e) : SiMathError("MatrixError : " + e ) {};
	};


  /**
		\class  Matrix 
		\brief  Class for storage of data and basic matrix operators
	 
		This is the implementation of the basic matrix class used in SiMath.
	 
		The following example shows how to create a matrix and how to manipulate it.
		\code 
		#include "SiMath.h"
	 
		SiMath::Matrix M(n,m,0.0);  // create new [nxm] matrix M with all elements equal to 0
		M += 1.0;                               // add 1.0 to all elements in M
		SiMath::Matrix A(M);        // use the copy constructor to make a copy of M
		for ( j=0; j<M.nbrColumns(); j++)
			A[i][j] = rand();                     // fill indivudual elements 
	 \endcode
	 
   */
  class Matrix
	{
		private:
		//@{
			unsigned int _nRows;   ///< Number of rows of matrix
			unsigned int _nCols;   ///< Number of columns of matrix
			double** _pMatrix;   ///< internal matrix is consists of an array of pointers to indices in the a matrix array
		//@}
	
		public:
		/// \name Structors
    //@{
			Matrix() : _nRows(0), _nCols(0), _pMatrix(NULL) {};                ///< Empty matrix
			Matrix(const unsigned int n, const unsigned int m);                ///< NxM matrix, no initial value
			Matrix(const unsigned int n, const unsigned int m, const double & v);   ///< NxM matrix, initial constant value v
			Matrix(const Matrix & M);                                          ///< Copy constructor
	
			/**
			  \brief Initialisation of NxM matrix from vector.
			  \throw MatrixError() if vector is not of length N*M
			 */
			Matrix(const unsigned int n, const unsigned int m, const Vector & vec);
	
			~Matrix();    ///< destructor
		//@}
	
		/// \name Dimension manipulation
		//@{
			/**
				\brief Reset dimension of matrix to [RxC] and fill matrix with 0
			 
			 \param r Number of rows
			 \param c Number of columns
			 \return void
			 */
			void reset (const unsigned int r, const unsigned int c);
			
			/**
				\brief Clear the content of the full matrix and make it empty.
			 */
			void clear();
			
		//@}
			
		/// \name Data Accessors
    //@{
			/** 
				\brief Number of rows
		   */
			inline unsigned int nbrRows() const { return _nRows;};
			
			/**
				\brief Number of columns
			 */
			inline unsigned int nbrColumns() const { return _nCols;};
	
			/**
			 \brief Get element [i,j] with checks on range of i and j 
	 
			 This is the slower but more robust method to get a specific element from the matrix. 
			 The indices i and j are checked and an error is thrown if they are out of range.
	 
			 \param  i  row index
			 \param  j  column index
			 \exception MatrixError()
			 */
			double getValueAt(const unsigned int i, const unsigned int j);
	
			/**
			 \brief Get element (i,j) with checks on range of i and j, const implementation
	 
			 This is the slower but more robust method to get a specific element from the matrix. 
			 The indices i and j are checked and an error is thrown if they are out of range.
	 
			 \param i row index
			 \param j column index
			 \exception MatrixError()
			 */
			const double getValueAt(const unsigned int i, const unsigned int j) const;
	
			/**
				\brief Get the i-th row as a SiMath::Vector

				This method creates a new Vector and copies the values of the i-th row into this Vector.
				Manipulation on the generated Vector do not affect the original matrix.
				\param i row index
			 */
			SiMath::Vector getRow(const unsigned int i) const;
			
			/**
				\brief Get the i-th column as a SiMath::Vector
				
				This method creates a new Vector and copies the values of the i-th column into this Vector
				Manipulation on the generated Vector do not affect the original matrix.
				\param i column index
			 */			
			SiMath::Vector getColumn(const unsigned int i) const;
			
		//@}			
			
		/// \name Data Adaptors
		//@{
			/**
				\brief Save method to set element (i,j) to value v
	 
			 This is the slower but more robust method to set a specific element in the matrix. 
			 The indices i and j are checked and an error is thrown if they are out of range.
	 
			 \param i row index
			 \param j column index
			 \param v value to be set
			 \exception MatrixError() 
			 \return void
			 */
			void setValueAt(const unsigned int i, const unsigned int j, double v);
	
			/**
				\brief Set the values of a row from a SiMath::Vector
				\param i row index
				\param src SiMath::Vector to be copied.
			 */
			void setRow(const unsigned int i, Vector & src);

			/**
				\brief Set the values of a column from a SiMath::Vector
				\param i column index
				\param src SiMath::Vector to be copied.
			 */
			void setColumn(const unsigned int i, Vector & src);
			
			/**
				\brief swap rows i and j in the matrix
				\param i row index
				\param j column index
				\return nothing
			 */
			void swapRows(unsigned int i, unsigned int j);

			/**
				\brief swap columns i and j in the matrix
				\param i row index
				\param j column index
				\return nothing
			 */
			void swapColumns(unsigned int i, unsigned int j);
			
		//@}
	
		///\name Arithmetic Operators
    //@{
			Matrix & operator= (const Matrix & M);      ///< copy assignment, resets the size of the matrix if needed
			Matrix & operator= (const double & v);      ///< set all elements in matrix to constant value
			Matrix & operator+= (const double & v);     ///< add constant value to all elements in matrix
			Matrix & operator+= (const Matrix & M);     ///< add full matrix element-wise
			Matrix & operator-= (const double & v);     ///< subtract constant value from all elements in matrix
			Matrix & operator-= (const Matrix & M);     ///< subtract full matrix element-wise
			Matrix & operator*= (const double & v);     ///< multiply all elements with a constant value
			Matrix & operator*= (const Matrix & M);     ///< multiply full matrix element-wise 
			Matrix & operator/= (const double & v);     ///< divide all elements with a constant value
			Matrix & operator/= (const Matrix & M);     ///< divide full matrix element-wise 
			Matrix & operator- ();                      ///< change sign of all elements in matrix
	
			// element-wise operators
			Matrix operator+ (const Matrix & M) const;  ///< add two matrices element by element and store the result in a new matrix
			Matrix operator- (const Matrix & M) const;  ///< substract two matrices element by element and store the result in a new matrix
			Matrix operator* (const Matrix & M) const;  ///< multiply two matrices element by element and store the result in a new matrix
			Matrix operator/ (const Matrix & M) const;  ///< divide two matrices element by element and store the result in a new matrix
			
			//  indexing
			/**
				\brief Direct access to row through a pointer
				\param i Row index
				\return Pointer to i-th row
			 */
			inline double * operator[] (const unsigned int i) { return _pMatrix[i]; };

			/**
				\brief Direct access to row through a pointer (const implementation)
				\param i Row index
				\return Pointer to i-th row
			 */
			inline const double * operator[] (const unsigned int i) const { return _pMatrix[i]; };
		//@}
	
			std::ostream & print(std::ostream & os) const;    ///< Debugging method to print matrix in a readable format
	};	
	
	/** \brief Output stream operator
		
		\param os  std::ostream &
		\param A   Matrix reference
		\return std::ostream reference
		*/
	// std::ostream & operator<< (std::ostream & os, const SiMath::Matrix & A) {	return A.print(os); };
	std::ostream & operator<< (std::ostream & os, const SiMath::Matrix & A);
	
};  // end of namespace declaration


#endif   __SIMATH_MATRIX__
