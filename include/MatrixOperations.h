/*
 *  MatrixOperations.h
 *
 *  Created by Gert Thijs on 20/05/06.
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

#ifndef __SIMATH_MATRIXOPERATIONS__
#define __SIMATH_MATRIXOPERATIONS__

#include "Matrix.h"
#include "Vector.h"


namespace SiMath
{	
	/**
	 * \defgroup _matrix_ops_ Matrix Operations
	 * \brief Collection of specific matrix operations
	 */
	
	//----------------------------------------------------------------------------
	///\name Matrix statistics and inspection
	//@{
	
	/**
	 * \brief Get the maximal value of a matrix.
	 * \ingroup _matrix_ops_
	 *
	 * \param A Reference to the matrix
	 * \return The maximal value as a double
	 */
	double max(const SiMath::Matrix& A);
		
	/**
	 * \brief Get the maximal value of a given row.
	 * \ingroup _matrix_ops_
	 *
	 * \param A Reference to the matrix
	 * \param r Index of the row
	 * \return The maximal value of row \c r as a double
	 */
	double maxOfRow(const SiMath::Matrix& A, const unsigned int r);
	
	/**
	 * \brief Get the maximum values of all rows.
	 * \ingroup _matrix_ops_
	 *
	 * \param A Reference to the matrix
	 * \return Vector with for each row the corresponding maximum value
	 */
	SiMath::Vector maxOfAllRows(const SiMath::Matrix& A);
	
	/**
	 * \brief Get the maximal value of a given column.
	 * \ingroup _matrix_ops_
	 *
	 * \param A Reference to the matrix
	 * \param c Index of the column
	 * \return The maximal value of column \c c as a double
	 */
	double maxOfColumn(const SiMath::Matrix& A, const unsigned int c);
	
	/**
	 * \brief Get the maximum values of all columns.
	 * \ingroup _matrix_ops_
	 *
	 * \param A Reference to the matrix
	 * \return Vector with for each column the corresponding maximum value
	 */
	SiMath::Vector maxOfAllColumns(const SiMath::Matrix& A);
	
	/**
	 * \brief Get the minimal value of a matrix.
	 * \ingroup _matrix_ops_
	 *
	 * \param A Reference to the matrix
	 * \return The minimal value as a double
	 */
	double min(const SiMath::Matrix& A);
	
	/**
	 * \brief Get the minimal value of a given row.
	 * \ingroup _matrix_ops_
	 *
	 * \param A Reference to the matrix
	 * \param r Index of the row
	 * \return The minimal value of row \c r as a double
	 */
	double minOfRow(const SiMath::Matrix& A, const unsigned int r);
	
	/**
	 * \brief Get the minimal values of all rows.
	 * \ingroup _matrix_ops_
	 *
	 * \param A Reference to the matrix
	 * \return Vector with for each row the corresponding minimal value
	 */
	SiMath::Vector minOfAllRows(const SiMath::Matrix& A);
	
	/**
	 * \brief Get the minimal value of a given column.
	 * \ingroup _matrix_ops_
	 *
	 * \param A Reference to the matrix
	 * \param c Index of the column
	 * \return The minimal value of column \c c as a double
	 */
	double minOfColumn(const SiMath::Matrix& A, const unsigned int c);
	
	/**
	 * \brief Get the minimal values of all columns.
	 * \ingroup _matrix_ops_
	 *
	 * \param A Reference to the matrix
	 * \return Vector with for each column the corresponding minimal value
	 */
	SiMath::Vector minOfAllColumns(const SiMath::Matrix& A);
	
	/**
	 * \brief Get the sum over all elements of a matrix.
	 * \ingroup _matrix_ops_
	 *
	 * \param A Reference to the matrix
	 * \return The sum as a double
	 */
	double sum(const SiMath::Matrix& A);
	
	/**
	 * \brief Get the sum of a given row.
	 * \ingroup _matrix_ops_
	 *
	 * \param A Reference to the matrix
	 * \param r Index of the row
	 * \return The sum of row \c r as a double
	 */
	double sumOfRow(const SiMath::Matrix& A, const unsigned int r);
	
	/**
	 * \brief Get the sum of all rows.
	 * \ingroup _matrix_ops_
	 *
	 * \param A Reference to the matrix
	 * \return Vector with for each row the corresponding sum
	 */
	SiMath::Vector sumOfAllRows(const SiMath::Matrix& A);
	
	/**
	 * \brief Get the sum of a given column.
	 * \ingroup _matrix_ops_
	 *
	 * \param A Reference to the matrix
	 * \param c Index of the column
	 * \return The sum of column \c c as a double
	 */
	double sumOfColumn(const SiMath::Matrix& A, const unsigned int c);
	
	/**
	 * \brief Get the sum of all columns.
	 * \ingroup _matrix_ops_
	 *
	 * \param A Reference to the matrix
	 * \return Vector with for each column the corresponding sum
	 */
	SiMath::Vector sumOfAllColumns(const SiMath::Matrix& A);
	
	/**
	 * \brief Get the mean value of a matrix.
	 * \ingroup _matrix_ops_
	 *
	 * \param A Reference to the matrix
	 * \return The mean value as a double
	 */
	double mean(const SiMath::Matrix& A);
	
	/**
	 * \brief Get the mean value of a given row.
	 * \ingroup _matrix_ops_
	 *
	 * \param A Reference to the matrix
	 * \param r Index of the row
	 * \return The mean value of row \c r as a double
	 */
	double meanOfRow(const SiMath::Matrix& A, const unsigned int r);
	
	/**
	 * \brief Get the mean values of all rows.
	 * \ingroup _matrix_ops_
	 *
	 * \param A Reference to the matrix
	 * \return Vector with for each row the corresponding mean value
	 */
	SiMath::Vector meanOfAllRows(const SiMath::Matrix& A);
	
	/**
	 * \brief Get the mean value of a given column.
	 * \ingroup _matrix_ops_
	 *
	 * \param A Reference to the matrix
	 * \param c Index of the column
	 * \return The mean value of column \c c as a double
	 */
	double meanOfColumn(const SiMath::Matrix& A, const unsigned int c);
	
	/**
	 * \brief Get the mean values of all columns.
	 * \ingroup _matrix_ops_
	 *
	 * \param A Reference to the matrix
	 * \return Vector with for each column the corresponding mean value
	 */
	SiMath::Vector meanOfAllColumns(const SiMath::Matrix& A);
	
	/**
	 * \brief Get the standard deviation of a matrix.
	 * \ingroup _matrix_ops_
	 *
	 * \param A Reference to the matrix
	 * \return The standard deviation as a double
	 */
	double stDev(const SiMath::Matrix& A);
	
	/**
	 * \brief Get the standard deviation of a given row.
	 * \ingroup _matrix_ops_
	 *
	 * \param A Reference to the matrix
	 * \param r Index of the row
	 * \return The standard deviation of row \c r as a double
	 */
	double stDevOfRow(const SiMath::Matrix& A, const unsigned int r);
	
	/**
	 * \brief Get the standard deviation of all rows.
	 * \ingroup _matrix_ops_
	 *
	 * \param A Reference to the matrix
	 * \return Vector with for each row the corresponding standard deviation
	 */
	SiMath::Vector stDevOfAllRows(const SiMath::Matrix& A);
	
	/**
	 * \brief Get the standard deviation of a given column.
	 * \ingroup _matrix_ops_
	 *
	 * \param A Reference to the matrix
	 * \param c Index of the column
	 * \return The standard deviation of column \c c as a double
	 */
	double stDevOfColumn(const SiMath::Matrix& A, const unsigned int c);
	
	/**
	 * \brief Get the standard deviation of all columns.
	 \ingroup _matrix_ops_
	 *
	 * \param A Reference to the matrix
	 * \return Vector with for each column the corresponding stanard deviation
	 */
	SiMath::Vector stDevOfAllColumns(const SiMath::Matrix& A);
	
	//@}
	
	/** \name Matrix in place adaptors */
	//----------------------------------------------------------------------------

	//@{
	
	/**
		\brief Take the absolute value of all elements in the matrix
	 \ingroup _matrix_ops_
	  \param A reference to Matrix to be edited.
	 */
	void abs(SiMath::Matrix & A);
	
	/**
		\brief Take the d-th power of all elements in the matrix
	 \ingroup _matrix_ops_
	  \param A reference to SiMath::Matrix to be edited.
	  \param d double value 
	  \return Nothing
	 */
	void pow(SiMath::Matrix & A, double d);
	
	/**
		\brief Take the n-th power of all elements in the matrix
	 \ingroup _matrix_ops_
	  \param A reference to SiMath::Matrix to be edited.
	  \param n power  
	  \return Nothing
	 */
	void pow(SiMath::Matrix & A, unsigned int n);
	
	/**
		\brief Normalise each element per column according to a scale and offset
	 \ingroup _matrix_ops_
		
		Each element A[i][j] will be replace by (A[i][j] - offset[j])/scale[j]
	 
		\param A Matrix to normalise, will be changed in place 
		\param offset Vector with offset values should contains as many elements as there are columns in A
		\param scale Vector with scaling values should contains as many elements as there are columns in A
	 \return nothing
	 \exception Throws MatrixError if dimension of offset or scale differ from A	 
	 */
	void columnNormalise(SiMath::Matrix & A, const Vector & offset, const Vector & scale);

	/**
		\brief Normalise each element per row according to a scale and offset
	 \ingroup _matrix_ops_
	 
	 Each element A[i][j] will be replace by (A[i][j] - offset[j])/scale[j]
	 
	 \param A Matrix to normalise, will be changed in place 
	 \param offset Vector with offset values should contains as many elements as there are rows in A
	 \param scale Vector with scaling values should contains as many elements as there are rows in A
	 \return nothing
	 \exception Throws MatrixError if dimension of offset or scale differ from A
	 */
	void rowNormalise(SiMath::Matrix & A, const Vector & offset, const Vector & scale);
		//@}
	
	/** \name Matrix computations */
	//----------------------------------------------------------------------------

	//@{
	
	/**
		\brief Returns the transpose of the matrix. 
	 \ingroup _matrix_ops_
		\exception none.
		\return SiMath::Matrix
	*/ 
	SiMath::Matrix transpose(const SiMath::Matrix & A);
	
	/**
		\brief Function to compute matrix inversion based on Gauss-Jordan elimination.
	 \ingroup _matrix_ops_
		\param A SiMath::Matrix reference which is changed in place
		\exception SiMath::MatrixError parameter is out of range.
		\return void.
	 */ 
	void gaussJordanInverse(SiMath::Matrix & A);
	
	
	/**
		\brief Function to compute matrix inversion based on LU-decomposition. 
	 \ingroup _matrix_ops_
	 
	 luInverse is a function to compute the inverse of a (square) matrix A. luCompose calls internally first luDecompose to compute the LU 
	 decomposition of the matrix A. This decomposition is then used to compute the inverse
	 of A by using the luSolve function with the unit matrix (with same dimension as A) as input.
	 
	 \param     A   SiMath::Matrix reference which is changed in place.
	 \exception SiMath::MatrixError() when matrix A is not square.
	 \exception SiMath::MatrixError() when matrix A is singular.
	 \return void.
	 */ 
	void luInverse(SiMath::Matrix & A);
	
	/**
		\brief Function to compute matrix inversion based on LU-decomposition. 
	 \ingroup _matrix_ops_
	 
	 luInverse is a function to compute the inverse of a (square) matrix A. luCompose calls internally first luDecompose to compute the LU 
	 decomposition of the matrix A. This decomposition is then used to compute the inverse
	 of A by using the luSolve function with the unit matrix (with same dimension as A) as input.
	 
	 \param     A   SiMath::Matrix<double> reference.
	 \param     B   SiMath::Matrix<double> reference which is updated with the inverse of A.
	 \exception SiMath::MatrixError() when matrix A is not square.
	 \exception SiMath::MatrixError() when matrix A is singular.
	 \return    Nothing
	 */ 
	void luInverse(SiMath::Matrix & A, SiMath::Matrix & B);
	
	
	/**
		\brief Function to compute LU-decomposition of matrix A.
	 \ingroup _matrix_ops_
	 
	 Compute the LU decomposition of a matrix A. De lower (L) and upper (U) triangular matrices are stored in the matrix A. The sequence of row permutations is stored in I. 
	 
	 \param     A SiMath::Matrix reference is changed in place and contains L and U after decomposition 
	 \param     I Index vector with order of permuted rows
	 \exception SiMath::MatrixError() when matrix A is singular.
	 \return    Nothing
	 */ 
	void luDecompose(SiMath::Matrix & A, std::vector<int> & I);
	
	
	/**
		\brief Function to solve the equation Ax = B.
	 \ingroup _matrix_ops_
	 
	 This function solves the equation Ax = B via LU-decomposition of A.
	 In the forward substitution step Ly = B is solved. Next in the backsubstitution step 
	 Ux = y is solved. The resulting solutions  are stored in the matrix B. The vector I contains the sequence of 
	 row permutations as compute with luDecompose()
	 
	 \param     A  SiMath::Matrix reference is changed in place and contains L and U after 
	 \param     I  Index vector with order of permuted rows
	 \param     B  Basic SiMath::Matrix reference which is changed in place and contains the resulting X. 
	 \exception SiMath::MatrixError() when parameter is out of range.
	 \return    Nothing
	 */ 
	void luSolve(SiMath::Matrix & A, std::vector<int> & I, SiMath::Matrix & B);
	
	/**
		\brief Function to compute the Jacobi Transformation of a symmetric matrix
		\param A Input Matrix, is changed in place
		\param eigD Vector to holds eigen values
		\param eigV Matrix to hold sorted eigenvectors
		\param nRot maximal number of rotations. Contains actual number of rotations after completion.
	 */
	void jacobi (Matrix & A, Vector & eigD, Matrix & eigV, unsigned int & nRot);

	
	/**
	 * \brief  Function to compute product between SiMath::Matrix A[mxn] and Vector U[nx1] 
	 * \ingroup _matrix_ops_
	 * 
	 * This computes the product between each row of A and the Vector U.
	 * 
	 * \param  A SiMath::Matrix of dimensions [mxn]
	 * \param  U Vector of size n
	 * \returns Vector of size m
	 * \throws SiMath::MatrixError exception if the dimensions of A and U do not match
	 */
	Vector rowProduct(const SiMath::Matrix & A, const Vector & U);
	
	/**
	 * \brief  Function to compute product between Vector U[1xm] and SiMath::Matrix A[mxn] 
	 * \ingroup _matrix_ops_
	 * 
	 * This computes the product between each column of SiMath::Matrix A and the Vector U.
	 * 
	 * \param  U reference to Vector of size m
	 * \param  A reference to SiMath::Matrix of dimensions [mxn]
	 * \returns Vector of size n
	 * \throws SiMath::MatrixError exception if the dimensions of A and U do not match
	 */
	Vector colProduct(const Vector & U, const SiMath::Matrix & A);
		
	/** 
	 * \brief Computes the product of two matrices A and B and returns the result
	 * \ingroup _matrix_ops_
	 *  \param A reference to [NxM] matrix.
	 *  \param B reference to [MxP] matrix.
	 *  \exception SiMath::MatrixError If dimension of A and B do not match for matrix product computation
	 *  \return Newly created product matrix of dimension [NxP].
	 */ 
	SiMath::Matrix product(const SiMath::Matrix & A, const SiMath::Matrix & B);
	

	/**
	 * \brief Computes the product of A and B and stores the results in C
	 * \ingroup _matrix_ops_
	 * \param A reference to [NxM] matrix.
	 * \param B reference to [MxP] matrix.
	 * \param C reference to matrix to store the resulting product. C will be of dimension [NxP]
	 * \exception SiMath::MatrixError If dimension of A and B do not match for matrix product computation
	 */
	void largeProduct(const SiMath::Matrix & A, const SiMath::Matrix & B, SiMath::Matrix & C);

	/**
	 * \brief Compute correlation between variables in A
	 * \ingroup _matrix_ops_
* 	
		If A is and [MxN] matrix, the correlation matrix will be a [NxN] matrix C with elements
		\f[
			C[i][j] = \frac{\sum_{k=1}^M (A[k][i] - M_i )*( A[k][j] - M_j )}{\sqrt{\sum_k(A[k][i] - M_i )^2 \sum_k( A[k][j] - M_j)^2}} 
		\f]
	 * with Mi and Mj the mean of the i-th and j-th column respectively.
	 * \param A reference to matrix
	 * \return Correlation matrix of dimension [NxN]
	 */
	 SiMath::Matrix correlation(const SiMath::Matrix & A);
	
	/**
	 * \brief Find the shortest path between all pairs of elements based on a distance matrix
	 * \ingroup _matrix_ops_
	 * 
	 * Function to compute the shortest path between all pairs of elements  
	 * when the distances between elements are stored in the input matrix A.
	 * The path information is stored in an additional matrix P.
	 * 
	 * \param     A  Square matrix which contains distance between elements with index i and j.
	 * \param     P  Final Path matrix
	 * \return    Final cost matrix contains the cost of the shortest path between elements with index i and j
	 * \exception Nothing
	 */
	SiMath::Matrix allPairShortestPath(const SiMath::Matrix & A, SiMath::Matrix & P);
	
	
	/**
	 * \brief Find the shortest path between all pairs of elements based on a distance matrix.
	 * \ingroup _matrix_ops_
	 * 
	 * Function to compute the all-pair shortest path from the distances between elements 
	 * stored in the matrix A.
	 * 
	 * \param     A Reference to a matrix which holds the distances.
	 * \return    Final cost matrix contains the cost of the shortest path between elements with index i and j
	 * \exception Nothing
	 */
	SiMath::Matrix allPairShortestPath(const SiMath::Matrix & A);
	
	//@}

};


#endif __SIMATH_MATRIXOPERATIONS__

