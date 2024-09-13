/*
 *  MatrixOperations.cpp
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

#include "Matrix.h"
#include "MatrixOperations.h"
#include "Utilities.h"


//----------------------------------------------------------------------------//
// max																																				//
//----------------------------------------------------------------------------//
double SiMath::max(const SiMath::Matrix& A)
{
	unsigned int r(A.nbrRows());
	unsigned int c(A.nbrColumns());
	double m(A[0][0]);
	
	for(int i=0 ; i<r ; ++i) {
		for(int j=0 ; j<c ; ++j) {
			if(A[i][j] > m)
				m = A[i][j];
		}
	}
	
	return m;
}

//----------------------------------------------------------------------------//
// maxOfRow																																		//
//----------------------------------------------------------------------------//
double SiMath::maxOfRow(const SiMath::Matrix& A, const unsigned int r)
{
	if(r >= A.nbrRows())
		throw MatrixError("Row index out of bounds.");

	double m(A[r][0]);
	unsigned int c(A.nbrColumns());
		
	for(int j=1 ; j<c ; ++j) {	
		if(A[r][j] > m) 
			m = A[r][j];
	}
	
	return m;
}

//----------------------------------------------------------------------------//
// maxOfAllRows																																//
//----------------------------------------------------------------------------//
SiMath::Vector SiMath::maxOfAllRows(const SiMath::Matrix& A)
{
	unsigned int r(A.nbrRows());
	unsigned int c(A.nbrColumns());
	
	Vector vec(r);
	for(int i=0 ; i<r ; ++i) {
		vec[i] = A[i][0];
	}
	
	for(int i=0 ; i<r ; ++i) {
		for(int j=1 ; j<c ; ++j) {
			if(A[i][j] > vec[i])
				vec[i] = A[i][j];
		}
	}
	
	return vec;
}

//----------------------------------------------------------------------------//
// maxOfColumn																																//
//----------------------------------------------------------------------------//
double SiMath::maxOfColumn(const SiMath::Matrix& A, const unsigned int c)
{
	if(c >= A.nbrColumns())
		throw MatrixError("Column index out of bounds.");

	double m(A[0][c]);
	unsigned int r(A.nbrRows());
		
	for(int i=1 ; i<r ; ++i) {	
		if(A[i][c] > m) 
			m = A[i][c];
	}
	
	return m;
}

//----------------------------------------------------------------------------//
// maxOfAllColumns																														//
//----------------------------------------------------------------------------//
SiMath::Vector SiMath::maxOfAllColumns(const SiMath::Matrix& A)
{
	unsigned int r(A.nbrRows());
	unsigned int c(A.nbrColumns());
	
	Vector vec(c);
	for(int j=0 ; j<c ; ++j) {
		vec[j] = A[0][j];
	}
	
	for(int j=0 ; j<c ; ++j) {
		for(int i=1 ; i<r ; ++i) {
			if(A[i][j] > vec[j])
				vec[j] = A[i][j];
		}
	}
	
	return vec;
}

//----------------------------------------------------------------------------//
// min																																				//
//----------------------------------------------------------------------------//
double SiMath::min(const SiMath::Matrix& A)
{
	unsigned int r(A.nbrRows());
	unsigned int c(A.nbrColumns());
	double m(A[0][0]);
	
	for(int i=0 ; i<r ; ++i) {
		for(int j=0 ; j<c ; ++j) {
			if(A[i][j] < m)
				m = A[i][j];
		}
	}
	
	return m;
}

//----------------------------------------------------------------------------//
// minOfRow																																		//
//----------------------------------------------------------------------------//
double SiMath::minOfRow(const SiMath::Matrix& A, const unsigned int r)
{
	if(r >= A.nbrRows())
		throw MatrixError("Row index out of bounds.");
	
	double m(A[r][0]);
	unsigned int c(A.nbrColumns());
		
	for(int j=1 ; j<c ; ++j) {	
		if(A[r][j] < m) 
			m = A[r][j];
	}
	
	return m;
}

//----------------------------------------------------------------------------//
// minOfAllRows																																//
//----------------------------------------------------------------------------//
SiMath::Vector SiMath::minOfAllRows(const SiMath::Matrix& A)
{
	unsigned int r(A.nbrRows());
	unsigned int c(A.nbrColumns());
	
	Vector vec(r);
	for(int i=0 ; i<r ; ++i) {
		vec[i] = A[i][0];
	}
	
	for(int i=0 ; i<r ; ++i) {
		for(int j=1 ; j<c ; ++j) {
			if(A[i][j] < vec[i])
				vec[i] = A[i][j];
		}
	}
	
	return vec;
}

//----------------------------------------------------------------------------//
// minOfColumn																																//
//----------------------------------------------------------------------------//
double SiMath::minOfColumn(const SiMath::Matrix& A, const unsigned int c)
{
	if(c >= A.nbrColumns())
		throw MatrixError("Column index out of bounds.");
	
	double m(A[0][c]);
	unsigned int r(A.nbrRows());
		
	for(int i=1 ; i<r ; ++i) {	
		if(A[i][c] < m) 
			m = A[i][c];
	}
	
	return m;
}

//----------------------------------------------------------------------------//
// minOfAllColumns																														//
//----------------------------------------------------------------------------//
SiMath::Vector SiMath::minOfAllColumns(const SiMath::Matrix& A)
{
	unsigned int r(A.nbrRows());
	unsigned int c(A.nbrColumns());
	
	Vector vec(c);
	for(int j=0 ; j<c ; ++j) {
		vec[j] = A[0][j];
	}
	
	for(int j=0 ; j<c ; ++j) {
		for(int i=1 ; i<r ; ++i) {
			if(A[i][j] < vec[j])
				vec[j] = A[i][j];
		}
	}
	
	return vec;
}

//----------------------------------------------------------------------------//
// sum																																				//
//----------------------------------------------------------------------------//
double SiMath::sum(const SiMath::Matrix& A)
{
	unsigned int r(A.nbrRows());
	unsigned int c(A.nbrColumns());
	double s(0.0);
	
	for(int i=0 ; i<r ; ++i) {
		for(int j=0 ; j<c ; ++j) {
			s += A[i][j];
		}
	}
	
	return s;
}

//----------------------------------------------------------------------------//
// sumOfRow																																		//
//----------------------------------------------------------------------------//
double SiMath::sumOfRow(const SiMath::Matrix& A, const unsigned int r)
{
	if(r >= A.nbrRows())
		throw MatrixError("Row index out of bounds.");
	
	double s(0.0);
	unsigned int c(A.nbrColumns());
		
	for(int j=0 ; j<c ; ++j) {	
		s += A[r][j];
	}
	
	return s;
}

//----------------------------------------------------------------------------//
// sumOfAllRows																																//
//----------------------------------------------------------------------------//
SiMath::Vector SiMath::sumOfAllRows(const SiMath::Matrix& A)
{
	unsigned int r(A.nbrRows());
	unsigned int c(A.nbrColumns());
	
	Vector vec(r);
	for(int i=0 ; i<r ; ++i) {
		vec[i] = A[i][0];
	}
	
	for(int i=0 ; i<r ; ++i) {
		for(int j=1 ; j<c ; ++j) {
			vec[i] += A[i][j];
		}
	}
	
	return vec;
}

//----------------------------------------------------------------------------//
// sumOfColumn																																//
//----------------------------------------------------------------------------//
double SiMath::sumOfColumn(const SiMath::Matrix& A, const unsigned int c)
{
	if(c >= A.nbrColumns())
		throw MatrixError("Column index out of bounds.");
	
	double s(0.0);
	unsigned int r(A.nbrRows());
		
	for(int i=0 ; i<r ; ++i) {	
		s += A[i][c];
	}
	
	return s;
}

//----------------------------------------------------------------------------//
// sumOfAllColumns																														//
//----------------------------------------------------------------------------//
SiMath::Vector SiMath::sumOfAllColumns(const SiMath::Matrix& A)
{
	unsigned int r(A.nbrRows());
	unsigned int c(A.nbrColumns());
	
	Vector vec(c);
	for(int j=0 ; j<c ; ++j) {
		vec[j] = A[0][j];
	}
	
	for(int j=0 ; j<c ; ++j) {
		for(int i=1 ; i<r ; ++i) {
			vec[j] += A[i][j];
		}
	}
	
	return vec;
}

//----------------------------------------------------------------------------//
// mean																																				//
//----------------------------------------------------------------------------//
double SiMath::mean(const SiMath::Matrix& A)
{
	unsigned int r(A.nbrRows());
	unsigned int c(A.nbrColumns());
	double m(0.0);
	
	for(int i=0 ; i<r ; ++i) {
		for(int j=0 ; j<c ; ++j) {
			m += A[i][j];
		}
	}
	
	return m / (r*c);
}

//----------------------------------------------------------------------------//
// meanOfRow																																	//
//----------------------------------------------------------------------------//
double SiMath::meanOfRow(const SiMath::Matrix& A, const unsigned int r)
{
	if(r >= A.nbrRows())
		throw MatrixError("Row index out of bounds.");

	double m(0.0);
	unsigned int c(A.nbrColumns());
	
	for(int j=0 ; j<c ; ++j) {	
		m += A[r][j];
	}
	
	return m/c;
}

//----------------------------------------------------------------------------//
// meanOfAllRows																															//
//----------------------------------------------------------------------------//
SiMath::Vector SiMath::meanOfAllRows(const SiMath::Matrix& A)
{	
	unsigned int r(A.nbrRows());
	unsigned int c(A.nbrColumns());
	SiMath::Vector vec(r, 0.0);
	
	for(int i=0 ; i<r ; ++i) {
		for(int j=0 ; j<c ; ++j) {
			vec[i] += A[i][j];
		}
	}
	
	for(int i=0 ; i<r ; ++i) {
		vec[i] /= c;
	}
	
	return vec;
}

//----------------------------------------------------------------------------//
// meanOfColumn																																//
//----------------------------------------------------------------------------//
double SiMath::meanOfColumn(const SiMath::Matrix& A, const unsigned int c)
{
	if(c >= A.nbrColumns())
		throw MatrixError("Column index out of bounds.");
	
	double m(0.0);
	unsigned int r(A.nbrRows());
	
	for(int i=0 ; i<r ; ++i) {	
		m += A[i][c];
	}
	
	return m/r;
}

//----------------------------------------------------------------------------//
// meanOfAllColumns																														//
//----------------------------------------------------------------------------//
SiMath::Vector SiMath::meanOfAllColumns(const SiMath::Matrix& A)
{	
	unsigned int r(A.nbrRows());
	unsigned int c(A.nbrColumns());
	SiMath::Vector vec(c, 0.0);
	
	for(int j=0 ; j<c ; ++j) {
		for(int i=0 ; i<r ; ++i) {
			vec[j] += A[i][j];
		}
	}
	
	for(int j=0 ; j<c ; ++j) {
		vec[j] /= r;
	}
	
	return vec;
}

//----------------------------------------------------------------------------//
// stDev																																			//
//----------------------------------------------------------------------------//
double SiMath::stDev(const SiMath::Matrix& A)
{
	unsigned int r(A.nbrRows());
	unsigned int c(A.nbrColumns());
	double sum(0.0);
	double std(0.0);
	
	for(int i=0 ; i<r ; ++i) {
		for(int j=0 ; j<c ; ++j) {
			sum += A[i][j];
			std += A[i][j] * A[i][j];
		}
	}
		
	int n(r*c);
	std -= sum * sum / n;
	std /= n-1;
	std = sqrt(std);
	
	return std;
}

//----------------------------------------------------------------------------//
// stDevOfRow																																	//
//----------------------------------------------------------------------------//
double SiMath::stDevOfRow(const SiMath::Matrix& A, const unsigned int r)
{
	if(r >= A.nbrRows())
		throw MatrixError("Row index out of bounds.");
	
	unsigned int c(A.nbrColumns());
	double sum(0.0);
	double std(0.0);
	
	for(int j=0 ; j<c ; ++j) {
		sum += A[r][j];
		std += A[r][j] * A[r][j];
	}
	
	std -= sum*sum/c;
	std /= c-1;
	std = sqrt(std);
	
	return std;
}
	
//----------------------------------------------------------------------------//
// stDevOfAllRows																															//
//----------------------------------------------------------------------------//
SiMath::Vector SiMath::stDevOfAllRows(const SiMath::Matrix& A)
{
	unsigned int r(A.nbrRows());
	unsigned int c(A.nbrColumns());
	SiMath::Vector sumVec(r, 0.0);
	SiMath::Vector stdVec(r, 0.0);
	
	for(int i=0 ; i<r ; ++i) {
		for(int j=0 ; j<c ; ++j) {
			sumVec[i] += A[i][j];
			stdVec[i] += A[i][j] * A[i][j];
		}
	}
	
	for(int i=0 ; i<r ; ++i) {
		stdVec[i] -= sumVec[i]*sumVec[i]/c;
		stdVec[i] /= c-1;
		stdVec[i] = sqrt(stdVec[i]);
	}
	
	return stdVec;
}
		
//----------------------------------------------------------------------------//
// stDevOfColumn																															//
//----------------------------------------------------------------------------//
double SiMath::stDevOfColumn(const SiMath::Matrix& A, const unsigned int c)
{
	if(c >= A.nbrColumns())
		throw MatrixError("Column index out of bounds.");
	
	unsigned int r(A.nbrRows());
	double sum(0.0);
	double std(0.0);
	
	for(int i=0 ; i<r ; ++i) {
		sum += A[i][c];
		std += A[i][c] * A[i][c];
	}
	
	std -= sum*(sum/r);
	std /= r-1;
	std = sqrt(std);
	
	return std;
}

//----------------------------------------------------------------------------//
// stDevOfAllColumns																													//
//----------------------------------------------------------------------------//
SiMath::Vector SiMath::stDevOfAllColumns(const SiMath::Matrix& A)
{
	unsigned int r(A.nbrRows());
	unsigned int c(A.nbrColumns());
	SiMath::Vector sumVec(c, 0.0);
	SiMath::Vector stdVec(c, 0.0);

	for(int j=0 ; j<c ; ++j) {
		for(int i=0 ; i<r ; ++i) {
			sumVec[j] += A[i][j];
			stdVec[j] += A[i][j] * A[i][j];
		}
	}
	
	for(int j=0 ; j<c ; ++j) {
		stdVec[j] -= sumVec[j]*sumVec[j]/r;
		stdVec[j] /= r-1;
		stdVec[j] = sqrt(stdVec[j]);
	}
	
	return stdVec;
}


//------------------------------------------------------------------------------------------------------------------------
SiMath::Matrix 
SiMath::transpose(const SiMath::Matrix & A)
{	
	SiMath::Matrix T(A.nbrColumns(),A.nbrRows());
	for( unsigned int i=0; i < A.nbrRows(); ++i )
	{
		for ( unsigned int j=0; j < A.nbrColumns(); ++j )
		{
			T[j][i] = A[i][j]; 
		}
	}
	
	return T;
}


//------------------------------------------------------------------------------------------------------------------------

void
SiMath::abs(SiMath::Matrix & A)
{
	for ( unsigned int i=0; i<A.nbrRows(); ++i)
		for ( unsigned int j=0; j<A.nbrColumns(); ++j)
			A[i][j] = fabs(A[i][j]);
	return;
}


//------------------------------------------------------------------------------------------------------------------------

void
SiMath::pow(SiMath::Matrix & A, double d)
{
	for ( unsigned int i=0; i<A.nbrRows(); ++i)
		for ( unsigned int j=0; j<A.nbrColumns(); ++j)
			A[i][j] = std::pow(A[i][j],d);
	return;
}


void
SiMath::pow(SiMath::Matrix & A, unsigned int d)
{
	for ( unsigned int i=0; i<A.nbrRows(); ++i)
		for ( unsigned int j=0; j<A.nbrColumns(); ++j)
			A[i][j] = SiMath::powi(A[i][j],d);
	return;
}


//------------------------------------------------------------------------------------------------------------------------
void 
SiMath::columnNormalise(SiMath::Matrix & A, const Vector & offset, const Vector & scale)
{
	if ( A.nbrColumns() != offset.size() || A.nbrColumns() != scale.size() )
		throw(SiMath::MatrixError("Dimension mismatch between matrix and offset and/or scaling"));
	
	for ( int j=0; j<A.nbrColumns(); ++j)
	{
		double m(offset[j]);
		double s(scale[j]);
		
		for ( int i=0; i<A.nbrRows(); ++i)
			A[i][j] = (A[i][j] - m)/s;
	}
	
	return;
}

//------------------------------------------------------------------------------------------------------------------------
void 
SiMath::rowNormalise(SiMath::Matrix & A, const Vector & offset, const Vector & scale)
{
	if ( A.nbrRows() != offset.size() || A.nbrRows() != scale.size() )
		throw(SiMath::MatrixError("Dimension mismatch between matrix and offset and/or scaling"));

	for ( int i=0; i<A.nbrRows(); ++i)
	{
		double m(offset[i]);
		double s(scale[i]);
		
		for ( int j=0; j<A.nbrColumns(); ++j)
			A[i][j] = (A[i][j] - m)/s;
	}
		
	return;
}




//------------------------------------------------------------------------------------------------------------------------------------
void 
SiMath::gaussJordanInverse(SiMath::Matrix & A)
{
	int n = A.nbrRows();
	int rowI = 0;
	int colI = 0;
	int j, ii, jj;
	// if ( A.NbrColumns() != n )
		// throw std::out_of_range;
	double dummy, maxVal, pivot;
	
	std::vector<int> pivotIndex(n,0), // initial pivot elements are zero
		rowIndex(n,0), 
		colIndex(n,0);
		
	for ( j=0; j<n; ++j)  // loop over columns to be reduced
	{
		maxVal = 0.0;
		// find pivot element 
		for ( ii=0; ii<n; ++ii)  // loop over rows to find next highest element  
		{
			if ( pivotIndex[ii] != 1 )  // check if row has been handled before
			{
				for ( jj=0; jj<n; ++jj)   // loop over colums
				{
					if ( pivotIndex[jj] == 0 && fabs(A[ii][jj]) >= maxVal ) // if 
					{
						maxVal = fabs(A[ii][jj]);
						rowI = ii;
						colI = jj;
					}
				}
			}
		}
		++(pivotIndex[colI]); // augment pivot index at position colI
		
		// now swap rows if needed, ie. the pivot element lies off the diagonal
		if ( rowI != colI )
			A.swapRows(rowI,colI);
		
		// store indices of swapped rows
		rowIndex[j] = rowI;
		colIndex[j] = colI;
		
		// check if it is a singular matrix
		//if ( A.GetValueAt(colI,colI) == 0.0 )
		//throw std::out_of_range;
		
		// compute the pivot element
		pivot = 1.0/A[colI][colI];
		A[colI][colI] = 1.0;
		
		// update row
		for ( jj=0; jj<n; ++jj)
			A[colI][jj] *= pivot;
		
		// reduce rows
		for ( ii=0; ii<n; ++ii )
		{
			if ( ii != colI )
			{
				dummy = A[ii][colI];
				A[ii][colI] = 0.0;
				for (jj=0; jj<n; ++jj)
					A[ii][jj] -= dummy*A[colI][jj];
			}
		}
		
	}
	
	// unscramble the matrix
	for ( j=n-1; j>=0; j--)
	{
		if ( rowIndex[j] != colIndex[j] )
			A.swapColumns(rowIndex[j],colIndex[j]);
	}
	
	return;
}	



void 
SiMath::luDecompose(SiMath::Matrix & A,std::vector<int> & I)
{
	int n = A.nbrRows();
	int i, j, k, kMax, iMax;
	std::vector<double> vScales(n,0);
	double maxVal = 0, dummy =0;
	double * pRowi = NULL;
	
	// first find the highest pivot element in each row and store it for implicit scaling
	for ( i=0; i<n; ++i )
	{
		maxVal = 0.0;
		for ( j=0; j<n; ++j)
		{
			if ( (dummy=fabs(A[i][j])) > maxVal) 
				maxVal = dummy;
		}
		if ( maxVal <= 1e-6 ){
			throw(SiMath::MatrixError("singular matrix"));
		}

		// std::cerr << "max val: " << maxVal << std::endl;
		vScales[i] = 1.0/maxVal;
	}
	
	double colJ[n]; // variable to store local copy of column
	
	// loop over columns 
	for ( j=0; j<n; ++j )
	{
		// make a local copy of column j
		for ( i=0; i<n; ++i )
			colJ[i] = A[i][j];
		
		for ( i=0; i<n; ++i )
		{
			pRowi = A[i];
			dummy = pRowi[j];
			kMax = i < j ? i : j;
			for ( k=0; k<kMax; ++k )
				dummy -= pRowi[k] * colJ[k];
			colJ[i] = dummy;
			pRowi[j] = colJ[i];
		}
		
		// search largest pivot element
		maxVal = 0.0;
		iMax = j;
		for ( i=j+1; i<n; ++i )
		{
			if ( (dummy = fabs(colJ[i]) * vScales[i]) >= maxVal )
			{
				maxVal = dummy;
				iMax = i;
			}
		}
		
		// check if we need to interchange rows
		if ( j != iMax ) // if current column index is not the maximal row index we need to interchange
		{
			// std::cerr << "Swap rows: " << iMax << " <-> " << j << std::endl;
			A.swapRows(iMax,j);
			vScales[iMax] = vScales[j];
		}
		// store row index in I
		I[j] = iMax;
		
		// finally divide by the pivot element
		if ( j != n-1 )
		{
			dummy = 1.0/A[j][j]; // A.GetValueAt(j,j);
			for ( i=j+1; i<n; ++i)
				A[i][j] *= dummy;
		}
		
		
	} // next column
	
	return;
}


void 
SiMath::luSolve(SiMath::Matrix & A, std::vector<int> & I, SiMath::Matrix & B)
{
	int n = A.nbrRows();
	int m = B.nbrColumns();
	int i, j, k;
	
	for (int i=0; i<n; ++i)
		B.swapRows(i,I[i]);
	
	// forward substitution pass
	for ( k=0; k<n; ++k )
		for ( i=k+1; i<n; ++i )
			for ( j=0; j<m; ++j ) // do this for all columns of B
				B[i][j] -= A[i][k] * B[k][j];
	
	// do the backsubstitution
	for ( i=n-1; i>=0; --i )
	{
		for ( j=0; j<m; ++j)  // divide all elements 
			B[i][j] /= A[i][i]; 
		
		for ( k=0; k<i; ++k )
			for ( j=0; j<m; ++j )
				B[k][j] -= A[k][i] * B[i][j];		
	}
	
	return;
}


void
SiMath::luInverse(SiMath::Matrix & A)
{
	int n = A.nbrRows();
	if ( n != A.nbrColumns() )
		throw(MatrixError("A is not square."));
	
	int i;
	std::vector<int> iPerm(n,0);  // vector to store permutation inidices
	
	// compute LU Decomposition
	SiMath::luDecompose(A,iPerm);
	
	// do the backsubstition for the unit matrix B
	SiMath::Matrix B(n,n,0);
	for ( i=0; i<n; ++i)
		B[i][i] = 1;
		
	// solve backsubstitution on b
	SiMath::luSolve(A,iPerm,B);
		
	// place results from B in A
	A = B;
	
	return;
	
}

void
SiMath::luInverse(SiMath::Matrix & A, SiMath::Matrix & B)
{
	unsigned int n = A.nbrRows();
	if ( n != A.nbrColumns() )
		throw(MatrixError("A is not square."));
	unsigned int i;
	std::vector<int> iPerm(n,0);  // vector to store permutation inidices
	
	if ( B.nbrRows() != n || B.nbrColumns() != n )
		B.reset(n,n);
	
	// compute LU Decomposition
	SiMath::luDecompose(A,iPerm);
	
	// do the backsubstition for the unit matrix B
	for ( i=0; i<n; ++i )
		B[i][i] = 1;
		
	// solve backsubstitution on b
	SiMath::luSolve(A,iPerm,B);
		
	return;
	
}

//--------------------------PRODUCT----------------------------

//--- compute A * U ----
SiMath::Vector
SiMath::rowProduct(const SiMath::Matrix & A, const SiMath::Vector & U)
{
	if ( U.size() != A.nbrColumns() )
		throw(MatrixError("Dimensions of matrix and vector are not compatible."));
	
	Vector v(A.nbrRows(),0.0);
	
	for ( unsigned int i=0; i<A.nbrRows(); ++i)
	{
		double s(0.0);
		for ( unsigned int j=0; j<A.nbrColumns(); ++j)
		{
			s += A[i][j] * U[j];
		}
		v[i] = s;
	}
	return v;  
}


//--- compute U * A ----
SiMath::Vector
SiMath::colProduct(const SiMath::Vector & U, const SiMath::Matrix & A)
{
	if ( U.size() != A.nbrRows() )
		throw(MatrixError("Dimensions of matrix and vector are not compatible."));
	
	Vector v(A.nbrColumns(),0.0);
	
	for ( unsigned int i=0; i<A.nbrColumns(); ++i)
	{
		double s(0.0);
		for ( unsigned int j=0; j<A.nbrRows(); ++j)
		{
			s += U[j] * A[j][i];
		}
		v[i] = s;
	}
	return v;  
}



SiMath::Matrix
SiMath::product(const SiMath::Matrix & A, const SiMath::Matrix & B)
{
	if ( A.nbrColumns() != B.nbrRows() )
		throw(MatrixError("Column dimension of A does not match row dimension of B."));
		
	SiMath::Matrix C(A.nbrRows(),B.nbrColumns(),0);
		
	for ( int i=0; i < A.nbrRows(); ++i )
	{
		for ( int j=0; j < B.nbrColumns();  ++j )
		{
			for ( int k=0; k < A.nbrColumns(); ++k )
			{	
				C[i][j] += A[i][k] * B[k][j];
			}
		}			
	}
	return C;
}


//
void 
SiMath::largeProduct(const SiMath::Matrix & A, const SiMath::Matrix & B, SiMath::Matrix & C)
{
	if ( A.nbrColumns() != B.nbrRows() )
		throw(MatrixError("Column dimension of A does not match row dimension of B."));
		
	// update C
	C.reset(A.nbrRows(),B.nbrColumns());
		
	for ( int i=0; i < A.nbrRows(); ++i )
	{
		for ( int j=0; j < B.nbrColumns();  ++j )
		{
			for ( int k=0; k < A.nbrColumns(); ++k )
			{	
				C[i][j] += A[i][k] * B[k][j];
			}
		}			
	}
	return;	
}

// correlation matrix
SiMath::Matrix
SiMath::correlation(const SiMath::Matrix & A)
{
	// Variables
	unsigned int n(A.nbrColumns());
	unsigned int m(A.nbrRows());

	// correlation matrix [nxn]
	SiMath::Matrix cor(n,n,1.0);
	
	for ( unsigned int k1=0; k1<n-1; ++k1 )
	{
		for ( unsigned int k2=k1+1; k2<n; ++k2 )
		{
			double sumSqrtX = 0.0;
			double sumSqrtY = 0.0;
			double sumXY = 0.0;
			double mX = A[0][k1];
			double mY = A[0][k2];
		
			for (unsigned int i=1; i<n; ++i ){
				double sweep = (i - 1.0) / i;
				double dX = A[i][k1] - mX;
				double dY = A[i][k2] - mY;
				
				sumSqrtX += dX * dX * sweep;
				sumSqrtY += dY * dY * sweep;
				sumXY += dX * dY * sweep;
		
				mX += dY / i;
				mY += dY / i;
			}
	
			// normalise
			sumSqrtX = sqrt( sumSqrtX / n );
			sumSqrtY = sqrt( sumSqrtY / n );
	
			if ( sumSqrtX == 0 || sumSqrtY == 0 ) {
				cor[k1][k2] = 0.0;
				cor[k2][k1] = 0.0;
			}else{
				cor[k1][k2] =  sumXY/(n * sumSqrtX * sumSqrtY );	
				cor[k2][k1] = cor[k1][k2];
			}
		}
	}
	return cor;
}


//--------------------------------------------------------------------------------------------------------------
SiMath::Matrix
SiMath::allPairShortestPath(const SiMath::Matrix & A, SiMath::Matrix & P)
{
	// initialize cost and path matrix
	int n = A.nbrColumns();
	SiMath::Matrix B(A);              // B is a copy of A
	P.reset(n,n); 
	P = -1;     // P is set to all -1
	for ( int k=0; k < n; ++k )
	{
		for ( int i=0; i < n; ++i )
		{
			for ( int j=0; j < n; ++j )
			{
				if ( B[i][k] + B[k][j] < B[i][j] )
				{
					B[i][j] = B[i][k] + B[k][j];
					P[i][j] = k;
				}
			}
		}
	}
	return B;
}


SiMath::Matrix
SiMath::allPairShortestPath(const SiMath::Matrix & A)
{
	SiMath::Matrix B(A); // make a copy of A
	
	// initialize cost and path matrix
	int n = A.nbrColumns();
	double dummy;
	for ( int k=0; k < n; ++k )
	{
		for ( int i=0; i < n; ++i )
		{
			for ( int j=0; j < n; ++j )
			{
				if ( (dummy = B[i][k] + B[k][j]) < B[i][j] )
				{
					B[i][j] = dummy; // update if path is shorter 
				}
			}
		}
	}
	return B;
}



// JACOBI TRANSFORMATION
void 
SiMath::jacobi (Matrix & A, Vector & eigD, Matrix & eigV, unsigned int & nRot)
{
	double onorm, dnorm;
	double b, dma, q, t, c, s;
	double atemp, vtemp, dtemp;
	int i, j, k, l;
	
	unsigned int nRows = A.nbrRows();

	if ( nRows != A.nbrColumns() )
	{
		throw(SiMathError("A not square. Unable to perfrom Jacobi transformation."));
	}
	
	// reset the eigen values and eigen vectors
	eigD.reset(nRows);
	eigV.reset(nRows,nRows);	
	for (i = 0; i < nRows; ++i) {
		eigV[i][i] = 1.0;
		eigD[i] = A[i][i];
	}
	
	// loop for maximal number of iterations
	for (l = 0; l < nRot; ++l) {
		// size of eigen values and off-diagonal elements
		dnorm = 0.0;
		onorm = 0.0;
		for (i = 0; i < nRows-1; ++i) {
			dnorm = dnorm + fabs(eigD[i]);
			for (j = i+1; j < nRows; ++j) {
				onorm += fabs(A[i][j]);
			}
		}

		std::cerr << l << " ";
		// check convergence
		if( (onorm/dnorm) <= 1.0e-12) 
			break;
		
		for ( i = 0; i < nRows - 1; ++i ) {
			for (j = i+1; j < nRows; ++j) {
				b = A[i][j];
				if(fabs(b) > 0.0) {
					dma = eigD[j] - eigD[i];
					if((fabs(dma) + fabs(b)) <=  fabs(dma)) {
						t = b / dma;
					}
					else {
						q = 0.5 * dma / b;
						t = 1.0/(fabs(q) + sqrt(1.0+q*q));
						if(q < 0.0) {
							t = -t;
						}
					}
					c = 1.0/sqrt(t * t + 1.0);
					s = t * c;
					A[i][j] = 0.0;
					for (k = 0; k < i; k++) {
						atemp = c * A[k][i] - s * A[k][j];
						A[k][j] = s * A[k][i] + c * A[k][j];
						A[k][i] = atemp;
					}
					for (k = i+1; k < j; k++) {
						atemp = c * A[i][k] - s * A[k][j];
						A[k][j] = s * A[i][k] + c * A[k][j];
						A[i][k] = atemp;
					}
					for (k = j+1; k < nRows; k++) {
						atemp = c * A[i][k] - s * A[j][k];
						A[j][k] = s * A[i][k] + c * A[j][k];
						A[i][k] = atemp;
					}
					for (k = 0; k < nRows; k++) {
						vtemp = c * eigV[k][i] - s * eigV[k][j];
						eigV[k][j] = s * eigV[k][i] + c * eigV[k][j];
						eigV[k][i] = vtemp;
					}
					dtemp = c*c*eigD[i] + s*s*eigD[j] - 2.0*c*s*b;
					eigD[j] = s*s*eigD[i] + c*c*eigD[j] +  2.0*c*s*b;
					eigD[i] = dtemp;
				}  // end if
			}
		}
	}
	// store number of rotations
	nRot = l;

	std::cerr << " => nbr rotations: " << nRot << std::endl; 
	// sort eigen values
	for ( i = 0; i < nRows - 1; ++i ) {
		std::cerr << " " << i;
		k = i;
		dtemp = eigD[k];
		for ( j = i; j < nRows; ++j ) {
			if( eigD[j] >= dtemp) {
				k = j;
				dtemp = eigD[k];
			}
		}

		// check if order is changed
		if( k != i ) {
			// switch eigenvalues k en i
			eigD[k] = eigD[i];
			eigD[i] = dtemp;
			// switch columns k and i 
			for ( j = 0; j < nRows; j++) {
				dtemp = eigV[j][i];
				eigV[j][i] = eigV[j][k];
				eigV[j][k] = dtemp;
			}
		}
	}
}

