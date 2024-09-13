/*
 *  Matrix.cpp
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

#include "Matrix.h"
#include "Utilities.h"


//------------------------------------------------------------------------------------------------------------------------
// CONSTRUCTORS
//------------------------------------------------------------------------------------------------------------------------
SiMath::Matrix::Matrix(const unsigned int n, const unsigned int m) :
	_nRows(n),
	_nCols(m),
	_pMatrix(0)
{
	if ( n && m )
	{
		double* dummy = new double[n*m];  // data
		_pMatrix = new double*[n];             // row pointers
		for ( unsigned int i=0; i<n; ++i)
		{
			_pMatrix[i] = dummy;
			dummy += m;
		}
	}
}


SiMath::Matrix::Matrix(const unsigned int n, const unsigned int m, const double & v) :
	_nRows(n),
	_nCols(m),
	_pMatrix(0)
{
	if ( n && m )
	{
		double * dummy = new double[n*m];
		_pMatrix = new double*[n];
		for ( unsigned int i=0; i<n; ++i)
		{
			_pMatrix[i] = dummy;
			dummy += m;
		}
	
		for ( unsigned int i=0; i<n; ++i)
			for ( unsigned int j=0; j<m; ++j )
				_pMatrix[i][j] = v; 
	}
}

SiMath::Matrix::Matrix(const unsigned int n, const unsigned int m, const SiMath::Vector & vec) :
	_nRows(n),
	_nCols(m),
	_pMatrix(0)
{
	if(vec.size() != n*m) {
		throw MatrixError("incorrect size of the input vector");
	}
	
	double * dummy(new double[n*m]);
	_pMatrix = new double*[n];
	for(unsigned int i=0 ; i<n ; ++i) {
		_pMatrix[i] = dummy;
		dummy += m;
	}
	for(unsigned int i=0 ; i<n ; ++i) {
		for(unsigned int j=0 ; j<m ; ++j) {
			_pMatrix[i][j] = vec[i*m+j];
		}
	}
}

SiMath::Matrix::Matrix(const SiMath::Matrix & src) : 
	_nRows(src._nRows),
	_nCols(src._nCols),
	_pMatrix(0)
{
	if ( _nRows && _nCols )
	{
		double * dummy(new double[_nRows * _nCols]);
		_pMatrix = new double*[_nRows];
		for ( unsigned int i=0; i<_nRows; ++i)
		{
			_pMatrix[i] = dummy;
			dummy += _nCols;
		}
	
		for ( unsigned int i=0; i<_nRows; ++i)
			for ( unsigned int j=0; j<_nCols; ++j )
				_pMatrix[i][j] = src[i][j]; 
	}
}


//---------------------------------------------------------------------------------------------------------------
SiMath::Matrix::~Matrix()
{
	// std::cerr << "deleting matrix ..." << _pMatrix << std::endl;
	if ( _pMatrix != NULL )
	{
		if ( _pMatrix[0] != NULL )
			delete[] (_pMatrix[0]);
		delete[] (_pMatrix);
	}
	_pMatrix = NULL;
}


//------------------------------------------------------------------------------------------------------------------------
// ACCESS
double
SiMath::Matrix::getValueAt(const unsigned int i, const unsigned int j)
{
	if ( i < 0 || j < 0 || i >= _nRows || j >= _nRows )
		throw(MatrixError("Indices i or j are out of range."));
		
	return _pMatrix[i][j];
}


const double 
SiMath::Matrix::getValueAt(const unsigned int i, const unsigned int j) const
{
	if ( i < 0 || j < 0 || i >= _nRows || j >= _nRows )
		throw(MatrixError("Indices i or j are out of range."));
		
	return _pMatrix[i][j];
}

SiMath::Vector
SiMath::Matrix::getRow(const unsigned int i) const
{
	if ( i < 0 || i >= _nRows )
		throw(MatrixError("Row index i out of range."));
	
	Vector v(_nCols); 
	for ( unsigned int j=0; j<_nCols; ++j )
		v[j] = _pMatrix[i][j];
	
	return v;
}


SiMath::Vector
SiMath::Matrix::getColumn(const unsigned int i) const
{
	if ( i < 0 || i >= _nCols )
		throw(MatrixError("Column index i out of range."));
	
	Vector v(_nRows);
	for ( unsigned int j=0; j<_nRows; ++j )
		v[j] = _pMatrix[j][i];
	
	return v;
}


//------------------------------------------------------------------------------------------------------------------------
inline void
SiMath::Matrix::setValueAt(const unsigned int i, const unsigned int j, double v)
{
	if ( i < 0 || j < 0 || i >= _nRows || j >= _nRows )
		throw(MatrixError("Indices i or j are out of range."));

	_pMatrix[i][j] = v;
}


void 
SiMath::Matrix::setRow(const unsigned int i, SiMath::Vector & src)
{
	if ( src.size() != _nCols )
		throw(MatrixError("Dimension mismatch between matrix and vector"));
	
	if ( i < 0 || i >= _nRows )
		throw(MatrixError("Row index i out of range."));
	
	for ( unsigned int j=0; j<_nCols; ++j )
		_pMatrix[i][j] = src[j];
}



void 
SiMath::Matrix::setColumn(const unsigned int i, SiMath::Vector & src)
{
	if ( src.size() != _nRows )
		throw(MatrixError("Dimension mismatch between matrix and vector"));
	
	if ( i < 0 || i >= _nCols )
		throw(MatrixError("Column index i out of range."));
	
	for ( unsigned int j=0; j<_nRows; ++j )
		_pMatrix[j][i] = src[j];
}


//------------------------------------------------------------------------------------------------------------------------
// OPERATORS
//
//double * 
//SiMath::Matrix::operator[] (const unsigned int i)
//{
//	return _pMatrix[i];
//}

//const double * 
//SiMath::Matrix::operator[] (const unsigned int i) const
//{
//	return _pMatrix[i];
//}


//------------------------------------------------------------------------------------------------------------------------
SiMath::Matrix &
SiMath::Matrix::operator= (const SiMath::Matrix & M)
{
	// check dimensions
	if ( _nRows != M.nbrRows() ||	_nCols != M.nbrColumns() )
	{
		// throw(std::runtime_error("wrong dimensions"));
		
		if ( _nRows && _pMatrix != 0 )
		{
			// delete old matrix
			if ( _nCols && _pMatrix[0] != NULL )
				delete[] _pMatrix[0];
			delete[] _pMatrix;
		}
		_pMatrix = NULL;
			
		// create a new matrix
		_nRows = M.nbrRows();
		_nCols = M.nbrColumns();
		_pMatrix = new double*[_nRows];
		_pMatrix[0] = new double[_nRows*_nCols];
		for (unsigned int i=1; i<_nRows; ++i)
			_pMatrix[i] = _pMatrix[i-1] + _nCols;
		
	}
	
	// fill in all new values	
	for (unsigned int i=0; i<_nRows; ++i)
		for (unsigned int j=0; j<_nCols; ++j)
			_pMatrix[i][j] = M[i][j];
	
	return *this;
}


//------------------------------------------------------------------------------------------------------------------------
///\brief Assigns constant value to all elements in matrix

SiMath::Matrix & 
SiMath::Matrix::operator= (const double & v) 
{
	for ( unsigned int i=0; i<_nRows; ++i)
		for ( unsigned int j=0; j<_nCols; ++j)
			_pMatrix[i][j] = v;
	
	return *this;
}


//------------------------------------------------------------------------------------------------------------------------
///\brief Add constant value to all elements in matrix

SiMath::Matrix & 
SiMath::Matrix::operator+= (const double & v)
{
	for ( int i=0; i<_nRows; i++)
				for ( int j=0; j<_nCols; j++)
					_pMatrix[i][j] += v;
	return *this;
}

//------------------------------------------------------------------------------------------------------------------------
///\brief Add all elements in matrix	

SiMath::Matrix & 
SiMath::Matrix::operator+= (const SiMath::Matrix & M)
{
	if ( M.nbrRows() != _nRows || M.nbrColumns() != _nCols )
		throw(MatrixError("Dimensions of two matrices do not match while adding them."));
	
	for ( unsigned int i=0; i<_nRows; ++i)
		for ( unsigned int j=0; j<_nCols; ++j)
			_pMatrix[i][j] += M[i][j];
	return *this;
}

//------------------------------------------------------------------------------------------------------------------------
///\brief Subtractconstant value to all elements in matrix

SiMath::Matrix & 
SiMath::Matrix::operator-= (const double & v)
{
	for ( int i=0; i<_nRows; i++)
		for ( int j=0; j<_nCols; j++)
			_pMatrix[i][j] -= v;
	return *this;
}

//------------------------------------------------------------------------------------------------------------------------
///\brief Subtract all elements in matrix	

SiMath::Matrix & 
SiMath::Matrix::operator-= (const SiMath::Matrix & M)
{
	if ( M.nbrRows() != _nRows || M.nbrColumns() != _nCols )
		throw(MatrixError("Dimensions of two matrices do not match while adding them."));
	
	for ( unsigned int i=0; i<_nRows; ++i)
		for ( unsigned int j=0; j<_nCols; ++j)
			_pMatrix[i][j] -= M[i][j];
	return *this;
}

//------------------------------------------------------------------------------------------------------------------------
///\brief Multiply all elements in matrix with a constant value

SiMath::Matrix & 
SiMath::Matrix::operator*= (const double & v)
{
	for ( unsigned int i=0; i<_nRows; ++i )
				for ( unsigned int j=0; j<_nCols; ++j )
					_pMatrix[i][j] *= v;
	return *this;
}


//------------------------------------------------------------------------------------------------------------------------

SiMath::Matrix & 
SiMath::Matrix::operator*= (const SiMath::Matrix & M)
{
	if ( M.nbrRows() != _nRows || M.nbrColumns() != _nCols )
		throw(MatrixError("Dimensions of two matrices do not match while multiplying them element-wise."));
	
	for ( unsigned int i=0; i<_nRows; ++i)
		for ( unsigned int j=0; j<_nCols; ++j)
			_pMatrix[i][j] *= M[i][j];
	return *this;
}

//------------------------------------------------------------------------------------------------------------------------
///\brief Multiply all elements in matrix with a constant value

SiMath::Matrix & 
SiMath::Matrix::operator/= (const double & v)
{
	for ( unsigned int i=0; i<_nRows; ++i )
		for ( unsigned int j=0; j<_nCols; ++j )
			_pMatrix[i][j] /= v;
	return *this;
}


//------------------------------------------------------------------------------------------------------------------------

SiMath::Matrix & 
SiMath::Matrix::operator/= (const SiMath::Matrix & M)
{
	if ( M.nbrRows() != _nRows || M.nbrColumns() != _nCols )
		throw(MatrixError("Dimensions of two matrices do not match while multiplying them element-wise."));
	
	for ( unsigned int i=0; i<_nRows; ++i)
		for ( unsigned int j=0; j<_nCols; ++j)
			_pMatrix[i][j] /= M[i][j];
	return *this;
}

//------------------------------------------------------------------------------------------------------------------------

SiMath::Matrix & 
SiMath::Matrix::operator- ()
{
	for ( unsigned int i=0; i<_nRows; ++i)
		for ( unsigned int j=0; j<_nCols; ++j)
			_pMatrix[i][j] = -_pMatrix[i][j];
	
	return *this;
}

//------------------------------------------------------------------------------------------------------------------------
SiMath::Matrix
SiMath::Matrix::operator+ (const SiMath::Matrix & M) const
{
	SiMath::Matrix B(M);
	for ( unsigned int i=0; i<_nRows; ++i)
		for ( unsigned int j=0; j<_nCols; ++j)
			B[i][j] = _pMatrix[i][j] + M[i][j];
	
	return B;
}

SiMath::Matrix
SiMath::Matrix::operator- (const SiMath::Matrix & M) const
{
	SiMath::Matrix B(M);
	for ( unsigned int i=0; i<_nRows; ++i)
		for ( unsigned int j=0; j<_nCols; ++j)
			B[i][j] = _pMatrix[i][j] - M[i][j];
	
	return B;
}

SiMath::Matrix
SiMath::Matrix::operator* (const SiMath::Matrix & M) const
{
	SiMath::Matrix B(M);
	for ( unsigned int i=0; i<_nRows; ++i)
		for ( unsigned int j=0; j<_nCols; ++j)
			B[i][j] = _pMatrix[i][j] * M[i][j];
	
	return B;
}

SiMath::Matrix
SiMath::Matrix::operator/ (const SiMath::Matrix & M) const
{
	SiMath::Matrix B(M);
	for ( unsigned int i=0; i<_nRows; ++i)
		for ( unsigned int j=0; j<_nCols; ++j)
			B[i][j] = _pMatrix[i][j] / M[i][j];
	
	return B;
}



//------------------------------------------------------------------------------------------------------------------------

std::ostream & 
SiMath::Matrix::print(std::ostream & os) const 
{
		for ( unsigned int i=0; i<_nRows; ++i)
		{
			os << "(";
			for ( unsigned int j=0; j<_nCols; ++j)
				os << " " << _pMatrix[i][j];
			os << " )" << std::endl;
		}
		os << std::endl;
		return os;
}


//------------------------------------------------------------------------------------------------------------------------
void 
SiMath::Matrix::swapRows(unsigned int i, unsigned int j)
{
	double dummy;
	for (unsigned int k=0; k<_nCols; ++k)         // loop over all columns
	{
		dummy = _pMatrix[i][k];            // store original element at [i,k]
		_pMatrix[i][k] = _pMatrix[j][k];   // replace [i,k] with [j,k]
		_pMatrix[j][k] = dummy;            // replace [j,k] with element originally at [i,k] 
	}
	return;
}


//------------------------------------------------------------------------------------------------------------------------
void 
SiMath::Matrix::swapColumns(unsigned int i, unsigned int j)
{
	double dummy;
	for (unsigned int k=0; k<_nRows; ++k)         // loop over all rows
	{
		dummy = _pMatrix[k][i];            // store original element at [k,i]
		_pMatrix[k][i] = _pMatrix[k][j];   // replace [k,i] with [k,j]
		_pMatrix[k][j] = dummy;            // replace [k,j] with element orignally at [k,i]
	}
	return;
}


//------------------------------------------------------------------------------------------------------------------------
void
SiMath::Matrix::reset (const unsigned int r, const unsigned int c)
{
	// check dimensions
	if ( _nRows != r ||	_nCols != c )
	{
		if ( _nRows != 0 && _nCols != 0 && _pMatrix != 0 )
		{
			// delete old matrix
			if ( _pMatrix[0] != NULL )
				delete[] _pMatrix[0];
			delete[] _pMatrix;
		}
				
		// create a new matrix
		_nRows = r;
		_nCols = c;
		if ( _nRows == 0 || _nCols == 0 )
		{
			_pMatrix = NULL;
			return;
		}
			
		_pMatrix = new double*[_nRows];
		_pMatrix[0] = new double[_nRows*_nCols];
		for ( unsigned int i=1; i<_nRows; ++i)
			_pMatrix[i] = _pMatrix[i-1] + _nCols;
	}
	
	// fill in all new values	
	for ( unsigned int i=0; i<_nRows; ++i)
		for ( unsigned int j=0; j<_nCols; ++j)
			_pMatrix[i][j] = 0;
	
}


void 
SiMath::Matrix::clear()
{
	// delete old matrix
	if ( _pMatrix != NULL )
	{
		if ( _pMatrix[0] != NULL )
			delete[] _pMatrix[0];
		delete[] _pMatrix;
	}
	_pMatrix = NULL;
	_nRows = 0;
	_nCols = 0;
}


std::ostream & 
SiMath::operator<< (std::ostream & os, const SiMath::Matrix & A) 
{
	return A.print(os);
}
