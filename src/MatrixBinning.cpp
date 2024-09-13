/*
 *  MatrixBinning.cpp
 *
 *  Created by Gert Thijs on 14/06/06.
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

#include "MatrixBinning.h"

SiMath::MatrixBinning::MatrixBinning(const Vector & m1, const Vector & m2, const Vector & b, bool byRows) :
	_minValues(m1),
	_maxValues(m2),
	_binSizes(b),
	_boundaries(),
	_byRows(byRows),
	_bBins(false)
{
	if ( m1.size() != m2.size() || m1.size() != b.size() )
	{
		throw(SiMath::MatrixBinningError("dimenion mismatch between min, max and sizes of bins."));
	}
		
	_boundaries.resize(m1.size());
		
	for (unsigned int i=0; i<m1.size(); ++i )
	{
		if ( b[i] < 2.0 )
			throw(SiMath::MatrixBinningError("number of bins smaller than 2."));
		
		double step = (m2[i] - m1[i])/b[i];
		Vector bounds((unsigned int)(b[i]+1));
		bounds[0] = m1[i];
		for ( unsigned int j=1; j<=b[i]; ++j )
			bounds[j] = bounds[j-1] + step;
		
		_boundaries[i] = bounds;
	}
	// bins have been set
	_bBins = true;
}

void
SiMath::MatrixBinning::setMaximalValues(const Vector & m)
{
	_maxValues = m;
	_bBins = false; // bins no longer valid
}

void
SiMath::MatrixBinning::setMinimalValues(const Vector & m)
{
	_minValues = m;
	_bBins = false; // bins no longer valid
}

void
SiMath::MatrixBinning::setBinSizes(const Vector & b)
{
	_binSizes = b;
	_bBins = false; // bins no longer valid
}

void
SiMath::MatrixBinning::updateBoundaries()
{
	if ( _minValues.size() != _maxValues.size() || _minValues.size() != _binSizes.size() )
	{
		throw(SiMath::MatrixBinningError("dimenion mismatch between min, max and sizes of bins."));
	}
	
	if ( _boundaries.size() != _minValues.size() )
		_boundaries.resize(_minValues.size());
		
	for (unsigned int i=0; i<_minValues.size(); ++i )
	{
		if ( _binSizes[i] < 2.0 )
			throw(SiMath::MatrixBinningError("number of bins smaller than 2."));
		
		double step = (_maxValues[i] - _minValues[i])/_binSizes[i];
		SiMath::Vector bounds((unsigned int)(_binSizes[i]+1));
		bounds[0] = _minValues[i];
		for ( unsigned int j=1; j<=_binSizes[i]; ++j )
			bounds[j] = bounds[j-1] + step;
		
		_boundaries[i] = bounds;
	}
	
	// bins have been set
	_bBins = true;
	
	return;
}

SiMath::Matrix
SiMath::MatrixBinning::calculate(const Matrix & A)
{
	if ( !_bBins )
		updateBoundaries();
	
	if ( _byRows )
	{
		if ( A.nbrRows() != _minValues.size() )
			throw(MatrixBinningError("Number of rows incompatible with the binning scheme."));
		
		// make a copy of A
		SiMath::Matrix B(A);
		for ( unsigned j=0; j<A.nbrRows(); ++j )
		{
			for ( unsigned i=0; i<A.nbrColumns(); ++i )
				B[j][i] = _findBin(j,A[j][i]);
		}		
		return B;
	}
	else
	{
		if ( A.nbrColumns() != _minValues.size() )
			throw(MatrixBinningError("Number of columns incompatible with the binning scheme."));
		
		// make a copy of A
		SiMath::Matrix B(A);
		for ( unsigned i=0; i<A.nbrColumns(); ++i )
		{
			for ( unsigned j=0; j<A.nbrRows(); ++j )
				B[j][i] = _findBin(i,A[j][i]);
		}
		return B;
		
	}
}


unsigned int 
SiMath::MatrixBinning::_findBin(unsigned int c, double v)
{
	if ( c > _boundaries.size() )
		throw(MatrixBinningError("Property index out of range"));
	
	unsigned int b(0);
	if ( v < _boundaries[c][0] )
		return 0;
	
	unsigned int e(_boundaries[c].size());
	if ( v >= _boundaries[c][e-1] )
		return e-1;
	
	unsigned int i(0);
	while ( e-b != 1 )
	{
		i = b + ( e - b )/2;
		if ( v < _boundaries[c][i] )
		{
			e = i;
		}
		else
		{
			b = i;
		}
	}
	
	return b;
}

