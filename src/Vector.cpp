/*
 *  Vector.cpp
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

#include "Vector.h"

//------------------------------------------------------------------------------------------------------------------
// CONSTRUCTORS
SiMath::Vector::Vector(const unsigned int n, const double * v) :
	_n(n), 
	_pVector(n)
{
	for ( unsigned int i=0; i<_n; ++i)
		_pVector[i] = v[i];
}


SiMath::Vector::Vector(const std::vector<double> & v) :
	_n(v.size()), 
	_pVector(_n)
{
	for ( unsigned int i=0; i<_n; ++i)
		_pVector[i] = v[i];
}


SiMath::Vector::Vector(const Vector & v) : 
	_n(v._n), 
	_pVector(_n)
{
	for ( unsigned int i=0; i<_n; ++i)
	 _pVector[i] = v._pVector[i];
}


SiMath::Vector::~Vector()
{
	_pVector.clear();
}

//-------------------------------------------------------------------------------------------------------------------
// size manipulation
void 
SiMath::Vector::clear()
{
	_pVector.clear();
	_n = 0;
}

void
SiMath::Vector::reset(unsigned int n)
{
	if ( _n !=  n ) // only reset the vector itself if the new size is larger
		_pVector.resize(n);
		
	_n = n;
	for ( unsigned int i=0; i<_n; ++i )
		_pVector[i] = 0;
}


void
SiMath::Vector::resize(unsigned int n)
{
	if ( _n !=  n )
		_pVector.resize(n);
	
	_n = n;
}


//------------------------------------------------------------------------------------------------------------------------
//
std::ostream & 
SiMath::Vector::print(std::ostream & os)
{
		os << _pVector[0];
		for ( unsigned int j=1; j<_n; ++j)
			os << " " << _pVector[j];
		os << std::endl;
		
		return os;
}


//------------------------------------------------------------------------------------------------------------------------
// STATISTICS
double
SiMath::Vector::getValueAt(const unsigned int i)
{
	if ( i >= _n )
		throw(VectorError("Index is larger than size of vector."));
	return _pVector[i];
}

 	
double 
SiMath::Vector::getValueAt(const unsigned int i) const
{
		if ( i >= _n )
			throw(VectorError("Index is larger than size of vector."));
		return _pVector[i];
}

 	
double 
SiMath::Vector::max() const
{
	double d = _pVector[0];
	for ( unsigned int i=1; i<_n; ++i)
	{
		if ( _pVector[i] > d )
		{
			d = _pVector[i];
		}
	}
	return d;
}

double
SiMath::Vector::max(unsigned int & index) const
{
	double d = _pVector[0];
	for ( unsigned int i=1; i<_n; ++i)
	{
		if ( _pVector[i] > d )
		{
			d = _pVector[i];
			index = i;
		}
	}
	return d;
}

double 
SiMath::Vector::min() const
{
	double d = _pVector[0];
	for ( unsigned int i=1; i<_n; ++i)
	{
		if ( _pVector[i] < d )
		{
			d = _pVector[i];
		}
	}
	return d;
}

double
SiMath::Vector::min(unsigned int & index) const 
{
	double d = _pVector[0];
	for ( unsigned int i=1; i<_n; ++i)
	{
		if ( _pVector[i] > d )
		{
			d = _pVector[i];
			index = i;
		}
	}
	return d;
}


double 
SiMath::Vector::sum() const
{
	double m(0.0);
	for ( unsigned int i=0; i<_n; ++i)
		m += _pVector[i];
	return m;
}


double 
SiMath::Vector::mean() const
{
	double m(0.0);
	for ( unsigned int i=0; i<_n; ++i)
		m += _pVector[i];
	return m/_n;
}

double 
SiMath::Vector::stDev() const
{
	double m(0.0);
	for ( unsigned int i=0; i<_n; ++i)
		m += _pVector[i];
		
	double s(0.0);
	for ( unsigned int i=0; i<_n; ++i)
		s += (m-_pVector[i])*(m-_pVector[i]);
		
	return sqrt(s/(_n-1));
}

double
SiMath::Vector::stDev(double m) const
{
	double s(0.0);
	for ( unsigned int i=0; i<_n; ++i)
		s += (m-_pVector[i])*(m-_pVector[i]);
		
	return sqrt(s/(_n-1));
}


//------------------------------------------------------------------------------------------------------------------------
// OPERATORS
SiMath::Vector & 
SiMath::Vector::operator=(const SiMath::Vector & src)
{
	if ( _n != src._n )
	{	
		_n = src._n;
		_pVector.resize(_n);
	}
	for ( unsigned int i=0; i<_n; ++i)
		_pVector[i] = src._pVector[i];
	
	return *this;
}


SiMath::Vector & 
SiMath::Vector::operator= (const double & v)
{
	for ( unsigned int i=0; i<_n; ++i)
		_pVector[i] = v;
	return *this;
}


SiMath::Vector &
SiMath::Vector::operator+= (const double & v)
{
	for ( unsigned int i=0; i<_n; ++i)
		_pVector[i] += v;
	return *this;
}


SiMath::Vector & 
SiMath::Vector::operator+= (const SiMath::Vector & V)
{
	if ( _n != V._n )
		throw(VectorError("Vectors are of different size."));
	
	for ( unsigned int i=0; i<_n; ++i)
		_pVector[i] += V._pVector[i];
	return *this;
}


SiMath::Vector &
SiMath::Vector::operator-= (const double & v)
{
	for ( unsigned int i=0; i<_n; ++i)
		_pVector[i] -= v;
	return *this;
}


SiMath::Vector & 
SiMath::Vector::operator-= (const SiMath::Vector & V)
{
	if ( _n != V._n )
		throw(VectorError("Vectors are of different size."));
	for ( unsigned int i=0; i<_n; ++i)
		_pVector[i] -= V._pVector[i];
	return *this;
}


SiMath::Vector &
SiMath::Vector::operator*= (const double & v)
{
	for ( unsigned int i=0; i<_n; ++i)
		_pVector[i] *= v;
	return *this;
}


SiMath::Vector & 
SiMath::Vector::operator*= (const SiMath::Vector & V)
{
	if ( _n != V._n )
		throw(VectorError("Vectors are of different size."));
	
	for ( unsigned int i=0; i<_n; ++i)
		_pVector[i] *= V._pVector[i];
		
	return *this;
}


SiMath::Vector &
SiMath::Vector::operator/= (const double & v)
{
	for ( unsigned int i=0; i<_n; ++i)
		_pVector[i] /= v;
	return *this;
}


SiMath::Vector &
SiMath::Vector::operator/= (const Vector & V)
{
	if ( _n != V._n )
		throw(VectorError("Vectors are of different size."));
		
	for ( unsigned int i=0; i<_n; ++i)
		_pVector[i] /= V._pVector[i];
	return *this;
}

 	
SiMath::Vector &
SiMath::Vector::operator- ()
{
	for ( unsigned int i=0; i<_n; ++i)
		_pVector[i] = -_pVector[i];
	return *this;
}

  
SiMath::Vector 
SiMath::Vector::operator+ (const Vector & V) const
{
	if ( _n != V._n )
		throw(VectorError("Vectors are of different size."));
	
	SiMath::Vector r(_n);
	for ( unsigned int i=0; i<_n; ++i)
		r[i] = _pVector[i] + V._pVector[i];
	
	return r;
}

SiMath::Vector 
SiMath::Vector::operator- (const Vector & V) const
{
	if ( _n != V._n )
		throw(VectorError("Vectors are of different size."));
	
	SiMath::Vector r(_n);
	for ( unsigned int i=0; i<_n; ++i)
		r[i] = _pVector[i] - V._pVector[i];
	
	return r;
}

SiMath::Vector 
SiMath::Vector::operator* (const Vector & V) const
{
	if ( _n != V._n )
		throw(VectorError("Vectors are of different size."));
	
	SiMath::Vector r(_n);
	for ( unsigned int i=0; i<_n; ++i)
		r[i] = _pVector[i] * V._pVector[i];
	
	return r;
}

SiMath::Vector 
SiMath::Vector::operator/ (const Vector & V) const
{
	if ( _n != V._n )
		throw(VectorError("Vectors are of different size."));
	
	SiMath::Vector r(_n);
	for ( unsigned int i=0; i<_n; ++i)
		r[i] = _pVector[i] / V._pVector[i];
	
	return r;
}

bool  
SiMath::Vector::operator== (const Vector & V) const
{
	if ( _n != V._n )
		throw(VectorError("Vectors are of different size."));
	
	for ( unsigned int i=0; i<_n; ++i)
	{	
		if ( _pVector[i] != V._pVector[i] ) return false;		
	}
	
	return true;	
}

bool  
SiMath::Vector::operator!= (const Vector & V) const
{
	if ( _n != V._n )
		throw(VectorError("Vectors are of different size."));

	for ( unsigned int i=0; i<_n; ++i)
	{	
		if ( _pVector[i] != V._pVector[i] ) return true;
	}
	
	return false;	
}


double 
SiMath::Vector::dotProd(const Vector & v){
	if ( v.size() != _n ){
		throw(VectorError(""));
	}
	
	double d(0.0);
	for (unsigned int i=0; i<_n; ++i )
	{
		d += _pVector[i] * v[i];
	}
	
	return d;
}



//------------------------------------------------------------------------------------
void 
SiMath::Vector::swap(const unsigned int i, const unsigned int j)
{
		if ( i >= _n || j >= _n)
			throw(VectorError("Indices to swap are larger than the vector size."));
		
		double dummy = _pVector[i];
		_pVector[i] = _pVector[j];
		_pVector[j] = dummy;
		
		return;
}



std::ostream & 
SiMath::operator<< (std::ostream & os, SiMath::Vector & A)
{
	return A.print(os); 
}
