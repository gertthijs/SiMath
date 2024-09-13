/*
 *  Distribution.cpp
 *
 *  Created by Gert Thijs on 17/05/06.
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

#include "Distribution.h"


//--------------------------------------------------------------------------------------
//  NORMAL DISTRIBUTION
//--------------------------------------------------------------------------------------
SiMath::NormalDist::NormalDist(double s, double m, double sh) : 
	_c1(0.0),
	_c2(0.0),
	Distribution(s,m,sh)
{
	if ( s <= 0.0 )
		throw(DistributionError("Scale parameter of normal distribution should be greater than 0."));
		
	_scale = s;
	_offset = m;
	_shape = sh;
			
	_c1 = 1.0/(sqrt(2.0 * PI) * _scale);
	_c2 = 0.5/(_scale * _scale);
}


void 
SiMath::NormalDist::setScale(double v)
{
	if ( v <= 0.0 )
		throw(SiMath::DistributionError("Scale parameter of normal distribution should be greater than 0."));
		
	_scale = v;
	_c1 = 1.0/(sqrt(2.0 * PI) * _scale);
	_c2 = 0.5/(_scale * _scale);
}

double 
SiMath::NormalDist::pdf(double x)
{
	return _c1 * exp( - _c2 * (x - _offset) * (x - _offset));
}


// numerical estimate of cumulative normal distribution
double 
SiMath::NormalDist::cdf(double x)
{
	// numerical approximation 
	// normalise input
	double z((x-_offset)/_scale);
	
	if ( z > 6.0 ) 
		return 1.0;
	if ( z < -6.0 )
		return 0.0;
	
	// constants
	const double B1 = 0.31938153;
	const double B2 = -0.356563782;
	const double B3 = 1.781477937;
	const double B4 = -1.82155978;
	const double B5 = 0.2316419;
	const double P = 0.2316419;
	const double C = 0.3989423;
	
	double t = 1.0/(1.0 + fabs(z) * P);
	double v = ((((B5 * t + B4)*t+B3)*t+B2)*t+B1)*t;
	double b = C * exp( (-z) * (z / 2.0));
	if ( z < 0.0 )
		return b * v;
	else
		return 1.0 - b * v;
}


double 
SiMath::NormalDist::sf(double x)
{
	throw(SiMath::DistributionError("Normal Survival Function not yet implemented."));
}

bool
SiMath::NormalDist::train(SiMath::Vector & data)
{
	_offset = data.mean();
	_scale = data.stDev(_offset);
	_c1 = 1.0/(sqrt(2.0 * PI) * _scale);
	_c2 = 0.5/(_scale * _scale);
	return true;
}


double
SiMath::NormalDist::sample()
{
  double x(0.0), y(0.0), r(0.0);
	
  do
	{
		// choose x,y in uniform from (-1,1)
		x = SiMath::randD(-1,1);
		y = SiMath::randD(-1,1);
		
		// radius 
		r = x * x + y * y;
	}
  while ( r > 1.0 || r == 0 );
	
  // Box-Muller transform
  return _offset + (_scale * y * sqrt (-2.0 * log (r) / r));
}

SiMath::Vector
SiMath::NormalDist::sample(unsigned int n)
{
	SiMath::Vector v(n,0.0);
  double x(0.0), y(0.0), r(0.0);
	
	for ( unsigned int i=0; i<n; ++i )
	{
		do
		{
			// choose x,y in uniform from (-1,1)
			x = SiMath::randD(-1,1);
			y = SiMath::randD(-1,1);
			
			// radius 
			r = x * x + y * y;
		}
		while ( r > 1.0 || r == 0 );
	
		// Box-Muller transform
		v[i] =  _offset + (_scale * y * sqrt (-2.0 * log (r) / r));
	}
	
	return v;
}


//--------------------------------------------------------------------------------------
//  GAMMA DISTRIBUTION
//--------------------------------------------------------------------------------------
SiMath::GammaDist::GammaDist(double s, double o, double sh) :
	Distribution(s,o,sh)
{
		if ( s <= 0.0 )
			throw(SiMath::DistributionError("Scale parameter of gamma distribution should be greater than 0."));
		
		_scale = s;
		
		if ( sh <= 0.0 )
			throw(SiMath::DistributionError("Shape parameter of gamma distribution should be greater than 0."));
		
		_shape = sh;
}

bool 
SiMath::GammaDist::train(SiMath::Vector & data)
{
	if ( data.min() < 0 )
		throw(SiMath::DistributionError("Unable to generate gamma distribution with negative data points."));
	
	double m = data.mean();
	double s = data.stDev(m);
	
	_shape = (m * m)/(s * s);
	_scale = (s * s)/m;
	_offset = 0.0;

	return true;	
}

void  
SiMath::GammaDist::setScale(double v)
{
	if ( v <= 0.0 )
		throw(SiMath::DistributionError("Scale parameter of gamma distribution should be greater than 0."));
		
	_scale = v;
}

void 
SiMath::GammaDist::setShape(double v)
{	
	if ( v <= 0.0 )
		throw(SiMath::DistributionError("Shape parameter of gamma distribution should be greater than 0."));
		
	_shape = v;
}

double
SiMath::GammaDist::pdf(double x)
{
	if ( x < _offset )
		return 0.0;
	double m = (x - _offset)/_scale;
	
	double d = SiMath::lnGamma(_shape) - m - log(_scale);
	d += (_shape - 1)*log(m);	
	return exp(d);
}

double
SiMath::GammaDist::cdf(double x)
{
	if ( x < _offset || x < 0 )
		return 0.0;

	return (SiMath::incompGamma(_shape, (x-_offset)/_scale))/_scale;
}

double
SiMath::GammaDist::sf(double x)
{
	if ( x < _offset )
		return 1.0;

	return 1.0 - (SiMath::incompGamma(_shape, (x-_offset)/_scale))/_scale;
}


//--------------------------------------------------------------------------------------
//  EXTREME VALUE DISTRIBUTION (maximum form)
//--------------------------------------------------------------------------------------
bool
SiMath::MaxExtremeValueDist::train(SiMath::Vector & data)
{
	_shape = 1.0;
	_scale = data.stDev() * PI / sqrt(6);
	_offset = data.mean() - 0.5772 * _scale;
	
	return true;
}


double 
SiMath::MaxExtremeValueDist::pdf(double x)
{
	double d = (_offset - x)/_scale;
	d = exp(d);
	
	return d * exp(-d) / _scale;
}

double 
SiMath::MaxExtremeValueDist::cdf(double x)
{
	double d = (_offset - x)/_scale;
	d = exp(d);
	
	return exp(-d);
}

double 
SiMath::MaxExtremeValueDist::sf(double x)
{
	double d = (_offset - x)/_scale;
	d = exp(d);
	
	return 1 - exp(-d);
}


double
SiMath::MaxExtremeValueDist::sample()
{
	// based on inverse cdf
	return log( -log(SiMath::randD(0,1)) )  * _scale + _offset;
}

SiMath::Vector
SiMath::MaxExtremeValueDist::sample(unsigned int n)
{
	Vector v(n,0.0);

	// based on inverse cdf
	for ( unsigned int i=0; i<n; ++i ){
		v[i] = log( -log(SiMath::randD(0,1)) ) * _scale + _offset;
	}
	
	return v;
}

