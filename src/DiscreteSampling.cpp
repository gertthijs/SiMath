/*
 *  DiscreteSampling.cpp
 *  SiMath
 *
 *  Created by Gert Thijs on 02/03/07.
 *  Copyright 2007 Silicos NV. All rights reserved.
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

 
 
#include "DiscreteSampling.h"
#include "Utilities.h"

SiMath::DiscreteSampling::DiscreteSampling(SiMath::Vector & v) : 
	_sum(0.0),
	_size(v.size()),
	_cumSum(_size)
{
	for ( unsigned int i=0; i<_size; ++i ){
		if ( v[i] < 0 )
			throw(SamplingError("Negative values in the input vector. Unable to sample from this vector."));
		_sum += v[i];
		_cumSum[i] = _sum;
	}
}


unsigned int 
SiMath::DiscreteSampling::sample() 
{
	if ( _sum == 0 )
		throw(SamplingError("Unable to sample from empty distribution"));

	unsigned int index(0);
	
	// get random number between 0 and _sum
	double v = SiMath::randD(0, _sum);
	
	while ( _cumSum[index] < v ){ // last element is _sum and should always be larger than v
		++index;
	}
	
	return index;
}


std::vector<unsigned int> 
SiMath::DiscreteSampling::sample(const unsigned int n) 
{
	if ( _sum == 0 )
		throw(SamplingError("Unable to sample from empty distribution"));

	std::vector<unsigned int> indices(n,0);
	
	for ( unsigned int i=0; i<n; ++i )
	{
		// get random number between 0 and _sum
		double v = SiMath::randD(0, _sum);
		
		indices[i] = 0;
		while ( _cumSum[indices[i]] < v ){ // last element is _sum and should always be larger than v
			indices[i] += 1;
		}
	}

	return indices;
}



// sampling without replacement
unsigned int 
SiMath::DiscreteSampling::sampleAndRemove() 
{
	if ( _sum == 0 )
		throw(SamplingError("Unable to sample from empty distribution"));

	unsigned int index(0);
	
	// get random number between 0 and _sum
	double v = SiMath::randD(0, _sum);
	
	while ( _cumSum[index] < v ){ // last element is _sum and should always be larger than v
		++index;
	}
	double w(0.0);
	if ( index != 0 ){
		w = _cumSum[index] - _cumSum[index - 1];
	}else{
		w = _cumSum[0];
	}
	
	for ( unsigned int j=index; j<_cumSum.size(); ++j ){
		_cumSum[j] -= w;
	}
	_sum = _cumSum[_size];
	
	return index;
}


std::vector<unsigned int> 
SiMath::DiscreteSampling::sampleAndRemove(const unsigned int n) 
{
	if ( _sum == 0 )
		throw(SamplingError("Unable to sample from empty distribution"));

	if ( n > _size )
		throw(SamplingError("Number of points is larger than actual vector. Unable to generate so many sample without replacement."));

	// indices storage
	std::vector<unsigned int> indices(n,0);
	
	for ( unsigned int i=0; i<n; ++i )
	{
		if ( _sum == 0 )
			throw(SamplingError("Unable to sample from empty distribution"));

		// get random number between 0 and _sum
		double v = SiMath::randD(0, _sum);
		
		unsigned int index = 0;
		while ( _cumSum[index] < v ){ // last element is _sum and should always be larger than v
			++index;
		}
		
		double w(0.0);
		if ( index != 0 ){
			w = _cumSum[index] - _cumSum[index - 1];
		}else{
			w = _cumSum[0];
		}
		
		for ( unsigned int j=index; j<_cumSum.size(); ++j ){
			_cumSum[j] -= w;
		}
		_sum = _cumSum[_size];

		indices[i] = index;
	}
	
	return indices;
}





