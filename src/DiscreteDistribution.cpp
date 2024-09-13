/*
 *  DiscreteDistribution.cpp
 *  SiMath
 *
 *  Created by Gert Thijs on 9/23/08.
 *  Copyright 2008 Silicos NV.. All rights reserved.
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


#include "DiscreteDistribution.h"

SiMath::BinomialDistribution::BinomialDistribution(double s) :
	_successRate(s)
{
	if ( s < 0.0 || s > 1.0 ){
		_successRate = 0.0;
		throw(DiscreteDistributionError("Success rate of binomial distribution should be in [0,1]."));
	}
}

SiMath::BinomialDistribution::BinomialDistribution(unsigned int p, unsigned int n) :
	_successRate(0.0)
{
	if ( n == 0 || p > n ){
		_successRate = 1.0;
	}else{
		_successRate = (double)(p/n);
	}
}


void 
SiMath::BinomialDistribution::setSuccessRate(double s){
	if ( s < 0.0 || s > 1.0 ){
		_successRate = 0.0;
		throw(DiscreteDistributionError("Success rate of binomial distribution should be in [0,1]."));
	}else{
		_successRate = s;
	}
}

double SiMath::BinomialDistribution::pdf(unsigned int p, unsigned int n){
	double v(0.0);
	if ( p == 0 || n == 0 || p > n ){
		return v;
	}

	v = p * log(_successRate);
	v += (n-p) * log(1 - _successRate);
	
	return  exp(log(SiMath::nchoosek(p,n)) + v);
	
}


SiMath::Vector
SiMath::BinomialDistribution::sample(unsigned int n){
	SiMath::Vector v(n);
	for ( unsigned int i=0; i<n; ++i) {
		if ( SiMath::randD(0,1) < _successRate ){
			v[i] = 1;
		}else{
			v[i] = 0;
		}
	}
	return v;
}


bool
SiMath::BinomialDistribution::sample(){
	if ( SiMath::randD(0,1) < _successRate ){
		return 1;
	}else{
		return 0;
	}
};

