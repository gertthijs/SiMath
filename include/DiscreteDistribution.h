/*
 *  DiscreteDistribution.h
 *  SiMath
 *
 *  Created by Gert Thijs on 9/23/08.
 *  Copyright 2008 Silicos NV.. All rights reserved.
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


#ifndef __SIMATH_DISCRETE_DISTRIBUTION__
#define __SIMATH_DISCRETE_DISTRIBUTION__

#include "Vector.h"
#include "Utilities.h"

namespace SiMath
{
	/**
	 \class DistributionError
	 \brief  Exception class
	 
	 This exception class is derived from the std::runtime_error class and is used to 
	 indicate that something went with a Distribution object.
	 \ingroup _error_handling_
	 */	
	class SIMATH_EXCEPTION DiscreteDistributionError : public SiMathError
	{	
	public:
		DiscreteDistributionError(const std::string & e) : SiMathError("DiscreteDistributionError : " + e ) {};
	};
	
	
	/**
	 \class BinomialDistribution
	 \brief Class for the representation of a binomial discrete distribution
	
	 */
	class BinomialDistribution
	{
		private:	
			double _successRate;
		
		public:
			BinomialDistribution(double rate);
			
			BinomialDistribution(unsigned int success, unsigned int failures);

			~BinomialDistribution(){ };
		
		double getSuccessRate(){ return _successRate; };
		
		void setSuccessRate(double s);
		
		double pdf(unsigned int p, unsigned int n);
		
		//double cdf(unsigned int success, unsigned int attempts);
		
		//double sdf(unsigned int success, unsigned int attempts);
		
		Vector sample(unsigned int n);
		
		bool sample();
	};
	
}


#endif __SIMATH_DISCRETE_DISTRIBUTION__
