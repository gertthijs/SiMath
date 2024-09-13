/*
 *  DiscreteSampling.h
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


#ifndef __SIMATH_DISCRETESAMPLING__
#define __SIMATH_DISCRETESAMPLING__

#include "Vector.h"

namespace SiMath{
	/**
		\class  SamplingError
		\brief  Exception class
	 
	 This exception class is derived from the std::runtime_error class and is used to 
	 indicate that something went wrong during the discrete sampling procedure.
	 \ingroup _error_handling_
	 */	
	class SIMATH_EXCEPTION SamplingError : public SiMathError
	{
		public:
			SamplingError(const std::string & e) : SiMathError("SamplingError : " + e) {};
	};


	/**
		\class DiscreteSampling
		\brief Class to draw samples from a set of points. 
	 
		The probability by which a point is selected is proportional to its value in the input vector. 
		The input vector should contain real positive values. 
		The sampling can be performed with replacement or without replacement. When used with 
		replacement each generated sample is independent of each each other and the internal state of
		the input vector is not changed. When sampling without replacement is performed the probability
		of a sample depends on the previous sampling steps. 
		
	 */
	class DiscreteSampling {
		public:
			/**
				\brief Constructor
				\param v Reference to the input vector of real positive values. Each value in the input vector represents a 
				\exception SamplingError if there are negative values in the input vector 
			 */
			DiscreteSampling(SiMath::Vector & v);
		
			/**
				\brief Destructor
			 */
			~DiscreteSampling() {};
			
			/**
				\brief Sample with replacement a point from the input vector 
			 \return The index of the sampled point.
				\exception SamplingError if the sum of the input vector is zero 
			 */
			unsigned int sample();
			
			/**
				\brief Sample with replacement n points from the input vector
				\param n The number of points to be sampled.
				\return A std::vector with the indices of the sampled points.
				\exception SamplingError if the sum of the input vector is zero 
			 */
			std::vector<unsigned int> sample(const unsigned int n);

			/**
				\brief Sample without replacement one point from the input vector
				\return The index of the sampled point.
				\exception SamplingError if the sum of the input vector is zero 
			 */
			unsigned int sampleAndRemove();
			
			/**
				\brief Sample n points without replacement from the input vector
				\param n The number of points to be sampled. This should be smaller than the size of the input vector.
				\return A std::vector with the indices of the sampled points.
				\exception SamplingError if the number of points asked for is larger than the size of the input vector
				\exception SamplingError if the sum of the input vector is zero 
			 */
			std::vector<unsigned int> sampleAndRemove(const unsigned int n);
			
			
		private:
			double _sum;
			unsigned int _size;
			std::vector<double> _cumSum;
	};

};

#endif __SIMATH_DISCRETESAMPLING__
