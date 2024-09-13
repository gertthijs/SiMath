/*
 *  Utilities.h
 *
 *  Created by Gert Thijs on 27/03/06.
 *  Copyright 2006 Silicos. All rights reserved.
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


#ifndef  __SIMATH_UTILITIES__
#define  __SIMATH_UTILITIES__

#include "Definitions.h"
#include "Vector.h"

namespace SiMath
{
	/**
		\brief Get random double in interval [a,b)
		\param a lower bound
		\param b upper bound
		\return double in the interval [a,b)
	 */
	double randD(double a, double b);
	
	/**
		\brief Get random integer in interval [a,b)
		\param a lower bound
		\param b upper bound
		\return integer in the interval [a,b)
	 */
	int randI(int a, int b);
	
	/**
		\brief Numerically safe way to compute sqrt(a^2 + b^2) 
		\param a
		\param b
	 */
	double triangle(const double & a, const double & b );							
	
	/**
		\brief Function to compute the (positive) integer power of a real number (x^n) 
		
		\param x  input value
	  \param n  power
		\return double
	 */
	double powi(double x, unsigned int n);

	/**
		\brief Function to compute the logarithm of the gamma function \\Gamma(x)
	 
	  \param x input value
	  \return double
	 */
	double lnGamma(double x);
	
	/**
		\brief Function to compute the incomplete Gamma function \\Gamma(a,x)
	 */
	double incompGamma(double a, double x);
	
	/**
		\brief Function to compute the beta function \\Beta(x,y)
		 
		\param x 
		\param y
		\return double
	 */
	double beta(double x, double y);
	
	/**
		\brief Function to compute the incomplete beta function \\Beta_x(a,b)
		 
		\param a
		\param b
		\param x
		\return double
		\exception SiMathError() if x is not within [0,1]
	 */
	double incompBeta(double a, double b, double x);
	
	/**
		\brief Helper function to compute the incomplete beta function
		 
		\param a
		\param b
		\param x
		\return double
		\exception SiMathError() when unable to complete computation within defined precision
	 */
	double cfBeta(double a, double b, double x);
	
  /**
		\brief Function to compute the F probability distribution
		 
	  \param x    Value at which the f probability distribution should be computed   
	  \param df1  Degrees of freedom in nominator
		\param df2  Degrees of freedom in denominator
		\return double
		\exception SiMathError()
	 */
	double FDist(double x , double df1, double df2);

	/**
		\brief Function to compute the sigmoid 
		
		\param x  input value
		\param A
		\param B
		\return double 
	 */
	double sigmoid(double x, double A, double B);

	/**
		\brief Function to fit a sigmoidal function through the given data points (x,y)
	 
	 \param values  Vector of input data 
	 \param targets  Vector of output data
	 \param A    reference to double which holds the first parameter
	 \param B    reference to double which holds the second parameter
	 \return nothing
	 \exception  SiMathError() if the size of input and output do not match
	 */	
	void fitSigmoid( const SiMath::Vector & values, const SiMath::Vector & targets, double & A, double & B);
	
	
	/**
		\brief Function to compute the number of combinations k from n.
	 */
	unsigned int nchoosek(unsigned int k, unsigned int n);
	
};


#endif   __SIMATH_UTILITIES__
