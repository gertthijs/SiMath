/*
 *  Histogram.h
 *
 *  Created by Gert Thijs on 11/9/06.
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

#ifndef __SIMATH_HISTOGRAM__
#define __SIMATH_HISTOGRAM__

#include "Matrix.h"
#include "Vector.h"

namespace SiMath
{
	/**
		\class HistogramError
	 \brief  Exception class
	 
	 This exception class is derived from the std::runtime_error class and is used to 
	 indicate that something went wrong while computing a histogram.
	 	\ingroup _error_handling_
	 */	
	class SIMATH_EXCEPTION HistogramError : public SiMathError
	{	
		public:
			HistogramError(const std::string & e) : SiMathError("HistogramError : " + e ) {};
	};

/**
\class Histogram
 \brief Class to calculate the histogram.
 
*/
class Histogram 
{
	public:
	/// \name Structors 
	//@{
	/**
	\brief Default constructor creates empty bins
	 */
	Histogram() : _minValue(), _maxValue(), _nbrBins(), _boundaries(), _bBins(false) {};
		
	/**
		\brief Constructor
	 \param m1 minimal value
	 \param m2 maximal value
	 \param bins  the number of bins
	 */
	Histogram(double m1, double m2, unsigned int bins);
	
	/**
				\brief Destructor
	 */
	~Histogram() {};
	//@}
	
	///\name Bin definition
	//@{
	/**
	 \brief Define the maximal value of the bins 
	 \param m Vector of maximal value for each variable
	 */
	void setMaximalValue(double m);
	
	/**
				\brief Define the minimal value of the bins
	 \param m Vector of minimal value  for each variable
			*/
	void setMinimalValue(double m);
	
	/**
		\brief Define the size of the bins
		\param n Number of bins
	 */
	void setNumberBins(unsigned int n);
	
	/**
				\brief Create the binning scheme from min, max and number of bins
	 */
	void updateBoundaries();
	
	//@}
	
	/**
		\brief Main method to create a histogram from a vector
		\param  A Reference to the input vector
		\return counts of each  matrix
		\exception MatrixBinningError if dimensions of A do not the dimension of the binning scheme
	 */
	Vector calculate(const Vector & A);
	
	///\name Inspectors
	//@{
	/**
		\brief Get the boundaries used in the histogram
		\return Vector of size number of bins + 1
	 */
	Vector getBoundaries() const {return _boundaries;};
	
	/**
				\brief Get the current set of maximal values
	 */
	double getMaximalValues() const { return _maxValues; };
	
	/**
				\brief Get the current set of minimal values
	 */
	double getMinimalValue() const { return _minValues; };
	
	/**
		\brief Get the current set of bin sizes
	 */
	Vector getBins() const { return _binSizes; };
	//@}
	
private:
		///\name Local variables
		//@{
		double  _minValue;  ///< holds the minimal value of the bins
		double  _maxValue;  ///< holds the maximal value of the bins
	  unsigned int  _nbrBins;  ///< holds the number of bins per property
	
		Vector _boundaries; ///< holds the boundaries of the bins
	
		bool _bBins;         ///< Status flag to indicate if the binning scheme has been compute or not
		//@}
	
	//@}
};	
};


#endif __SIMATH_MATRIXBINNING__
