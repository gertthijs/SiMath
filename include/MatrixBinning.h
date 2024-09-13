/*
 *  MatrixBinning.h
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


#ifndef __SIMATH_MATRIXBINNING__
#define __SIMATH_MATRIXBINNING__

#include "Matrix.h"
#include "Vector.h"

namespace SiMath
{
	/**
		\class MatrixBinningError
	 \brief  Exception class
	 
	 This exception class is derived from the std::runtime_error class and is used to 
	 indicate that something went with a MatrixBinning object.
	 	\ingroup _error_handling_
	 */	
	class SIMATH_EXCEPTION MatrixBinningError : public SiMathError
	{	
		public:
		MatrixBinningError(const std::string & e) : SiMathError("MatrixBinningError : " + e ) {};
	};

	/**
		\class MatrixBinning
		\brief Class to transform a matrix into a discrete form (either by column or by row)
		
		In the first step of the binning procedure, the boundaries of the bins are computed.
		The discrete bins are equidistant in the range from the minimal to the maximal value.
		For instance, if the min is 0 and the maximal value is 10 and the number of bins is 
		set to 5. The boundaries of the bins are [0,2), [2,4), [4,6), [6,8), and [8,10). A value of 1.3 
		will be assigned to the first bin get bin number 0. The value 8.9 will get bin number 4.
		
		Calling this class on matrix generates a discrete matrix on either the columns (default) 
		or the rows. 
		
		An example of binning a vector is given below (in pseudo code)
		\code 
		v = [0.1 0.5 1.2 4.5 5.1 3.1 2.1 1.1]
		max = 6
		min = 0
		size = 10
		boundaries = [0 0.6 1.2 1.8 2.4 3.0 3.6 4.2 4.8 5.4 6.0]
		bin_v = [0 0 2 7 8 5 3 1]
		\endcode
	 */
	class MatrixBinning 
	{
		public:
			/// \name Structors 
			//@{
			/**
				\brief Default constructor creates empty binning
			 */
			MatrixBinning() : _minValues(), _maxValues(), _binSizes(), _boundaries(), _byRows(false), _bBins(false) {};
		
			/**
				\brief Constructor
				\param m1 Vector with the minimal values for each variable
				\param m2 Vector with the maximal values for each variable
				\param b  Vector with the number of bins for each variable
				\param byRow Defines if the bins should be computed per row or per column (default is column)
			 */
			MatrixBinning(const Vector & m1, const Vector & m2, const Vector & b, bool byRow=false);

			/**
				\brief Destructor
			 */
			~MatrixBinning() {};
			//@}
			
			///\name Bin definition
			//@{
			/**
				\brief Define the maximal value of the bins 
				\param m Vector of maximal value for each variable
			 */
			void setMaximalValues(const Vector & m);

			/**
				\brief Define the minimal value of the bins
				\param m Vector of minimal value  for each variable
			*/
			void setMinimalValues(const Vector & m);

			/**
				\brief Define the sizes of the bins
				\param b Vector of sizes of each bin
			 */
			void setBinSizes(const Vector & b);

			/**
				\brief Create the binning scheme from min, max and number of bins
			 */
			void updateBoundaries();
			
			/**
				\brief Define if the binning is done on the rows or the columns
				\param b If true, each row is binned
			 */
			void setByRows(bool b) {_byRows = b; };
			//@}
			
			/**
				\brief Main method to create the matrix binning
				\param  A Reference to the input matrix
				\return Discretised version of the input matrix
				\exception MatrixBinningError if dimensions of A do not the dimension of the binning scheme
			 */
			Matrix calculate(const Matrix & A);

			///\name Inspectors
			//@{
			/**
				\brief Get the boundaries of column or row i
				\param i Index in the binning scheme
			 */
			Vector getBoundaries(unsigned int i) const;
			
			/**
				\brief Get the current set of maximal values
			 */
			Vector getMaximalValues() const { return _maxValues; };

			/**
				\brief Get the current set of minimal values
			 */
			Vector getMinimalValues() const { return _minValues; };

			/**
				\brief Get the current set of bin sizes
			 */
			Vector getBinSizes() const { return _binSizes; };
			//@}
			
		private:
			///\name Local variables
			//@{
			Vector  _minValues;  ///< Vector to hold the minimal value of the bins
			Vector  _maxValues;  ///< Vector to hold the maximal value of the bins
			Vector  _binSizes;   ///< Vector to hold the number of bins per property
 		
			std::vector<Vector> _boundaries; ///< std::vector to hold the boundaries of the bins per property
			
			bool _byRows;        ///< Indicates if matrix should be binned by row or by column
			bool _bBins;         ///< Status flag to indicate if the binning scheme has been compute or not
			//@}
			
			/// \name Local helper functions
			//@{
			/**
				\brief Local helper function to decide in which bin the value falls
				\param c Bin index (column or row)
				\param v Value to bin
				\return bin value
			 */
			unsigned int _findBin(unsigned int c, double v);
			//@}
	};	
};


#endif __SIMATH_MATRIXBINNING__
