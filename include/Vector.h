/*
 *  Vector.h
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

#ifndef   __SIMATH_VECTOR__
#define   __SIMATH_VECTOR__

#include "Definitions.h"
#include <vector>
#include <iostream>

namespace SiMath
{
	/**
		\class VectorError
	 \brief  Exception class
	 
	 This exception class is derived from the std::runtime_error class and is used to 
	 indicate that something went wrong when manipulating a Vector.
	 \ingroup _error_handling_
	 */	
	class SIMATH_EXCEPTION VectorError : public SiMathError
	{	
		public:
			VectorError(const std::string & e) : SiMathError("VectorError : " + e ) {};
	};

	
	/**
		\class Vector
		\brief Class to hold a vector of double values
	 */
  class Vector
	{
		private:
		//@{
			unsigned int _n;                   ///< Number of data points in vector
			std::vector<double> _pVector;	     ///< std::vector to hold all values
		//@}
	
		public:
		//@{
			Vector() : _n(0), _pVector(0) {};                                    ///< Empty vector
			Vector(const unsigned int n) : _n(n), _pVector(n) {};                ///< vector of length n, no initial value
			Vector(const unsigned int n, const double & v) : _n(n), _pVector(n, v) {};   ///< vector of length n, initial constant value v
			Vector(const unsigned int n, const double * v);                           ///< Copy of data stored in an array double[n]
			Vector(const std::vector<double> & v);                                    ///< Copy of data stored in std::vector<double>
			
			Vector(const Vector & src);               ///< copy constructor

			~Vector();    ///< destructor
			
			/**
				\brief Create an empty vector of size 0.
			 */
			void clear();
			
			/**
				\brief Reset the vector to have size n and all values equal to 0.
				\param n New length of the vector.
			*/
			void reset(unsigned int n);
			 
			/**
				\brief Resize the vector to have length n.
			 
				If n is smaller than the current size the last points will be deleted 
				and the first n points will be the same.
			  If n is larger than the current, new entries with value 0.0 will be 
				added at the end of the vector.
			 
				\param n New length of the vector.
			 */
			void resize(unsigned int n);
		//@}
		
		//@{
			/**
				\brief Save method to get the value at position i
				\param i index
				\return value stored at position i
				\exception VectorError when i is out of range
			 */
			double getValueAt(const unsigned int i);

			/**
				\brief Save method to get the value at position i
			 \param i index
			 \return value stored at position i
			 \exception VectorError when i is out of range
			 */
			double getValueAt(const unsigned int i) const;
			
			/** 
				\brief Get maximal value of the vector
				\return Maximal value
			 */
			double max() const;

			/** 
				\brief Get maximal value of the vector and store the position
				\param index Position of the maximal value in the vector
				\return Maximal value
				*/
			double max(unsigned int & index) const;
			
			/** 
				\brief Get minimal value of the vector
				\return Minimal value
				*/
			double min() const;

			/** 
				\brief Get minimal value of the vector
				
				\return Minimal value
			 */
			double min(unsigned int & index) const; 
			
			/**
				\brief Compute the sum of all values in the vector
			 */
			double sum() const;
			
			/**
				\brief Compute the mean of all values in the vector
			 */
			double mean() const;
			
			/**
				\brief Compute the standard deviation of all values in the vector
			 */
			double stDev() const;
			
			/**
				\brief Compute the standard deviation of all values in the vector given the mean.
				\param m Mean of the vector
			 */
			double stDev(double m) const;
			
			/**
				\brief Get the size of the vector
			 */
			unsigned int size() const { return _n; };
		//@}

		///\group Operators
		//@{
			Vector & operator= (const Vector & src);    ///< copy assignment, resets the size of the Vector if needed
			Vector & operator= (const double & v);      ///< set all elements in Vector to constant value
			Vector & operator+= (const double & v);     ///< add constant value to all elements in Vector
			Vector & operator+= (const Vector & V);     ///< add full Vector element-wise
			Vector & operator-= (const double & v);     ///< subtract constant value to all elements in Vector
			Vector & operator-= (const Vector & V);     ///< subtract full Vector element-wise
			Vector & operator*= (const double & v);     ///< multiply all elements with a constant value
			Vector & operator*= (const Vector & V);     ///< multiply full Vector element-wise 
			Vector & operator/= (const double & v);     ///< divide all elements with a constant value
			Vector & operator/= (const Vector & V);     ///< divide full Vector element-wise 
			Vector & operator- ();                      ///< change sign of all elements in Vector

			Vector operator+ (const Vector & V) const;  ///< operator to write C = A + B
			Vector operator- (const Vector & V) const;  ///< operator to write C = A - B
			Vector operator* (const Vector & V) const;  ///< operator to write C = A * B
			Vector operator/ (const Vector & V) const;  ///< operator to write C = A / B
			
			bool operator== (const Vector & V) const;   ///< check if two vectors are the same, which is only true if all elements are the same
			bool operator!= (const Vector & V) const;   ///< check if two vectors are different
			
			inline double & operator[] (const unsigned int i) { return	_pVector[i];};       ///< set i-th element from vector
			inline double operator[] (const unsigned int i) const { return	_pVector[i];};  ///< get i-th element from vector (const implementation)
		//@}

		//@{
			/**
				\brief Swap elements i and j in the Vector
				\param i 
				\param j
				\return nothing
				\exception VectorError when i or j is out of range
			 */
			void swap(const unsigned int i, const unsigned int j);
			
			double dotProd(const Vector & v);
			
			std::ostream & print(std::ostream & os);      ///< Print vector in a readable format, for debugging purposes

			const double * getArrayPointer() const { return &(_pVector[0]); }; ///< direct access to the data
 			double * getArrayPointer() { return &(_pVector[0]); };             ///< direct access to the data
		//@}
	};
		
	// output generation of a vector
	// std::ostream & operator<< (std::ostream & os, SiMath::Vector & A){ return A.print(os); };
	std::ostream & operator<< (std::ostream & os, SiMath::Vector & A);
	
};


#endif    __SIMATH_VECTOR__
