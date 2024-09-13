/*
 *  InformationGain.h
 *
 *  Created by Jonatan Taminau
 *  Copyright (c) 2006 - Silicos NV. All Rights Reserved.
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

#ifndef __SIMATH_INFORMATIONGAIN_H__
#define __SIMATH_INFORMATIONGAIN_H__

#include "Matrix.h"

#include <list>
#include <set>
#include <map>


namespace SiMath
{
	/**
		\class  InformationGainError
		\brief  Exception class
	 
		This exception class is derived from the std::runtime_error class and is used to 
		indicate that something went wrong during the information gain algorithm.
	 \ingroup _error_handling_
	 */	
	class SIMATH_EXCEPTION InformationGainError : public SiMathError
	{
		public:
			InformationGainError(const std::string & e) : SiMathError("InformationGainError : " + e) {};
	};

	/**
	  \class InformationGain
	  \brief Class to compute the information gain from a contingency table
	 
		This class is called on a discretised data matrix and a set of class labels
		to identify those columns in a matrix that contain the most information
		to predict or explain the class labels. 
		This class is useful as a preprocessing step where one likes to reduce the dimensionality
		of the problem by feature extraction.
	 
	 */
	class InformationGain
	{
		public:
			
			/**
			  \brief Default constructor
			 */
			InformationGain() {};
		
			/**
		    \brief Main calculation function for multiple variables
			 
				\param data    Matrix in discrete (binned) format
				\param labels  Vector with class labels
				\exception Throws InformationGainError if dimensions of A and labels do not match
				\return Vector with gain for each variable (=column) in the data matrix wrt the class labels.
			 */
			Vector calculate(const Matrix & data, const Vector & labels);
		
			/**
				\brief Main calculation function for one variable
			 
				\param data    Matrix in discrete (binned) format
				\param labels  Vector with class labels
				\exception Throws InformationGainError if dimensions of A and labels do not match
				\return Gain for the given variable wrt the class labels.
			 */
			double calculate(const Vector & data, const Vector & labels);
			
		private:
		
			/**
				\brief Local helper function to compute entropy
			 */
			double _entropy(const Vector & y);
	};
};

#endif __SIMATH_INFORMATIONGAIN_H__
