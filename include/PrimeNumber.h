/*
 *  PrimeNumber.h
 *
 *  Created by Gert Thijs on 27/03/06.
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

#include "Definitions.h"
#include <vector>

#ifndef  __SIMATH_PRIMENUMBER__
#define  __SIMATH_PRIMENUMBER__

namespace SiMath
{
	/**
		\class PrimeNumber
		\brief Class to generate prime numbers. 
			
		The numbers are compute as requested and stored internally as a std::vector of unsigned integers.
		The numbers can be retrieved by index or using an iterator-like function.
	 
		To compute the next prime number, increasing numbers of the form 6n-1 and 6n+1 starting from the last
	  prime found are created. Next, they are checked if they can be divided by the current set of prime numbers.
		If not, the generated number is the new prime number. 
	 */
	class PrimeNumber
	{
		public:
			/**
				\brief Constructor initialises internal array of integers.
			 */
			PrimeNumber();
		
			/**
				\brief Get the first prime number (2) and reset the array iterator to 0
			 */
			unsigned int first();
			
			/**
				\brief Augments the array iterator and returns the next prime number.
			 
				To start from the beginning of the list always call first(). 
				If at the end of the list, this function will generate the next prime number.
				As such this function can called as long as there is no overflow.
			 */
			unsigned int next();
			
			/**
				\brief Get the last prime number from the current list of prime numbers.
			 */
			unsigned int last();

			/** 
				\brief Get the i-th prime number.
				
				If i is outside the current range of computed prime numbers, new numbers will be added 
				untill the i-th element is found.
				\param i Index of the requested prime number.
			 */
			unsigned int get(unsigned int i);
			
		private:
			unsigned int _primeCounter;               ///< Internal counter to keep track of prime number generation
			unsigned int _whichPrime;                 ///< Internal iterator index to access the prime numbers with first(), last() and next() 
			std::vector<unsigned int> _primeNumbers;  ///< Vector of already generated prime numbers
			void _addPrimeNumber();                   ///< Helper function to generate the next prime number
			
	};
};


#endif   __SIMATH_PRIMENUMBER__
