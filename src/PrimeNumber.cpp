/*
 *  PrimeNumber.cpp
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

#include "PrimeNumber.h"


// constructor
SiMath::PrimeNumber::PrimeNumber() : 
	_primeCounter(2), _primeNumbers(0), _whichPrime(0)
{
	unsigned int i[16] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 59 };
	for ( int j=0; j<16; ++j )
		_primeNumbers.push_back(i[j]);
	
}


// get the i-th prime number
unsigned int 
SiMath::PrimeNumber::get(unsigned int i)
{
	while ( i >= _primeNumbers.size() )
		_addPrimeNumber();
	
	return _primeNumbers[i];
}

unsigned int 
SiMath::PrimeNumber::next()
{
	++_whichPrime;
	if ( _whichPrime >= _primeNumbers.size() )
		_addPrimeNumber();
	
	return _primeNumbers[_whichPrime];
}

unsigned int
SiMath::PrimeNumber::first()
{
	_whichPrime = 0;
	return _primeNumbers[_whichPrime];
}

unsigned int
SiMath::PrimeNumber::last()
{
	return _primeNumbers[_primeNumbers.size() - 1];
}

//------------------------------------------------------------------------------------------------------------------
void 
SiMath::PrimeNumber::_addPrimeNumber()
{
	bool found = false;
	unsigned int n = _primeNumbers[_primeNumbers.size()-1]; // last prime number
	while ( !found )
	{
		n += _primeCounter;  // next number to be searched is of the form 6i+1 or 6i-1 
												 // std::cerr << n << " ";
		
		found = true;
		for (unsigned int j=2; j<_primeNumbers.size(); ++j )
		{
			if ( ( n % _primeNumbers[j] ) == 0 )  // check if remainder of division by prime number is 0
			{
				found = false;
				break;
			}		
		}
				
		_primeCounter = ( _primeCounter == 2 ) ? 4 : 2; // move to next number
	}
	
	// add prime number to list
	_primeNumbers.push_back(n);
	
	return;
}

