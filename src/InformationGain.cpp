/*
 *  InformationGain.cpp
 *
 *  Created by Jonatan Taminau
 *  Copyright 2006, Silicos NV. All rights reserved.
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

#include "InformationGain.h"
#include <iomanip>

//----------------------------------------------------------------------------//
// calculate																																	//
//----------------------------------------------------------------------------//
SiMath::Vector
SiMath::InformationGain::calculate(const SiMath::Matrix & A, const SiMath::Vector & labels)
{
	// dimensions  of A [nxm]
	unsigned int n(A.nbrRows());
	unsigned int m(A.nbrColumns());
			
	if ( n != labels.size() )
		throw(SiMath::InformationGainError("Number of labels not compatible with size of data matrix."));
	
	std::map<double, int> bins;
	std::map<double, int>::iterator itM;
	
	// create initial gain vector
	SiMath::Vector gain(m,_entropy(labels));
	// std::cerr << "GAIN: " << gain;
	
	// loop over all columns
	for ( unsigned int j=0; j<m; ++j )
	{
		// count the number of entries per bin
		bins.clear();
		for( unsigned int i=0 ; i<n ; ++i) 
		{
			itM = bins.find(A[i][j]);
			if(itM==bins.end()) 
			{
				bins.insert(std::make_pair(A[i][j], 1));
			}
			else 
			{
				++itM->second;
			}
		}
		
		for( itM=bins.begin() ; itM!=bins.end() ; ++itM )
		{
			std::map<double, int> classes;
			std::map<double, int>::iterator itC;
			int nbrRows(0);
			// std::cerr << "next bin: " << itM->first << " = " << itM->second << std::endl;
			for(int i=0 ; i<n ; ++i) {
				if( A[i][j]==itM->first )
				{
					itC = classes.find(labels[i]);
					if(itC==classes.end()) {
						classes.insert(std::make_pair(labels[i], 1));
					}
					else {
						++itC->second;
					}
					++nbrRows;
				}
			}
			
			double e(0.0);
			for(itC=classes.begin() ; itC!=classes.end() ; ++itC)
			{
				double p(itC->second / (double)nbrRows);
				e += -p * log2(p);
			}
			double contribution((nbrRows/(double)n)*e);
			gain[j] -= contribution;
			
			//std::cerr << "* entropy of bin " << itM->first << ": " << std::setw(8) << std::left 
			//	<< e << " (scaled contribution: -" << std::setw(8) << std::left 
			//	<< contribution << ")  bin size: " << itM->second << std::endl;
		}
	}
	return gain;
}

//------------------------------------------------------------------------------------------------------------------------
double
SiMath::InformationGain::calculate(const SiMath::Vector & A, const SiMath::Vector & labels)
{
	unsigned int n(A.size());
	
	if ( n != labels.size() )
		throw(SiMath::InformationGainError("Number of labels not compatible with size of data matrix."));
	
	std::map<double, int> bins;
	std::map<double, int>::iterator itM;
	
	// create initial gain vector
	double gain(_entropy(labels));
	// std::cerr << "GAIN: " << gain;
	// count the number of entries per bin
	bins.clear();
	for( unsigned int i=0 ; i<n ; ++i) 
	{
		itM = bins.find(A[i]);
		if(itM==bins.end()) 
		{
			bins.insert(std::make_pair(A[i], 1));
		}
		else 
		{
			++itM->second;
		}
	}
		
	for( itM=bins.begin() ; itM!=bins.end() ; ++itM )
	{
		std::map<double, int> classes;
		std::map<double, int>::iterator itC;
		int nbrRows(0);
		for(int i=0 ; i<n ; ++i) {
			if( A[i] == itM->first )
			{
				itC = classes.find(labels[i]);
				if(itC==classes.end())
				{
					classes.insert(std::make_pair(labels[i], 1));
				}
				else 
				{
					++itC->second;
				}
				++nbrRows;
			}
		}
			
		double e(0.0);
		for(itC=classes.begin() ; itC!=classes.end() ; ++itC)
		{
			double p(itC->second / (double)nbrRows);
			e += -p * log2(p);
		}
		double contribution((nbrRows/(double)n)*e);
		gain -= contribution;
		
		//std::cerr << "* entropy of bin " << itM->first << ": " << std::setw(8) << std::left 
		//			<< e << " (scaled contribution: -" << std::setw(8) << std::left 
		//			<< contribution << ")  bin size: " << itM->second << std::endl;
	}
	return gain;
}


//----------------------------------------------------------------------------//
// _entropy	calculation            																						//
//----------------------------------------------------------------------------//
double 
SiMath::InformationGain::_entropy(const SiMath::Vector & y)
{
	double e(0.0);
	std::map<double, int> classes;
	std::map<double, int>::iterator itM;
	unsigned int r(y.size());
	
	for(int i=0 ; i<r ; ++i) 
	{
		itM = classes.find(y[i]);
		if(itM==classes.end()) 
		{
			classes.insert(std::make_pair(y[i], 1));
		}
		else 
		{
			++itM->second;
		}
	}
	
	for(itM=classes.begin() ; itM!=classes.end() ; ++itM) 
	{
		double p(itM->second / (double)r);
		if( p != 0 ) 
			e += -p * log2(p);
	}
	return e;
}
