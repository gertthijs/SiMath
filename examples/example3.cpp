/*
 *  example3.cpp
 *  SiMath
 *
 *  Created by Gert Thijs on 10/25/06.
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


/**
	\file example3.cpp
	\brief Procedure to select the most diverse set of n points from a data matrix
 */

#include <time.h>
#include <fstream>
#include <string>
#include <vector>
#include <iostream>

#include "SiMath/SiMath.h"

void splitLine(const std::string & l, std::vector<double> & ls);

int
main(int argc, char * const argv[]) 
{
	//initialise random number generator
	srandom(time(NULL));
	
	SiMath::NormalDist dd(5,2);
	
	std::cerr << dd.sample();
	for ( unsigned int i=0; i<1000; ++i )
		std::cerr << "," << dd.sample();
	std::cerr << std::endl;
	
	// open tab-delimited text file for reading input data
	// in this example the first column contains the target values
	std::string inFile(argv[1]);
	unsigned int nCompounds = strtol(argv[2],NULL,10);
	std::ifstream ifs(inFile.c_str());
	std::string line;
	unsigned int c(0); // line counter
	std::vector<double> values;
	while( !ifs.eof() )
	{
		// read next line
		getline(ifs, line);
		
		if ( line.empty() ) continue;
		c++;
		
		// store elements in local vector
		splitLine(line, values);
	}
	ifs.close();
	ifs.clear();
	
	// c is the number of lines (==rows) in the data matrix
	// number of columns is values.size / c
	if ( (values.size() % c) != 0 )
	{
		std::cerr << "ERROR: Number of data points is not consistent." << std::endl;
		exit(101);
	}
	unsigned int n = values.size() / c; // number of variables

	// create data matrix
	SiMath::Matrix x(c,n);
	for ( unsigned int i=0; i<c; ++i )
	{
		for ( unsigned int j=0; j<n; ++j )
			x[i][j] = values[i*n + j];
	}
	values.clear();

	// print correlation matrix 
	SiMath::Matrix CC = SiMath::correlation(x);
	std::cerr << CC;
		
	// use mx to store current center
	SiMath::Vector mx = meanOfAllColumns(x);
	std::vector<bool> seen(c,false);
	
	// find first data point most distant to center
	double d(0.0);
	double dFar(0.0);
	unsigned int index(0);
	for ( unsigned int k=0; k<c; ++k )
	{
		if ( seen[k] ) continue;
		
		d = SiMath::similarityByEuclidean(n,x[k], mx.getArrayPointer());
		if ( d > dFar ){
			dFar = d;
			index = k;
		}
	}
	// store current 	
	mx = x.getRow(index);
	std::cout << index << ": " << mx;
	seen[index] = true;
	
	// create vector to hold the distance of each compound to the closest compound in the current selection
	SiMath::Vector minDist(c,HUGE_VALF);
	for ( unsigned int i=1; i<nCompounds; ++i ){
		if ( i >= c ) // when the number of data points asked is greater than the complete data set
			break;
		
		d = 0.0;
		dFar = 0.0;
		index = 0;
		for ( unsigned int k=0; k<c; ++k )
		{
			if ( seen[k] ) continue;
			
			d = SiMath::similarityByEuclidean(n,x[k], mx.getArrayPointer());
			// store the closest distance to the current set
			if ( d < minDist[k] )
				minDist[k] = d;
			
			// check if this distance is farther away
			if ( minDist[k] > dFar ){
				dFar = minDist[k];
				index = k;
			}
		}
		// store current
		mx = x.getRow(index);
		std::cout << index << ": " << mx;
		seen[index] = true;		
	}
	
	// done
	exit(0);
}



void
splitLine(const std::string & line, std::vector<double> & v)
{
	const std::string::size_type len = line.length();
	std::string::size_type i = 0;
	std::string delimiters(" \t\n\r\f");
	std::string sub="";
	
	while ( i < len )
	{
		// eat leading whitespace
		i = line.find_first_not_of(delimiters, i);
		if (i == std::string::npos)
		{
			return;   // nothing left but white space
		}
				
		// find the end of the token
		std::string::size_type j = line.find_first_of(delimiters, i);
				
		// push token onto container
		if (j == std::string::npos) // end of string
		{
			sub = line.substr(i);
			v.push_back(strtod(sub.c_str(),NULL));
			return;
		} 
		else // 
		{
			sub = line.substr(i,j-i);
			v.push_back(strtod(sub.c_str(),NULL));
		}
		// move positions to next substring
		i = j + 1;
	}
}

