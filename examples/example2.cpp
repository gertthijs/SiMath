/*
 *  example1.cpp
 *  SiMath
 *
 *  Created by Gert Thijs on 10/11/06.
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

	// open tab-delimited text file for reading input data
	// in this example the first column contains the target values
	std::string inFile(argv[1]);
	unsigned int clusters = strtol(argv[2],NULL,10);
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

	std::cerr << "read " << values.size() << "  values\n";
	
	// c is the number of lines (==rows) in the data matrix
	// number of columns is values.size / c
	
	if ( (values.size() % c) != 0 )
	{
		std::cerr << "ERROR: Number of data points is not consistent." << std::endl;
		exit(101);
	}

	unsigned int n = values.size() / c;
	std::cerr << "Number of rows = " << c << std::endl;
	std::cerr << "Number of columns = " << n << std::endl;
	
	SiMath::Vector y(c);     //
	SiMath::Matrix x(c,n-1); // 
	
	// create data matrix by rows from the input values
	for ( unsigned int i=0; i<c; ++i )
	{
		y[i] = values[i*n];
		for ( unsigned int j=0; j<n-1; ++j )
			x[i][j] = values[i*n + 1 + j];
	}
	values.clear();
	
	// normalisation
	SiMath::Vector mx = meanOfAllColumns(x);
	std::cerr << "mean = " << mx;
	SiMath::Vector sdx = stDevOfAllColumns(x);
	std::cerr << "stdev = " << sdx;
	
	// destructive normalisation
	columnNormalise(x,mx,sdx);
	
	// do pca analysis => compute svd of x
	SiMath::SVD pca(x);
	
	// eigen vectors
	SiMath::Matrix v(pca.getV());
	
	// singular values
	SiMath::Vector s(pca.getSingularValues());
	// normalise 
	s /= sqrt(c-1);
	
	std::cerr << "singular values: " << s << std::endl;
	
	// transform x according to eigen vectors
	SiMath::Matrix xrot = SiMath::product(x, v);

	// create a k-means clustering
	SiMath::KMeansParameters p;
	p.nbrClusters = clusters;
	p.dimension = 4;
	
	SiMath::KMeans model(p);
	std::vector<unsigned int> labels = model.cluster(xrot);

	std::cerr << "cluster labels: " << std::endl;
	for ( unsigned int k=0; k<labels.size()/3; ++k )
	  std::cerr << " " << labels[k] << " " << labels[k+50] << " " << labels[k+100] << std::endl;
 
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
