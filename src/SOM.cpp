/*
 *  SOM.cpp
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

#include "SOM.h"
#include "Similarity.h"
#include "Matrix.h"
#include "Vector.h"

SiMath::SOM::SOM(const SomParameters & p) :
	_nbrVar(p.dimension),
	_xGridSize(p.xGridSize),
	_yGridSize(p.yGridSize),
	_maxStep(p.maxStepSize),
	_maxIterations(p.maxIterations),
	_neighbours(p.neighbourLevel),
	_tau(p.tau),
	_somCenters(p.xGridSize * p.yGridSize)
{
	unsigned int index(0);
	for ( unsigned int i=0; i<_xGridSize; ++i)
	{
		for ( unsigned int j=0; j<_yGridSize; ++j)
		{
			_somCenters[index].coord.reset(_nbrVar);
			_somCenters[index].members = 0;
			_somCenters[index].x = i;
			_somCenters[index].y = j;			
			++index;
		}
	}				
}


SiMath::SOM::~SOM()
{
	
}


//-------------------------------------------------------------------------------------------------------
// TRAINING
std::vector<SiMath::SomGrid>
SiMath::SOM::cluster(const SiMath::Matrix & data, bool clusterData)
{
	unsigned int n = data.nbrRows();
	unsigned int m = data.nbrColumns();

	// create empty center matrix ==> add two columns to store grid indices
	if ( _nbrVar != m )
		throw(SOMError("Dimension of the data do not correspond with the dimension of the model."));
		
	// grid indices
	unsigned int xn = 0;
	unsigned int yn = 0;
	
	std::cerr << "Start SOM procedure." << std::endl;
	
	// add spectrophore to list
	unsigned int index(0);
	for ( xn=0; xn<_xGridSize; ++xn )
	{
		for ( yn=0; yn<_yGridSize; ++yn )
		{
			// set the data block of the center
			double s = 0.0;
			for ( unsigned int j=0; j<m; ++j)
			{	
				double d = -1.0 + 2.0 * ((double)rand()/RAND_MAX);
				_somCenters[index].coord[j] = d;
				s += (d*d);
			}
			// normalize
			s = sqrt(s/m);
			for ( unsigned int j=0; j<m; ++j)
				_somCenters[index].coord[j] /= s;
			
			
			// augment index
			++index;
		}
	}
	// std::cerr << "Initial centers:\n" << _somCenters;
	
	
	// iteration parameters
	unsigned int levelSteps = _maxIterations/(_neighbours+1);
	int level = _neighbours;
	unsigned int dataIndex = 0;	
	
	// create a boolean matrix to check if a cell has been visited in an iteration
	std::vector<bool> cellSeen(_xGridSize * _yGridSize,0);
	
	// linear decreasing learning rate
	double tauStep = _tau/_maxIterations;
	double tau = _tau + tauStep;
	
	std::cerr << "  initial learning rate: " << tau << std::endl;
	
	for ( unsigned int nIter = 0; nIter<_maxIterations; ++nIter )
	{
		// move to next data point
		unsigned int step = 1 + (unsigned int)(((double)rand()/RAND_MAX) * _maxStep);
		dataIndex += step;
		while ( dataIndex >= n )
			dataIndex -= n;

		// make a copy of the current row
	  Vector currentRow(data.getRow(dataIndex));

		// update the learning parameter tau
		tau -= tauStep;
		
		// find cell closest to data point (si)
		unsigned int xBest = 0;
    unsigned int yBest = 0;
		unsigned int index = 0;
		double closest = HUGE_VALF;
		for ( xn=0; xn<_xGridSize; ++xn )
		{ 
			for ( yn=0; yn<_yGridSize; ++yn )
			{ 
				index = (xn * _yGridSize) + yn;
				double d = SiMath::similarityByEuclidean(_somCenters[index].coord,currentRow);
				if (d < closest)
				{ 
					xBest = xn;
					yBest = yn;
					closest = d;
				}
				cellSeen[index] = 0;
			}
		}     
		
		// update cells
		if ( (nIter % levelSteps) == 0 )  
		{
			level = _neighbours - nIter/levelSteps;
			if ( level == -1 )
				break;
			std::cerr << "Setting level of neighbours at iteration " << nIter << " => " << level << std::endl;
			std::cerr << "    learning constant: " << tau << std::endl;
		}
		// std::cerr << "level = " << level << " (" << xBest << "," << yBest << "):"<< std::endl;
		for ( int i=-level; i<=level; ++i )
		{
			int x = xBest + i;
			if ( x < 0 )
			{
				xn = x + _xGridSize;
			}
			else if (  x >= _xGridSize )
			{
				xn = x - _xGridSize;
			}
			else
			{	
				xn = x;
			}
			
			for ( int j=-level; j<=level; ++j )
			{
				int y = yBest + j;
				if ( y < 0 )
				{
					yn = y + _yGridSize;
				}
				else if ( y >= _yGridSize )
				{
					yn = y - _yGridSize;
				}
				else
				{
					yn = y;
				}
				
				// get cell index from 
				index = (xn * _yGridSize) + yn;
				
				if ( cellSeen[index] == 1 )
					continue;
				
				cellSeen[index] = 1;
				
				
				// update cell center
				for ( unsigned int j=0; j<m; j++ )
				{ 
					_somCenters[index].coord[j] += tau * (currentRow[j] - _somCenters[index].coord[j]);
				}
				
				// normalize centers
				double s = 0.0;
				for ( unsigned int j=0; j<m; j++ )
				{ 
					double d = _somCenters[index].coord[j];
					s += d*d;
				}
				if ( s > 0 )
				{ 
					s = sqrt(s/m);
					for ( unsigned int j=0; j<m; j++ )
						_somCenters[index].coord[j] /= s;
				}
      }
    }
	}
	
	// create a matrix to store the cluster indices for each data point
	std::vector<SomGrid>  cl;	
	if ( clusterData )
	{
		for ( unsigned int i=0; i<n; ++i )
		{
			Vector v(data.getRow(i));
			
			// compute mapping to cluster
			unsigned int xBest = 0;
			unsigned int yBest = 0;
			double closest(HUGE_VALF);
			
			// find cell closest to data point (si)
			for ( xn=0; xn<_xGridSize; ++xn )
			{ 
				for ( yn=0; yn<_yGridSize; ++yn )
				{ 
					unsigned int index = (xn * _yGridSize) + yn; 
					double d = SiMath::similarityByEuclidean(_somCenters[index].coord, v);
					if (d < closest)
					{	
						xBest = xn;
						yBest = yn;
						closest = d;
					}
				}
			}
			SomGrid sg;
			sg.x = xBest;
			sg.y = yBest;
			sg.d = closest;
			cl.push_back(sg);
			// cl[i][0] = xBest;
			// cl[i][1] = yBest;
		}
	}
	std::cerr << "done SOM clustering." << std::endl;
	return cl;
}


std::vector<SiMath::SomGrid>
SiMath::SOM::update(const SiMath::Matrix & data, bool clusterData)
{
	unsigned int n = data.nbrRows();
	unsigned int m = data.nbrColumns();
	
	if ( _nbrVar != m )
		throw(SOMError("New data not compatible with current set of cluster centers."));
	
	// grid indices
	unsigned int xn = 0;
	unsigned int yn = 0;
	
	// iteration parameters
	unsigned int levelSteps = _maxIterations/(_neighbours+1);
	int level = _neighbours;
	unsigned int dataIndex = 0;	

	// create a boolean matrix to check if a cell has been visited in an iteration
	std::vector<bool> cellSeen(_xGridSize * _yGridSize,0);
	
	double tauStep = _tau/_maxIterations;
	double tau = _tau + tauStep;
	
	std::cerr << "  initial learning rate: " << tau << std::endl;
	
	for ( unsigned int nIter = 0; nIter<_maxIterations; ++nIter )
	{
		// move to next data point
		unsigned int step = 1 + (unsigned int)(((double)rand()/RAND_MAX) * _maxStep);
		dataIndex += step;
		if ( dataIndex >= n )
			dataIndex -= n;
		Vector currentRow(data.getRow(dataIndex));

		// update the learning parameter tau
		tau -= tauStep;
		
		// find cell closest to data point (si)
		unsigned int xBest = 0;
    unsigned int yBest = 0;
		unsigned int index = 0;
		double closest = HUGE_VALF;
		for ( xn=0; xn<_xGridSize; ++xn )
		{ 
			for ( yn=0; yn<_yGridSize; ++yn )
			{ 
				// index = (xn * _yGridSize) + yn; 
				double d = SiMath::similarityByEuclidean(_somCenters[index].coord, currentRow);
				if (d < closest)
				{ 
					xBest = xn;
					yBest = yn;
					closest = d;
				}
				cellSeen[index] = 0;
				++index;
			}
		}     
		
		// update cells
		if ( (nIter % levelSteps) == 0 )  
		{
			level = _neighbours - nIter/levelSteps;
			if ( level == -1 )
				break;
			std::cerr << "Setting level of neighbours at iteration " << nIter << " => " << level << std::endl;
			std::cerr << "    learning constant: " << tau << std::endl;
		}
		// std::cerr << "level = " << level << " (" << xBest << "," << yBest << "):"<< std::endl;
		for ( int i=-level; i<=level; ++i )
		{
			int x = xBest + i;
			if ( x < 0 )
			{
				xn = x + _xGridSize;
			}
			else if (  x >= _xGridSize )
			{
				xn = x - _xGridSize;
			}
			else
			{	
				xn = x;
			}
			
			for ( int j=-level; j<=level; ++j )
			{
				int y = yBest + j;
				if ( y < 0 )
				{
					yn = y + _yGridSize;
				}
				else if ( y >= _yGridSize )
				{
					yn = y - _yGridSize;
				}
				else
				{
					yn = y;
				}
				
				// get cell index from 
				index = (xn * _yGridSize) + yn;
				
				if ( cellSeen[index] == 1 )
					continue;
				
				cellSeen[index] = 1;
				
				// update cell center
				for ( unsigned int j=0; j<m; j++ )
				{ 
					_somCenters[index].coord[j] += tau * (currentRow[j] - _somCenters[index].coord[j]);
				}
				
				// normalize centers
				double s = 0.0;
				for ( unsigned int j=0; j<m; j++ )
				{ 
					double d = _somCenters[index].coord[j];
					s += d*d;
				}
				if ( s > 0 )
				{ 
					s = sqrt(s/m);
					for ( unsigned int j=0; j<m; j++ )
						_somCenters[index].coord[j] /= s;
				}
      }
    }
	}
	
	// create a matrix to store the cluster indices for each data point
	std::vector<SiMath::SomGrid> cl;	
	if ( clusterData ){
		std::cerr << " clustering the data ..." << std::endl;
		for ( unsigned int i=0; i<n; ++i )
		{	
			Vector currentRow(data.getRow(i));
			
			// compute mapping to cluster
			unsigned int xBest = 0;
			unsigned int yBest = 0;
			double closest(HUGE_VALF);
			// find cell closest to data point (si)
			unsigned int index(0);
			for ( xn=0; xn<_xGridSize; ++xn )
			{ 
				for ( yn=0; yn<_yGridSize; ++yn )
				{ 
					// unsigned int index = (xn * _yGridSize) + yn; 
					double d = SiMath::similarityByEuclidean(_somCenters[index].coord, currentRow);
					if (d < closest)
					{ 
						xBest = xn;
						yBest = yn;
						closest = d;
					}
					++index;
				}
			}
			SomGrid sg;
			sg.x = xBest;
			sg.y = yBest;
			sg.d = closest;
			cl.push_back(sg);
		}
	}
	
	return cl;
}


std::vector<SiMath::SomGrid>
SiMath::SOM::assign(const SiMath::Matrix & data )
{
	// check dimensions of the centers and the data
	if ( _nbrVar != data.nbrColumns() )
		throw(SOMError("New data not compatible with current set of cluster centers."));

	// create a matrix to store the cluster indices for each data point
	std::vector<SiMath::SomGrid> cl;	
	for ( unsigned int i=0; i<data.nbrRows(); ++i )
	{	
		Vector currentRow(data.getRow(i));

		// compute mapping to cluster
		unsigned int xBest = 0;
		unsigned int yBest = 0;
		double closest(HUGE_VALF);
		// find cell closest to data point (si)
		unsigned int index(0);
		for ( unsigned int xn=0; xn<_xGridSize; ++xn )
		{ 
			for ( unsigned int yn=0; yn<_yGridSize; ++yn )
			{ 
				// unsigned int index = (xn * _yGridSize) + yn; 
				double d = SiMath::similarityByEuclidean(_somCenters[index].coord, currentRow);
				if (d < closest)
				{ 
					xBest = xn;
					yBest = yn;
					closest = d;
				}
				++index;
			}
		}
		SomGrid sg;
		sg.x = xBest;
		sg.y = yBest;
		sg.d = closest;
		cl.push_back(sg);
	}

	return cl;
}


SiMath::SomGrid 
SiMath::SOM::assign(const SiMath::Vector & data )
{
	// check dimensions of the centers and the data
	if ( _nbrVar != data.size() )
		throw(SOMError("New data not compatible with current set of cluster centers."));
	
	// create a matrix to store the cluster indices for each data point
	SiMath::SomGrid sg;

	// compute mapping to cluster
	unsigned int xBest = 0;
	unsigned int yBest = 0;
	double closest(HUGE_VALF);
	// find cell closest to data point 
	unsigned int index=0;
	for ( unsigned int xn=0; xn<_xGridSize; ++xn )
	{ 
		for ( unsigned int yn=0; yn<_yGridSize; ++yn )
		{ 
			// unsigned int index = (xn * _yGridSize) + yn; 
			double d = SiMath::similarityByEuclidean(_somCenters[index].coord, data);
			if (d < closest)
			{ 
				xBest = xn;
				yBest = yn;
				closest = d;
			}
			++index;
		}
	}
	
	// define cluster labels
	sg.x = xBest;
	sg.y = yBest;
	sg.d = closest;
	
	return sg;
}


SiMath::Vector 
SiMath::SOM::getCenter(const unsigned int x,const unsigned int y)
{
	if ( x<0 || x >= _xGridSize )
		throw(SOMError("x index out of range."));
	if ( y<0 || y >= _yGridSize )
		throw(SOMError("y index out of range."));
	
	return _somCenters[y + x*_yGridSize].coord;
}


void 
SiMath::SOM::setMaxIterations(unsigned int i)
{
	_maxIterations = i;
}

void 
SiMath::SOM::setMaxStepSize(unsigned int n)
{
	_maxStep = n;
}

void 
SiMath::SOM::setTau(double t)
{
	if ( t < 0 || t > 1 )
		throw(SOMError("tau should be between 0 and 1."));
	
	_tau = t;
}

void 
SiMath::SOM::setNeighbourLevel(unsigned int n)
{
	_neighbours = n;
}


void 
SiMath::SOM::setCenter(const unsigned int x,const unsigned int y, Vector & v)
{
	if ( x<0 || x >= _xGridSize )
			throw(SOMError("x index out of range."));
	if ( y<0 || y >= _yGridSize )
		throw(SOMError("y index out of range."));
	if ( v.size() != _nbrVar )
		throw(SOMError("Center vector is not compatible with the number of variables defined."));
	
	unsigned int i = (x * _yGridSize) + y; 
	
	_somCenters[i].coord = v;
	_somCenters[i].x = x;
	_somCenters[i].y = y;
	_somCenters[i].members = 0;
	
	return;
}
