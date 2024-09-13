/*
 *  SOM.h
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

#ifndef __SIMATH_SOM__
#define __SIMATH_SOM__

// #include "SiMath.h"
#include "Vector.h"
#include "Matrix.h"

namespace SiMath
{

	/**
		\class SOMError
		\brief  Exception class
	 
		This exception class is derived from the std::runtime_error class and is used to 
		indicate that something went wrong during the SOM algorithm.
	 \ingroup _error_handling_
	 */	
	class SIMATH_EXCEPTION SOMError : public SiMathError
	{
		public:
			SOMError(const std::string & e) : SiMathError("SOMError : " + e ) {};
	};
	
	/**
		\struct SomParameters
		\brief  Set of necessary SOM parameters
	 \ingroup _model_parameters_

	 Usage:
	 \code
	 SiMath::SOMParameters p;
	 p.dimension = A.nbrColumns();
	 p.xGridSize = 10;
	 p.yGridSize = 10;
	 p.neighbourLevel = 4;
	 p.maxIterations = 5 * A.nbrRows(); // on average each datapoint is presented 5 times to the SOM 
	 \endcode
	 
	 */
	struct SomParameters
	{
		unsigned int dimension;         ///< Dimension of the data space
		unsigned int xGridSize;         ///< Dimension of the grid in x direction
		unsigned int yGridSize;         ///< Dimension of the grid in y direction
		unsigned int maxIterations;     ///< Total number of iterations
		unsigned int maxStepSize;       ///< Maximal step size to iterate over the data
		unsigned int neighbourLevel;    ///< NUmber of neighbouring cells to be updated in the first iterations
		double tau;                     ///< Learning rate
		
		/**
			\brief Default parameters:
				- dimension = 0,
				- xGridSize = 0 and yGridSize = 0, 
				- maxIterations = 9999, 
				- maxStepSize = 999, 
				- neighbourLevel = 2, 
				- tau =0.9
		 */			
		SomParameters() : dimension(0), xGridSize(0), yGridSize(0), maxIterations(9999), maxStepSize(999), neighbourLevel(2), tau(0.9) {};
	};
	
	/**
		\struct SomGrid
		\brief Structure to hold the grid coordinates of a data point and the distance to the grid reference
	 */
	struct SomGrid
	{
		unsigned int x;   ///< x coordinate
		unsigned int y;   ///< y coordinate
		double       d;   ///< distance to grid reference
		/**
			\brief Default coordinates are (0,0)
		 */
		SomGrid() : x(0), y(0), d(0.0) {};
	};
	
	
	/**
		\class SOM
		\brief Class to construct, train and manipulate a Self-organising Map
		\ingroup _clustering_
	 
		This SOM implementation maps the data onto a toroidal grid. This means that in a X-by-Y grid
		the cells at one end (eg. [0,2]) are neigbours of the cells at the other end ([N,2] in this
		example).

	  The algorithm runs for a predefined number of iterations and in each iteration one
	  data point is presented to the map. The number of neighbouring cells that is updated
		depends on the iteration. In the beginning the maximum number is updated and near the end
		only one cell is updated. The update of a cell also depends on the learning rate tau which
	  decreases linearly over all iterations between the start value and 0.
	 
	  The data is presented to the algorithm in a semi-randomised way. In each iteration a 
		a new data point is sampled in the interval [1,stepsize) from the current selected point.
		When at the end of the data matrix, sampling proceeds at the beginning.
	 
	 */
	class SOM
	{
		public:
			/**
				\brief Constructor based on predefined set of parameters
				\param mp SomParameters
			 */
			SOM(const SomParameters & mp);
		
			/**
				\brief Destructor
			 */
			~SOM();


			/**
				\brief Define the number of iterations
				\param i Number of iterations
			 */
			void setMaxIterations(unsigned int i);
			
			/**
				\brief Define the maximum stepsize to iterate over the data
				\param n Stepsize 
			 */
			void setMaxStepSize(unsigned int n);
			
			/**
				\brief Define the initial learning rate of the cell updating
				\param t Learning rate 
			 */
			void setTau(double t);
			
			/**
				\brief Define how many neighbours should be updated initially
					
				The neighbourhood level defines which cells around a selected will be updated. In the 
				following example the neighbourhood in a 5x5 toroidal grid is illustrated for two different 
				cells. All cells numbered 1 are at level 1 and the ones number 2 are at level 2 from the cell 
				with number 0.
							 
			 \code
			 [2 2 2 2 2     [1 1 1 2 2
			  2 1 1 1 2      1 0 1 2 2
			  2 1 0 1 2      1 1 1 2 2
			  2 1 1 1 2      2 2 2 2 2
			  2 2 2 2 2]     2 2 2 2 2]
			 \endcode
				
			 \param n Neighbourhood level
			 */
			void setNeighbourLevel(unsigned int n);
			
			/**
				\brief Define the coordinates of cell [x,y]
				\param x x grid index
				\param y y grid index
				\param v Vector of coordinates of cell representative 
			 */
			void setCenter(const unsigned int x,const unsigned int y, Vector & v);
			
			/**
				\brief Main procedure to create the SOM from a data matrix
				\param data Data matrix
				\param clusterData boolean to indicate if the data should be clustered and stored in the labels
				\return Vector of grid coordinates for every data point
			 */
			std::vector<SomGrid> cluster(const SiMath::Matrix & data, bool clusterData=true);

			/**
				\brief Procedure to update a SOM from a new data matrix
				\param data Data matrix
			 \param clusterData boolean to indicate if the data should be clustered and labels should be returned.
			 \return Vector of grid coordinates for every data point
			 */
			std::vector<SomGrid> update(const SiMath::Matrix & data, bool clusterData=true);
			
			/**
				\brief Method to map data onto the SOM
				\brief data Data matrix to map
				\return Vector of grid coordinates for every data point
				*/
			std::vector<SomGrid> assign(const SiMath::Matrix & data);
			
			/** 
				\brief Method to map a single data Vector onto the SOM
				\param data Single data vector
				\return Grid coordinates of the data point
			 */
			SomGrid assign(const Vector & data);
			
			// unsigned int getNbrVariables() const { return _nbrVar; };
			
			/**
				\brief Get the current number of cells in the x direction
			 */
			unsigned int getXGridSize() const { return _xGridSize; };

			/**
				\brief Get the current number of cells in the y direction
			 */
			unsigned int getYGridSize() const { return _yGridSize; };

			/**
				\brief Get current number of iterations
			 */
			unsigned int getMaxIterations() const { return _maxIterations; };

			/**
				\brief Get current step size
			 */
			unsigned int getMaxStepSize() const { return _maxStep;};
			
			/**
				\brief Get the representative of cell [x,y]
				\param x x grid index
				\param y y grid index
			 */
			Vector getCenter(const unsigned int x,const unsigned int y);
			
		private:
			unsigned int _nbrVar;            ///< Dimension of the data space
			unsigned int _xGridSize;         ///< Number of grid points in the x-direction
			unsigned int _yGridSize;         ///< Number of grid points in the y-direction
			unsigned int _maxIterations;     ///< Maximal number of iterations
			unsigned int _maxStep;           ///< Maximal setp size to iterate over the data
			unsigned int _neighbours;        ///< Maximal number of neighbours to update
			double _tau;                     ///< Learning rate
			
			/**
				\brief Internal cell representation
			 */
			struct Center {
				Vector coord;          ///< Coordinates  
				unsigned int members;  ///< Number of data points assigned to this cell during training
				unsigned int x;        ///< x grid index
				unsigned int y;        ///< y grid index
				
				Center(): coord(0), members(0), x(0), y(0) {};
			};

			std::vector<Center> _somCenters;  ///< Local storage of cell representatives
			
	};
	
};

#endif __SIMATH_SOM__
