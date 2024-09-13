/*
 *  KMeans.h
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

#ifndef  __SIMATH_KMEANS__
#define  __SIMATH_KMEANS__

#include "Definitions.h"
#include "Matrix.h"
#include "Vector.h"

namespace SiMath
{
	
	/**
		\class KMeansError
		\brief  Exception class
	 
		This exception class is derived from the std::runtime_error class and is used to 
		indicate that something went wrong when running the k-Means algroithm.
	 \ingroup _error_handling_
	 */	
	class SIMATH_EXCEPTION KMeansError : public SiMathError
	{
		public:
			/**
				\brief Main error message call
				\param e Error message
			 */
			KMeansError(const std::string & e) : SiMathError("KMeansError : " + e ) {};
	};


	/**
		\struct KMeansParameters
		\brief Definition of the parameters used in the k-Means algorithm.
		\ingroup _model_parameters_
	 
		Usage:
		\code
		SiMath::KMeansParameters p;
		p.dimension = A.nbrColumns();
		p.nbrClusters = 20;
		p.nbrIterations = 200;
		\endcode
	 
		*/
	struct KMeansParameters {
		unsigned int dimension;      ///< Dimension of the data space
		unsigned int nbrClusters;    ///< Maximal number of clusters to be computed
		unsigned int maxIterations;  ///< Maximal number of iterations
		double changeLevel;          ///< Convergence level
		
		/**
			\brief Default parameters:
			- dimension = 0, 
			- nbrClusters = 0, 
			- maxIterations = 1000, 
			- changeLevel = 0.0
		 */
		KMeansParameters() : dimension(0), nbrClusters(0), maxIterations(1000), changeLevel(0.0) {};
	};
	
	/**
		\class KMeans
		\brief Implementation of the k-Means algorithm for clustering
		\ingroup _clustering_
	 
		k-Means clustering is a classic partitioning algorithm. The algorithm starts with
		define k cluster centers. Next each point in the data set is assigned to the closest 
		cluster center. When all points have been assigned to a cluster the centers are recalculated
		as the mean profile of all points in the cluster. Then, the assignment procedure is repeated.
		The algorithm stops when in an iteration none of the points changes cluster.  

	 */
	class KMeans
	{
		public:
			
			/**
				\brief Constructor
				\param p Reference to a predefined set of parameters
			 */
			KMeans(const KMeansParameters & p);

			/**
				\brief Destructor
			 */
			~KMeans();
			
			/**
				\brief Main clustering step
			 
				\param data Reference to a Matrix that holds the data points
				\return std::vector with the cluster labels of all data points
			 */
			std::vector<unsigned int> cluster(const Matrix & data);

			/**
				\brief Clustering step that starts from the currently defined centers
			 
				\param data Reference to a Matrix that holds the data points
				\return A std::vector with the cluster labels of all data points
			 */
			std::vector<unsigned int> update(const Matrix & data);
			
			/**
				\brief Compute the cluster labels of new data
				
				\param data Reference to a Matrix that holds the data points
				\return std::vector with the cluster labels of all data points
			 */
			std::vector<unsigned int> assign(const Matrix & data);
			
			/**
				\brief Compute the cluster label of a new data point
			 
				\param data Reference to a Vector that holds the data point
				\return The cluster label of the data point
			 */
			unsigned int assign(const Vector & data);
			
			/**
				\brief Check if the algorithm is converged
				\return boolean
			 */
			bool isConverged() { return _isConverged; };
			
			/**
				\brief Get the current number of iterations 
			 */
			unsigned int getNbrOfIterations() { return _maxIterations; };

			/**
				\brief Get the current level of required change before convergence 
			 */
			double getChangeLevel() { return _changeLevel; };
			 
			/**
				\brief Get the current number of clusters
			 */
			unsigned int getNbrOfClusters() { return _nbrClusters; };
			
			/**
				\brief Get the center 
			 */
			Vector getCenter(const unsigned int i);

			/**
				\brief Define the center of a cluster
				\param i Cluster number
				\param v Vector with cluster center coordinates
			 */
			void setCenter(unsigned int i, const Vector & v);

			/**
				\brief Set the level of convergence required.
			 
			  This parameter is used to change the level of convergence required. This level is 
				defined by the fraction of points that have been assigned to a different 
				cluster in consecutive iterations. When this fraction is below the preset value
				the algorithm is considered to be converged. It is particulary useful 
				when clustering very large datasets, where full convergence can take some time.
			 
				By default, the level is set to 0 which means full convergence and no point is
				allowed to change cluster. 
			 
				\param d Value between 0 and 1
			 */
			void setChangeLevel(double d);
			
			/**
				\brief Set the number of clusters.
				
				This functions sets the expected number of clusters and 
			 also resets the current centers to zero.
			 
			 \param n Number of clusters
			 */
			void setNbrOfClusters(unsigned int n);
			
			/**
				\brief Set the maximum number of iterations of the clustering
				\param n Maximal number of iterations 
			 */
			void setMaxIterations(unsigned int n) { _maxIterations=n;};
			

		private:
				
			/**
				\struct Center
				\brief Struct to hold the cluster center representation
			 */
			struct Center {
				Vector coord;             ///< Coordinates of the cluster center
				unsigned int members;     ///< Number of members assigned to cluster
				
				Center() : coord(0), members(0) {};
			};
			
			std::vector<Center> _centers;   ///< Local cluster center storage
			unsigned int _nbrClusters;      ///< Total number of clusters
			unsigned int _nbrVar;           ///< Dimension of the data space
			unsigned int _maxIterations;    ///< Maximum number of iterations
			double _changeLevel;            ///< Minimal level of 
			bool _isConverged;              ///< Status variable to indicate if teh training hs converged 
			bool _dataSeen;                 ///< Status variable to indicate if data have been presented to the model
			
			void _pruneCenters();			      ///< Local helper function to prune empty centers. Can change the numbe rof clusters
	};
	
	
}; // end of namespace declaration


#endif   __SIMATH_KMEANS__
