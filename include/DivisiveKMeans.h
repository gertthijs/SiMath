/*
 *  DivisiveKMeans.h
 *
 *  Created by Gert Thijs on 10/05/06.
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


#ifndef  __SIMATH_DIVISIVE_KMEANS__
#define  __SIMATH_DIVISIVE_KMEANS__

#include "Definitions.h"
#include "Matrix.h"
#include "Vector.h"

#include <queue>
#include <map>

namespace SiMath
{
	/**
		\class  DivisiveKMeansError
		\brief  Exception class
	 
		This exception class is derived from the std::runtime_error class and is used to 
		indicate that something went wrong during the divisive k-means algorithm.
	 \ingroup _error_handling_
	 */	
	class SIMATH_EXCEPTION DivisiveKMeansError : public SiMathError
	{
		public:
			DivisiveKMeansError(const std::string & e) : SiMathError("DivisiveKMeansError : " + e ) {};
	};


	/** 
		\struct DivisiveKMeansParameters
		\brief Struct to hold all parameters for the divisive k-means algorithm
		\ingroup _model_parameters_

		Usage:
			\code
			SiMath::DivisiveKMeansParameters p;
			p.dimension = A.nbrColumns();
			p.maxMembers = 20;
			\endcode
		*/
	struct DivisiveKMeansParameters {
		unsigned int dimension;      ///< Dimension 
		unsigned int maxDepth;       ///< Maximal depth of the tree to build
		unsigned int maxIterations;  ///< Maximal number of iteration in each k-means step
		unsigned int maxMembers;     ///< Maximal number of members in a cluster

		/**
			\brief Default parameters:
			 - dimension = 0, 
		   - maxDepth = 8, 
			 - maxIterations = 1000, 
		   - maxMembers = 1
		 */
		DivisiveKMeansParameters() : dimension(0), maxDepth(8), maxIterations(1000), maxMembers(1) {};
	};

	/**
		\class DivisiveKMeans
		\brief Class to create a partioning of the data using a divisive k-means procedure
	 \ingroup _clustering_
	 
		The divisive k-means algorithm starts by dividing the data set in two distinct cluster
		using a winner-takes-all principle as in k-means. In an iterative procedure will each
	  cluster be divided again in two distinct clusters if it has more elements than the
		predefined maximum. When applying this method, a binary tree is generated. 
		Since the number of nodes in the tree grows exponential the user can impose additional
		treshold on the depth of the tree. 
	 
		\code
		DivisiveKMeansParameters p;
		p.dimension = A.nbrColumns();
		p.depth = 10;
		p.maxMembers = 30;
		
		DivisiveKMeans model(p);
		Vector labels = model.cluster(A);
		
		// only print out the final nodes
		for ( unsigned int n=1; n<=model.getNbrOfClusters(); ++n )
		{
			if ( model.getChild1(n) == 0 && model.getChild2(n) == 0 )
			{
				Vector v = model.getCenter(n);
				...
			}
		}
		\endcode
		
	 */
	class DivisiveKMeans
	{

		public:
			/**
				\brief Constructor from parameter set
			  
				\param p Set of predefined parameters.
			 */
			DivisiveKMeans(const DivisiveKMeansParameters & p);

			/**
				\brief Destructor
			 */
			~DivisiveKMeans();
	
			/**
				\brief Main procedure to generate the partitioning of the data
				\param data Input data
				\return std::vector with cluster labels for all data points
			 */
			std::vector<unsigned int> cluster(const Matrix & data);

			/**
				\brief Compute the cluster labels of new data
			 
			  \param data Reference to a Matrix that holds the data points
			  \return std::vector with the cluster labels of all data points
			 */
			std::vector<unsigned int> assign(const Matrix & data);

			/**
				\brief Compute the cluster labels of one single data vector
			 
			 \param data Reference to a Matrix that holds the data points
			 \return The cluster label
			 */
			unsigned int assign(const Vector & data);
	
			/**
				\brief Check if the algorithm is converged
			 */
			bool isConverged() { return _isConverged; };
			
			/**
				\brief Get the number of iterations
			 */
			unsigned int getNbrOfIterations() { return _maxIterations; };
			
			/**
				\brief Get the number of clusters found
			 */
			unsigned int getNbrOfClusters() { return _nbrClusters; };
			
			/**
				\brief Get the current maximal depth
			 */
			unsigned int getMaxDepth() { return _maxDepth;};
			
			/**
				\brief Get the current maximum number of members in a cluster
			 */
			unsigned int getMaxNbrOfmembers() { return _maxMembers;};
			
			/**
				\brief Get the coordinates of the i-th center
				\param i Cluster index
			 */
			Vector getCenter(const unsigned int i);

			/**
				\brief Get the cluster number of the parent cluster 
				\param i Cluster index
			 */
			unsigned int getParent(const unsigned int i);

			/**
				\brief Get the cluster number of the first child cluster 
				\param i Cluster index
				\return Cluster index. If 0 then there is no first child.
			 */
			unsigned int getChild1(const unsigned int i);
			
			/**
				\brief Get the cluster number of the second child cluster 
				\param i Cluster index
				\return Cluster index. If 0 then there is no second child.
			 */
			unsigned int getChild2(const unsigned int i);
			
			/**
				\brief Clear the complete internal state of the object
			 */
			void reset();
	
			/**
				\brief Add a cluster center to build a model 
				
				This function is made to upload a predefined clustering center by center.
				To do this parent nodes should always be added before child nodes. This
				implies that the cluster index i should always be greater than the parent number n.
				
				\param i Cluster index
				\param v Vector of coordinates
				\param n Parent number
			 */
			void setCenter(unsigned int i, const Vector & v, unsigned int n);
			
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
				\brief Define number of clusters (only needed when reconstructing a tree)
				\param n Number of clusters
				*/
			void setNbrOfClusters(unsigned int n);

			/**
				\brief Update the maximum number of iterations of each k-means run
				\param n Maximum number of iterations
			 */
			void setMaxIterations(unsigned int n) { _maxIterations=n;};

			/**
				\brief Set maximum number of members in a cluster
				\param n Maximum number
			 */
			void setMaxNbrOfMembers(unsigned int n) { _maxMembers=n;};
		
			/**
				\brief Update maximal depth of the tree
				\param n Maximal depth
			 */
			void setMaxDepth(unsigned int n) { _maxDepth=n;};
			

		private:
			/**
				\struct Center
				\brief Structure to hold the centers of the clusters
			*/
			struct Center {
				Vector coord;          ///< Coordinates of the center
				unsigned int number;   ///< Cluster number
				unsigned int depth;    ///< Depth at which the node is in the tree
				unsigned int members;  ///< Number of members assigned to node during clustering
				unsigned int parent;   ///< Cluster number of the parent node
				unsigned int child1;   ///< Cluster number of the first child node
				unsigned int  child2;  ///< Cluster number of the seconf child node
				
				Center() : coord(0), number(0), depth(0), members(0), parent(0), child1(0), child2(0) {};
			};
			
			
			std::vector<Center> _centers;                  ///< Local vcetor to hold all centers
			std::queue<Center> _processQueue;              ///< Local queue to hold all clusters that can be further splitted
			std::map<unsigned int, Center *> _centerMap;   ///< Local mapping of cluster number and actual center
			
			unsigned int _nbrClusters;           ///< Number of clusters found      
			unsigned int _nbrVar;                ///< Dimension of the data matrix
			unsigned int _maxMembers;            ///< Maximal number of members in a cluster 
			unsigned int _maxIterations;         ///< Maximal number of iterations
			unsigned int _maxDepth;              ///< Maximal depth of the clustering
			unsigned _latestNumber;              ///< Latest generated cluster number
			
			bool _isConverged;                   
			bool _dataSeen;
			
			void _splitCluster(const Matrix & data, std::vector<unsigned int> & labels, Center & c);
			void _createMap();
			
	};


}; // end of namespace declaration


#endif   __SIMATH_KMEANS__
